use std::path::PathBuf;
use structopt::StructOpt;
use rayon::prelude::*;
use itertools::{izip, Itertools};

mod measure_distance;
mod simple_math;
mod config;
mod binning;

use simple_math::{cell_volume, sphere_vol};
use config::Config;
use binning::{make_bins, BinResult, sample_file, add_bins};

/// Reads a Ge file and outputs
/// the pair correlation function on
/// stdout in format
/// domain g2 std_err(g2) poisson_std_err(g2) bin_width
/// Warning: program assumes all numeric
/// conversions are trivially valid, don't
/// put in absurdly large numbers
/// Currently assumes all configurations
/// have the same number density
#[derive(StructOpt, Debug)]
#[structopt(name = "pair_correlation")]
pub struct Opt {
   
    // From https://github.com/TeXitoi/structopt
    /// Verbosity of monitoring output
    #[structopt(short, long, parse(from_occurrences))]
    verbosity: u8,

    /// Gives the maximum radius to sample to
    #[structopt(short, long, default_value = "1.0")]
    cutoff: f64,

    /// Gives the number of bins
    #[structopt(short, long, default_value = "1")]
    nbins: usize,

    /// Gives an offset for every bin, i.e. makes a gap around the origin,
    /// to inhibit divide by s_1(r) effect
    #[structopt(short, long, default_value = "0.0")]
    offset: f64,

    /// implements a really basic autoscaling of bin size
    /// based on Poisson assumption, give the initial first
    /// bin width, overrides --nbins
    #[structopt(long)]
    autoscale: Option<f64>,

    /// Give the block sizes to compute for uncertainty analysis
    #[structopt(long)]
    blocks: Option<Vec<usize>>,

    /// Gives the list of files, interpreted as ensemble average
    #[structopt(parse(from_os_str))]
    files: Vec<PathBuf>,

}

// Takes bin count and converts it to an approx
// to g_2 by normalization
// Return sum of variances for block testing if needed
fn format_output(count: &[usize],
                 count2: &[usize],
                 domain: &[f64],
                 lower_limit: &[f64],
                 upper_limit: &[f64],
                 n_particles: usize,
                 dim: usize,
                 n_ens: usize,
                 rho: f64) -> f64
{
    let mut sum_of_variance = 0.0;
    for (c, c2, r, l, u) in izip!(count,
                                  count2,
                                  domain,
                                  lower_limit,
                                  upper_limit) 
    {
        let width = *u - *l;
        // 2*count/(n_particles*n_ens) = rho s1(r) g2(r) dr
        let f_n_ens = n_ens as f64;
        let ens_c = (*c as f64)/f_n_ens;
        // Sample variance reminder:
        // https://en.wikipedia.org/wiki/Variance
        let ens_var_c = (*c2 as f64)/(f_n_ens - 1.0) - ens_c*ens_c*f_n_ens/(f_n_ens - 1.0);
        // Strategy for normalization taken from Ge's code
        let g2_coeff = 2.0
            /((sphere_vol(dim, *u) - sphere_vol(dim, *l))
                *rho*(n_particles as f64)
            );
        let g2 = ens_c*g2_coeff;
        let var_g2 = ens_var_c*g2_coeff*g2_coeff;
        let std_err_g2 = var_g2.sqrt()/((n_ens as f64).sqrt());
        let poisson_est = (*c as f64).sqrt()*g2_coeff/f_n_ens;
        println!("{} {} {} {} {}", r, g2, std_err_g2, poisson_est, width);
        sum_of_variance += std_err_g2.powi(2);
    }
    sum_of_variance
}

// From 
// https://stackoverflow.com/questions/39204908/how-to-check-release-debug-builds-using-cfg-in-rust
// and
// https://vallentin.io/2019/06/06/versioning
#[cfg(debug_assertions)]
const BUILD_MODE: &str = "debug";
#[cfg(not(debug_assertions))]
const BUILD_MODE: &str = "release";

#[cfg(feature = "using_make")]
const USING_MAKE: &str = "true";
#[cfg(not(feature = "using_make"))]
const USING_MAKE: &str = "false";

fn main() {

    eprintln!("Build info for pair_correlation:");
    eprintln!("Using make: {}", USING_MAKE);
    eprintln!("Cargo version: {}", env!("C_VER"));
    eprintln!("Commit SHA: {}", env!("VERGEN_SHA"));
    eprintln!("Commit date: {}", env!("VERGEN_COMMIT_DATE"));
    eprintln!("Version: {}", env!("VERGEN_SEMVER"));
    eprintln!("Build: {}", BUILD_MODE);
    eprintln!("Build time: {}", env!("VERGEN_BUILD_TIMESTAMP"));
    eprintln!("Compiler version: {}", env!("V_RUSTC"));
    eprintln!("RUSTFLAGS: {}", env!("S_RUSTFLAGS"));
    eprintln!("Target: {}", env!("VERGEN_TARGET_TRIPLE"));
    eprintln!("Clean working directory for build: {}", env!("WD_IS_CLEAN"));
    eprintln!("Start program:\n");

    let opt = Opt::from_args();
    if opt.verbosity > 0 {
        eprintln!("{:?}", opt);
    }

    let first_config = Config::parse(&opt.files[0]);

    // Assume that every configuration has the same
    // density as the first
    let rho: f64 = (first_config.n_particles as f64)
                    / (cell_volume(first_config.dim, &first_config.unit_cell));

    // Make the bins
    // Bin values between lower_limit <= val < upper_limit
    // domain is given as midpoint
    let (domain, lower_limit, upper_limit) = make_bins(&opt, &first_config);
    let n_bins = domain.len();
    let n_ens = opt.files.len();

    // Bin each config separately
    let individual_results: Vec<BinResult> = opt.files.par_iter()
            .enumerate()
            .map(|(i, x)| { 
                if opt.verbosity > 1 { eprintln!("Working on {}", i); };
                sample_file(x, &lower_limit, &upper_limit) 
            })
            .collect();

    // Combine each configuration's results through either
    // a block averaging analysis or naive analysis
    // Block averaging and variance analysis
    if let Some(blocks) = opt.blocks.clone() {
        let mut blocked_results: Vec<BinResult> = Vec::new();
        for size in blocks.clone() {
            if opt.verbosity > 0 {
                println!("----- {} block -----", size);
            }
            blocked_results = individual_results.iter()
                    .chunks(size).into_iter() // Split into blocks
                    .map(|x| 
                         x.fold(BinResult {dim: 0, n_particles: 0, count: vec![0; n_bins], count2: vec![0; n_bins] },
                            |acc, x| add_bins(acc, x)))
                    .collect();

            // Original count2 is ignored and blocked value computed
            let final_result: BinResult = blocked_results.iter()
                    .fold(BinResult {dim: 0, n_particles: 0, count: vec![0; n_bins], count2: vec![0; n_bins] },
                        |acc, x| add_bins(acc, x));

            let sum_of_variance = format_output(&final_result.count,
                      &final_result.count2,
                      &domain,
                      &lower_limit,
                      &upper_limit,
                      final_result.n_particles*size,
                      final_result.dim,
                      n_ens/size,
                      rho
            );

            if opt.verbosity > 0 {
                println!("Summary uncertainty estimate (block size, statistic)");
                println!("sue {} {}", size, sum_of_variance);
            }
        }

        if opt.verbosity > 0 {
            println!("----- Drift analysis, using final block size -----");
            for result in blocked_results.iter() {
                format_output(&result.count,
                      &result.count2,
                      &domain,
                      &lower_limit,
                      &upper_limit,
                      result.n_particles,
                      result.dim,
                      *blocks.last().unwrap(),
                      rho
                );
            }
        }

    } else { // Simple averaging and variance analysis, assuming independence
        let summed_result: BinResult = individual_results.iter()
             .fold(BinResult {dim: 0, n_particles: 0, count: vec![0; n_bins], count2: vec![0; n_bins] },
                   |acc, x| add_bins(acc, x));

        format_output(&summed_result.count,
                      &summed_result.count2,
                      &domain,
                      &lower_limit,
                      &upper_limit,
                      summed_result.n_particles,
                      summed_result.dim,
                      n_ens,
                      rho
        );
    }
}
