use std::path::PathBuf;
use std::io::prelude::*;
use std::io::BufReader;
use std::fs::File;
use structopt::StructOpt;
use rayon::prelude::*;
use itertools::{izip, Itertools};

mod measure_distance;
mod simple_math;

use measure_distance::measure_distance;
use simple_math::{cell_volume, sphere_vol, sphere_radius};

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
struct Opt {
   
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

// Holds the result of parsing a configuration
pub struct Config {
    dim: usize,
    n_particles: usize,
    unit_cell: Vec<f64>,
    coords: Vec<f64>
}

impl Config {
    // Parses Ge's file format into a Config
    fn parse(path: &PathBuf) -> Config {
        let mut input = BufReader::new(File::open(path).unwrap());
        let mut line = String::new();
        input.read_line(&mut line).unwrap();

        let dim: usize = line.trim().parse().unwrap();
        let mut unit_cell = Vec::new();
        let mut coords = Vec::new();
        let mut n_particles = 0;
        
        for (i, line) in input.lines().enumerate() {
            if i < dim {
                for (j, token) in line.unwrap().split_whitespace().enumerate() {
                    if j < dim {
                        unit_cell.push(token.parse().unwrap());
                    }
                }
            } else {
                for (j, token) in line.unwrap().split_whitespace().enumerate() {
                    if j < dim {
                        coords.push(token.parse().unwrap());
                    }
                }
                n_particles += 1;
            }
        }

        Config { dim, n_particles, unit_cell , coords }
    }
}

// Holds the data from binning a configuration
struct BinResult {
    dim: usize,
    n_particles: usize,
    count: Vec<usize>,
    count2: Vec<usize>, // accumulator for square of count results
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

// use binary search
fn bin_distance(r: f64, count: &mut [usize], lower_limit: &[f64], upper_limit: &[f64]) {
    let length = count.len();
    let mut lower_idx = 0;
    let mut upper_idx = length - 1;
    let mut idx = (upper_idx - lower_idx)/2;
    loop {
        let l = lower_limit[idx];
        let u = upper_limit[idx];
        if r >= l && r < u {
            count[idx] += 1;
            break;
        } else { //continue search
            if lower_idx == upper_idx {
                break; // r not in binning range
            } else if r >= l {
                lower_idx = idx + 1;
                idx = lower_idx + (upper_idx - lower_idx)/2;
            } else {
                upper_idx = idx - 1;
                idx = lower_idx + (upper_idx - lower_idx)/2;
            }
        }
    }
}

fn sample_file(path: &PathBuf,
               lower_limit: &[f64],
               upper_limit: &[f64]
) -> BinResult
{
    let config = Config::parse(path);
    
    let mut count = vec![0; lower_limit.len()];

    for i in 0..config.n_particles {
        for j in 0..i {
            let r = measure_distance(&config, i, j);
            bin_distance(r, &mut count, lower_limit, upper_limit);
        }
    }

    BinResult {dim: config.dim, n_particles: config.n_particles, count, count2: Vec::new() }
}

// Takes the right term as the truth for dim and n_particles
// Sums the histogram
fn add_bins(mut acc: BinResult, x: &BinResult) -> BinResult {
    for (x, x2, y) in izip!(acc.count.iter_mut(), acc.count2.iter_mut(), x.count.iter()) {
        *x += y;
        *x2 += y*y;
    }
    acc.dim = x.dim;
    acc.n_particles = x.n_particles;
    acc
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

    let rho: f64 = (first_config.n_particles as f64)
                    / (cell_volume(first_config.dim, &first_config.unit_cell));

    let domain: Vec<f64>;
    let lower_limit: Vec<f64>;
    let upper_limit: Vec<f64>;
    
    if let Some(first_bin_width) = opt.autoscale { //Bins that scale with surface area of sphere
        
        let dim = first_config.dim;
        
        let bin_vol: f64 = sphere_vol(dim, opt.offset+first_bin_width)
            - sphere_vol(dim, opt.offset);

        let mut r_vec = vec![opt.offset, opt.offset+first_bin_width];
        while *r_vec.last().unwrap() < opt.cutoff {
            let outer_vol = bin_vol + sphere_vol(dim, *r_vec.last().unwrap());
            r_vec.push(sphere_radius(dim, outer_vol));
        }
        r_vec.pop(); // went one too high in previous loop

        // https://stackoverflow.com/questions/54273751/rust-and-vec-iterator-how-to-filter
        lower_limit = r_vec[0..r_vec.len()-1].iter().cloned().collect();
        upper_limit = r_vec[1..r_vec.len()].iter().cloned().collect();
        domain = (lower_limit.iter()).zip(upper_limit.iter())
                               .map(|(&l, &u)| (l+u)/2.0)
                               .collect();

        if opt.verbosity > 0 {
            eprintln!("Using {} bins", domain.len());
        }
    
    } else { // Equal width bins
        
        let step: f64 = opt.cutoff/(opt.nbins as f64);

        domain = (0..opt.nbins)
        .map(|x| step/2.0 + (x as f64)*step + opt.offset)
        .collect();
        
        lower_limit = (0..opt.nbins)
        .map(|x| (x as f64)*step + opt.offset)
        .collect();

        upper_limit = (0..opt.nbins)
        .map(|x| ((x+1) as f64)*step + opt.offset)
        .collect();
        
    }

    let n_ens = opt.files.len();

    let n_bins = domain.len();

    let individual_results: Vec<BinResult> = opt.files.par_iter()
            .enumerate()
            .map(|(i, x)| { 
                if opt.verbosity > 1 { eprintln!("Working on {}", i); };
                sample_file(x, &lower_limit, &upper_limit) 
            })
            .collect();


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

    } else {
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
