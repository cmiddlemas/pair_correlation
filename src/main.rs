use std::path::PathBuf;
use std::io::prelude::*;
use std::io::BufReader;
use std::fs::File;
use std::f64::consts::PI;
use structopt::StructOpt;
use rayon::prelude::*;
use itertools::izip;

/// Reads a Ge file and outputs
/// the pair correlation function on
/// stdout in format
/// domain g2 std_err(g2) poisson_std_err(g2) bin_width
/// Warning: program assumes all numeric
/// conversions are trivially valid, don't
/// put in absurdly large numbers
#[derive(StructOpt, Debug)]
#[structopt(name = "pair_correlation")]
struct Opt {
    
    /// Gives the maximum radius to sample to
    #[structopt(short, long)]
    cutoff: f64,

    /// Gives the number of bins
    #[structopt(short, long)]
    nbins: usize,

    /// Gives the value of rho to assume
    #[structopt(long, default_value = "1.0")]
    rho: f64,

    /// Gives an offset for every bin, i.e. makes a gap around the origin,
    /// to inhibit divide by s_1(r) effect
    #[structopt(short, long, default_value = "0.0")]
    offset: f64,

    /// implements a really basic autoscaling of bin size
    /// based on Poisson assumption
    #[structopt(long)]
    autoscale: bool,

    /// Gives the list of files, interpreted as ensemble average
    #[structopt(parse(from_os_str))]
    files: Vec<PathBuf>,

}

// Holds the result of parsing a configuration
struct Config {
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

// Strategy for normalization taken from Ge's code
fn sphere_vol(dim: usize, r: f64) ->f64 {
    match dim {
        1 => 2.0*r,
        2 => PI*r*r,
        3 => 4.0*PI*r*r*r/3.0,
        _ => panic!("Dimension not implemented!"),
    }
}

// Takes bin count and converts it to an approx
// to g_2 by normalization
fn format_output(count: Vec<usize>,
                 count2: Vec<usize>,
                 domain: Vec<f64>,
                 step: f64,
                 n_particles: usize,
                 dim: usize,
                 n_ens: usize,
                 rho: f64)
{
    let half_width = step/2.0;
    for (c, c2, r) in izip!(count.iter(), count2.iter(), domain.iter()) {
        // 2*count/(n_particles*n_ens) = rho s1(r) g2(r) dr
        let f_n_ens = n_ens as f64;
        let ens_c = (*c as f64)/f_n_ens;
        // Sample variance reminder:
        // https://en.wikipedia.org/wiki/Variance
        let ens_var_c = (*c2 as f64)/(f_n_ens - 1.0) - ens_c*ens_c*f_n_ens/(f_n_ens - 1.0);
        let g2_coeff = 2.0
            /((sphere_vol(dim, r + half_width) - sphere_vol(dim, r - half_width))
                *rho*(n_particles as f64)
            );
        let g2 = ens_c*g2_coeff;
        let var_g2 = ens_var_c*g2_coeff*g2_coeff;
        let std_err_g2 = var_g2.sqrt()/((n_ens as f64).sqrt());
        let poisson_est = (*c as f64).sqrt()*g2_coeff/f_n_ens;
        println!("{} {} {} {} {}", r, g2, std_err_g2, poisson_est, step);
    }
}

// Measures distance on torus defined by config.unit_cell
fn measure_distance(config: &Config,
                    i: usize,
                    j: usize
) -> f64
{
    // Comparing vs r^2 is more efficient in most dim
    let mut r2 = std::f64::MAX;
    match config.dim {
        1 => {
            for l in -1..2 {
                let candidate_x = config.coords[i]
                    - (config.coords[j] + (l as f64)*config.unit_cell[0]);
                let candidate_r2 = candidate_x.powi(2);
                if candidate_r2 < r2 {
                    r2 = candidate_r2;
                }
            }
        }
        2 => {
            for l in -1..2 {
                for m in -1..2 {
                    let candidate_x = config.coords[2*i]
                        - (config.coords[2*j] 
                            + (l as f64) * config.unit_cell[0]
                            + (m as f64) * config.unit_cell[2]
                        );
                    let candidate_y = config.coords[2*i + 1]
                        - (config.coords[2*j + 1]
                            + (l as f64) * config.unit_cell[1]
                            + (m as f64) * config.unit_cell[3]
                        );
                    let candidate_r2 = candidate_x.powi(2) + candidate_y.powi(2);
                    if candidate_r2 < r2 {
                        r2 = candidate_r2;
                    }
                }
            }
        }
        3 => {
            for l in -1..2 {
                for m in -1..2 {
                    for n in -1..2 {
                        let candidate_x = config.coords[3*i]
                            - (config.coords[3*j] 
                                + (l as f64) * config.unit_cell[0]
                                + (m as f64) * config.unit_cell[3]
                                + (n as f64) * config.unit_cell[6]
                            );
                        let candidate_y = config.coords[3*i + 1]
                            - (config.coords[3*j + 1]
                                + (l as f64) * config.unit_cell[1]
                                + (m as f64) * config.unit_cell[4]
                                + (n as f64) * config.unit_cell[7]
                            );
                        let candidate_z = config.coords[3*i + 2]
                            - (config.coords[3*j + 2]
                                + (l as f64) * config.unit_cell[2]
                                + (m as f64) * config.unit_cell[5]
                                + (n as f64) * config.unit_cell[8]
                            );
                        let candidate_r2 = candidate_x.powi(2) + candidate_y.powi(2) + candidate_z.powi(2);
                        if candidate_r2 < r2 {
                            r2 = candidate_r2;
                        }
                    }
                }
            }
        }
        _ => {
            panic!("Dimension not implemented");
        }
    }
    
    r2.sqrt()
}

fn sample_file(path: &PathBuf,
               step: f64,
               nbins: usize,
               offset: f64
) -> BinResult
{
    let config = Config::parse(path);
    
    let mut count = vec![0; nbins];

    for i in 0..config.n_particles {
        for j in 0..i {
            let r = measure_distance(&config, i, j);
            let bin = ((r-offset)/step).floor() as isize;
            if bin >= 0 && bin < (nbins as isize) {
                count[bin as usize] += 1;
            }
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

    let step: f64 = opt.cutoff/(opt.nbins as f64);

    let domain: Vec<f64> = (0..opt.nbins)
        .map(|x| step/2.0 + (x as f64)*step + opt.offset)
        .collect();

    let n_ens = opt.files.len();

    let summed_result: BinResult = opt.files.par_iter()
             .enumerate()
             .map(|(i, x)| { eprintln!("Working on {}", i); sample_file(x, step, opt.nbins, opt.offset) })
             .collect::<Vec<BinResult>>().iter()
             .fold(BinResult {dim: 0, n_particles: 0, count: vec![0; opt.nbins], count2: vec![0; opt.nbins] },
                   |acc, x| add_bins(acc, x));

    format_output(summed_result.count,
                  summed_result.count2,
                  domain,
                  step,
                  summed_result.n_particles,
                  summed_result.dim,
                  n_ens,
                  opt.rho
    );
}
