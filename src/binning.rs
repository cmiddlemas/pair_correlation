use std::path::PathBuf;
use crate::config::Config;
use crate::measure_distance::measure_distance;
use itertools::izip;

// Holds the data from binning a configuration
pub struct BinResult {
    pub dim: usize,
    pub n_particles: usize,
    pub count: Vec<usize>,
    pub count2: Vec<usize>, // accumulator for square of count results
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

pub fn sample_file(path: &PathBuf,
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
pub fn add_bins(mut acc: BinResult, x: &BinResult) -> BinResult {
    for (x, x2, y) in izip!(acc.count.iter_mut(), acc.count2.iter_mut(), x.count.iter()) {
        *x += y;
        *x2 += y*y;
    }
    acc.dim = x.dim;
    acc.n_particles = x.n_particles;
    acc
}

