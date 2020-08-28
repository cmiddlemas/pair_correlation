use std::path::PathBuf;
use crate::Opt;
use crate::config::Config;
use crate::measure_distance::measure_distance;
use crate::simple_math::{sphere_vol, sphere_radius, cell_volume};
use itertools::izip;

// Returns (domain, lower_limit, upper_limit),
// where domain is given as the midpoint of the bin
pub fn make_bins(opt: &Opt, first_config: &Config)
    -> (Vec<f64>, Vec<f64>, Vec<f64>)
{
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
    
    } else if opt.logarithm {
        
        let log_interval = (opt.cutoff - opt.offset).log10();
        let log_step = log_interval/(opt.nbins as f64);

        domain = (0..opt.nbins)
            .map(|x| 10.0f64.powf(log_step/2.0 + (x as f64)*log_step) + opt.offset)
            .collect();
        lower_limit = (0..opt.nbins)
            .map(|x| 10.0f64.powf((x as f64)*log_step) + opt.offset)
            .collect();
        upper_limit = (0..opt.nbins)
            .map(|x| 10.0f64.powf(((x+1) as f64)*log_step) + opt.offset)
            .collect();

    } else { // Equal width bins
        
        let step: f64 = (opt.cutoff - opt.offset)/(opt.nbins as f64);

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
    (domain, lower_limit, upper_limit)
}

// Holds the data from binning a configuration
pub struct BinResult {
    pub dim: usize,
    pub n_obs: usize,
    pub rho: f64,
    pub count: Vec<usize>,
    pub count2: Vec<usize>, // accumulator for square of count results
}

// use binary search
// will panic if r is not in any bin
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
            if r >= l {
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
               upper_limit: &[f64],
               opt: &Opt
) -> BinResult
{
    let config = if opt.asc {
        Config::parse_asc(path)
    } else {
        Config::parse(path)
    };
    
    let mut count = vec![0; lower_limit.len()];

    for i in 0..config.n_particles {
        for j in 0..i {
            let r = measure_distance(&config, i, j);
            if r >= lower_limit[0] && r <= upper_limit[upper_limit.len() - 1] {
                bin_distance(r, &mut count, lower_limit, upper_limit);
            }
        }
    }

    let c_vol = cell_volume(config.dim, &config.unit_cell);

    BinResult {
        dim: config.dim,
        n_obs: config.n_particles,
        rho: (config.n_particles as f64)/c_vol,
        count,
        count2: Vec::new()
    }
}

// Takes the right term as the truth for dim and n_particles
// Sums the histogram
pub fn add_bins(mut acc: BinResult, b: &BinResult) -> BinResult {
    for (x, x2, y) in izip!(acc.count.iter_mut(), acc.count2.iter_mut(), b.count.iter()) {
        *x += y;
        *x2 += y*y;
    }
    acc.dim = b.dim;
    acc.n_obs += b.n_obs;
    acc.rho += b.rho;
    acc
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[should_panic]
    fn binary_search1() {
        let lower = [1.0, 2.0, 3.0, 4.0];
        let upper = [2.0, 3.0, 4.0, 5.0];
        let r = 0.9999999;
        let mut count = [0usize, 0usize, 0usize];
        bin_distance(r, &mut count, &lower, &upper);
        assert!(count == [0, 0, 0]);
    }
    
    #[test]
    #[should_panic]
    fn binary_search2() {
        let lower = [1.0, 2.0, 3.0];
        let upper = [2.0, 3.0, 4.0];
        let r = 5.0;
        let mut count = [0usize, 0usize, 0usize];
        bin_distance(r, &mut count, &lower, &upper);
        assert!(count == [0, 0, 0]);
    }
    
    #[test]
    fn binary_search3() {
        let lower = [1.0, 2.0, 3.0];
        let upper = [2.0, 3.0, 4.0];
        let r = 3.5;
        let mut count = [0usize, 0usize, 0usize];
        bin_distance(r, &mut count, &lower, &upper);
        assert!(count == [0, 0, 1]);
    }

    #[test]
    #[should_panic]
    fn binary_search4() {
        let lower: Vec<f64> = (1..101).map(|x| x as f64).collect();
        let upper: Vec<f64> = (2..102).map(|x| x as f64).collect();
        let r = 0.9;
        let mut count = vec![0usize; 100];
        bin_distance(r, &mut count, &lower, &upper);
    }
    
    #[test]
    #[should_panic]
    fn binary_search5() {
        let lower: Vec<f64> = (1..103).map(|x| x as f64).collect();
        let upper: Vec<f64> = (2..104).map(|x| x as f64).collect();
        let r = 104.1;
        let mut count = vec![0usize; 102];
        bin_distance(r, &mut count, &lower, &upper);
    }

    #[test]
    fn binary_search6() {
        let lower = [1.0, 2.0, 3.0];
        let upper = [2.0, 3.0, 4.0];
        let r = 1.5;
        let mut count = [0usize, 0usize, 0usize];
        bin_distance(r, &mut count, &lower, &upper);
        assert!(count == [1, 0, 0]);
    }

}
