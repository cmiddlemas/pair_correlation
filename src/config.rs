use std::path::PathBuf;
use std::io::prelude::*;
use std::io::BufReader;
use std::fs::File;

// Holds the result of parsing a configuration
pub struct Config {
    pub dim: usize,
    pub n_particles: usize,
    pub unit_cell: Vec<f64>,
    pub coords: Vec<f64>
}

impl Config {
    // Parses Ge's file format into a Config
    pub fn parse(path: &PathBuf) -> Config {
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

        Config { dim, n_particles, unit_cell, coords }
    }

    pub fn parse_asc(path: &PathBuf) -> Config {
        let mut input = BufReader::new(File::open(path).unwrap());
        let mut line = String::new();
        input.read_line(&mut line).unwrap();

        let dim: usize = line.trim()
                             .split_whitespace()
                             .next().unwrap()
                             .parse().unwrap();
        input.read_line(&mut line).unwrap();
        let unit_cell = line.trim()
                            .split_whitespace()
                            .map(|x| x.parse().unwrap())
                            .collect();
        
        let mut coords = Vec::new();
        let mut n_particles = 0;
        for line in input.lines() {
            coords.extend(
                line.unwrap().trim().split_whitespace().map(|x| x.parse::<f64>().unwrap())
            );
            n_particles += 1;
        }
        
        Config { dim, n_particles, unit_cell, coords }
    }
}
