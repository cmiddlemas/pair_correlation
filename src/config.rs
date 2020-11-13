use std::path::PathBuf;
use std::io::prelude::*;
use std::io::BufReader;
use std::fs::File;
use nalgebra::{Matrix2, Vector2, Matrix3, Vector3};

// Holds the result of parsing a configuration
pub struct Config {
    pub dim: usize,
    pub n_particles: usize,
    pub unit_cell: Vec<f64>,
    pub coords: Vec<f64>,
    pub diameter: Option<f64>,
}

impl Config {
    // Parses Ge's file format into a Config
    pub fn parse(path: &PathBuf) -> Config {
        let mut input = BufReader::new(File::open(path).unwrap());
        let mut buf = String::new();
        input.read_line(&mut buf).unwrap();

        let dim: usize = buf.trim().parse().unwrap();
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

        Config { dim, n_particles, unit_cell, coords, diameter: None }
    }

    pub fn parse_asc(path: &PathBuf) -> Config {
        let mut input = BufReader::new(File::open(path).unwrap());
        let mut buf = String::new();
        
        input.read_line(&mut buf).unwrap();
        let dim: usize = buf.trim()
                            .split_whitespace()
                            .next().unwrap()
                            .parse().unwrap();
        buf.clear();

        input.read_line(&mut buf).unwrap();
        let unit_cell: Vec<f64> = buf.trim()
                           .split_whitespace()
                           .map(|x| x.parse().unwrap())
                           .collect();
        
        let mut coords = Vec::new();
        let mut n_particles = 0;
        match dim {
            2 => {
                let cell_mat = Matrix2::from_column_slice(&unit_cell);
                for line in input.lines() {
                    let rel_coord = Vector2::from_iterator(
                                        line.unwrap().trim()
                                            .split_whitespace()
                                            .map(|x| x.parse::<f64>().unwrap())
                                            .take(dim)
                    );
                    let global_coord = cell_mat*rel_coord;
                    coords.extend_from_slice(global_coord.as_slice());
                    n_particles += 1;
                }
            }
            3 => {
                let cell_mat = Matrix3::from_column_slice(&unit_cell);
                for line in input.lines() {
                    let rel_coord = Vector3::from_iterator(
                                        line.unwrap().trim()
                                            .split_whitespace()
                                            .map(|x| x.parse::<f64>().unwrap())
                                            .take(dim)
                    );
                    let global_coord = cell_mat*rel_coord;
                    coords.extend_from_slice(global_coord.as_slice());
                    n_particles += 1;
                }
            }
            _ => unimplemented!(),
        }

        Config { dim, n_particles, unit_cell, coords, diameter: None }
    }

    // Only handles 3d data for now
    pub fn parse_donev(path: &PathBuf) -> Config {
        let dim: usize = 3;
        let mut input = BufReader::new(File::open(path).unwrap());
        let mut buf = String::new();

        input.read_line(&mut buf).unwrap();
        input.read_line(&mut buf).unwrap();
        buf.clear();
        input.read_line(&mut buf).unwrap();
        let n_particles: usize = buf.trim().parse().unwrap();
        buf.clear();
        input.read_line(&mut buf).unwrap();
        let diam: f64 = buf.trim().parse().unwrap();
        buf.clear();
        input.read_line(&mut buf).unwrap();
        let unit_cell: Vec<f64> = buf.trim()
                                     .split_whitespace()
                                     .map(|x| x.parse().unwrap())
                                     .collect();
        input.read_line(&mut buf).unwrap();

        let mut coords = Vec::new();
        for line in input.lines() {
            let one: Vec<f64> = line.unwrap()
                                    .trim()
                                    .split_whitespace()
                                    .map(|x| x.parse::<f64>().unwrap())
                                    .collect();
            coords.extend_from_slice(&one);
        }

        Config { dim, n_particles, unit_cell, coords, diameter: Some(diam) }
    }
}
