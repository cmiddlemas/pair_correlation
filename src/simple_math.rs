use nalgebra::{Matrix2, Matrix3};
use std::f64::consts::PI;

pub fn cell_volume(dim: usize, unit_cell: &[f64]) -> f64 {
    match dim {
        1 => unit_cell[0],
        2 => {
            let cell_mat = Matrix2::from_column_slice(unit_cell);
            cell_mat.determinant()
        }
        3 => {
            let cell_mat = Matrix3::from_column_slice(unit_cell);
            cell_mat.determinant()
        }
        _ => panic!("Dimension not implemented!"),
    }
}

pub fn sphere_vol(dim: usize, r: f64) -> f64 {
    match dim {
        1 => 2.0*r,
        2 => PI*r*r,
        3 => 4.0*PI*r*r*r/3.0,
        _ => panic!("Dimension not implemented!"),
    }
}

pub fn sphere_radius(dim: usize, vol: f64) -> f64{
    match dim {
        1 => 0.5*vol,
        2 => (vol/PI).sqrt(),
        3 => (3.0*vol/(4.0*PI)).powf(1.0/3.0),
        _ => panic!("Dimension not implemented!"),
    }
}
