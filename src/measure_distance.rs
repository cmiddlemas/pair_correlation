use crate::Config;

// Measures distance on torus defined by config.unit_cell
pub fn measure_distance(config: &Config,
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
