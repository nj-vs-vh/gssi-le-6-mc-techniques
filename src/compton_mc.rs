use bmp::{Image, Pixel};
use rand::prelude::*;

use ndhistogram::{axis::UniformNoFlow, ndhistogram, Histogram};

use crate::{
    ex7::sample_exponential,
    ex8::sample_compton_scattering,
    plot::{plot_lines_3d, AxLim},
};

#[derive(Debug, Clone)]
struct Vector3 {
    pub x: f32, // cm
    pub y: f32, // cm
    pub z: f32, // cm
}

impl Vector3 {
    fn normalize(&mut self) {
        let length = (self.x.powi(2) + self.y.powi(2) + self.z.powi(2)).sqrt();
        self.x /= length;
        self.y /= length;
        self.z /= length;
    }

    fn mult(&self, k: f32) -> Vector3 {
        Vector3 {
            x: k * self.x,
            y: k * self.y,
            z: k * self.z,
        }
    }

    fn add(&self, other: &Vector3) -> Vector3 {
        Vector3 {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

struct Cuboid {
    pub min: Vector3,
    pub max: Vector3,
}

impl Cuboid {
    fn contains(&self, point: &Vector3) -> bool {
        (self.min.x < point.x && point.x < self.max.x)
            && (self.min.y < point.y && point.y < self.max.y)
            && (self.min.z < point.z && point.z < self.max.z)
    }
}

struct Material {
    pub geometry: Cuboid,
    pub a: u16,       // atomic mass
    pub z: u16,       // atomic number
    pub density: f32, // g/cm^2
}

struct Photon {
    energy: f32, // KeV
    r: Vector3,
    direction: Vector3,
}

impl Photon {
    fn k(&self) -> f32 {
        self.energy / 0.511
    }
}

const R_E_SQUARED: f32 = 79.4; // mb

fn compton_integral_cross_section(photon: &Photon, material: &Material) -> f32 {
    let k = photon.k();
    let geometric_factor = std::f32::consts::PI
        * if k < 0.2 {
            (8.0 / 3.0)
                * (1.0 / (1.0 + 2.0 * k).powi(2))
                * (1.0 + 2.0 * k + 6.0 * k.powi(2) / 5.0 - k.powi(3) / 2.0 + 2.0 * k.powi(4) / 7.0
                    - 6.0 * k.powi(5) / 35.0
                    + 8.0 * k.powi(6) / 105.0
                    + 4.0 * k.powi(7) / 105.0)
        } else {
            2.0 * ((((1.0 + k) / k.powi(2))
                * ((2.0 + 2.0 * k) / (1.0 + 2.0 * k) - (1.0 + 2.0 * k).ln() / k))
                + (1.0 + 2.0 * k).ln() / (2.0 * k)
                - (1.0 + 3.0 * k) / (1.0 + 2.0 * k).powi(2))
        };
    return (material.z as f32) * R_E_SQUARED * geometric_factor; // mb
}

fn mean_free_path(photon: &Photon, material: &Material) -> f32 {
    // 1660.6 = 1 / (N_A * 1 mb)
    1660.6 * (material.a as f32)
        / ((material.z as f32)
            * material.density
            * compton_integral_cross_section(photon, material))
}

enum StepResult {
    Continue,
    ParticleExited,
    ParticleRangedOut,
}

fn compton_mc_step(
    rng: &mut ThreadRng,
    photon: &mut Photon,
    material: &Material,
    min_energy: &f32,
) -> StepResult {
    // throwing out low-energy photons that would be absorbed through photoeffect
    if photon.energy < *min_energy {
        return StepResult::ParticleRangedOut;
    }
    // moving to the next interaction point
    let lambda = mean_free_path(photon, material);
    let step_size = sample_exponential(rng, 1.0 / lambda);
    photon.direction.normalize();
    let next_interaction_point = photon.r.add(&photon.direction.mult(step_size));
    if !material.geometry.contains(&next_interaction_point) {
        return StepResult::ParticleExited;
    }
    photon.r = next_interaction_point;
    // simulating scattering
    let scattered_photon = sample_compton_scattering(rng, &photon.k());
    photon.energy *= scattered_photon.energy_ratio;
    let local_theta = scattered_photon.theta;
    let local_phi = 2.0 * std::f32::consts::PI * rng.gen::<f32>();
    // scattered particle direction in local frame, where inial particle went along Ox
    // standard vector pointing at (theta, phi)  but with z <-> x
    let (local_x, local_y, local_z) = (
        local_theta.cos(),
        local_theta.sin() * local_phi.sin(),
        local_theta.sin() * local_phi.cos(),
    );
    // rotating new direction from local frame back to lab, defined by yaw and pitch angles of initial direction
    // phi = yaw angle
    // alpha = pitch angle (theta - pi/2)
    let sinalpha = -photon.direction.z;
    let cosalpha = (1.0 - sinalpha.powi(2)).sqrt();
    let cosphi = photon.direction.x / cosalpha;
    let sinphi = photon.direction.y / cosalpha;
    // coordinatwize matrix multiplication
    photon.direction.x =
        cosphi * cosalpha * local_x - sinphi * local_y + cosphi * sinalpha * local_z;
    photon.direction.y =
        sinphi * cosalpha * local_x + cosphi * local_y + sinphi * sinalpha * local_z;
    photon.direction.z = -sinalpha * local_x + cosalpha * local_z;
    StepResult::Continue
}

pub fn main() {
    // no absorption = photons do a long random walk before going out
    run_compton_mc(100, 500.0, 0.0);
    // more or less realistic
    run_compton_mc(1_000_000, 200.0, 1.0);
    run_compton_mc(1_000_000, 800.0, 1.0);
    run_compton_mc(1_000_000, 3200.0, 1.0);
    run_compton_mc(1_000_000, 6400.0, 1.0);
    run_compton_mc(1_000_000, 12800.0, 1.0);
}

fn run_compton_mc(total_photons: usize, initial_energy: f32, min_energy: f32) {
    let material = Material {
        geometry: Cuboid {
            min: Vector3 {
                x: 0.0,
                y: -2.0,
                z: -2.0,
            },
            max: Vector3 {
                x: 4.0,
                y: 2.0,
                z: 2.0,
            },
        },
        z: 82,
        a: 208,
        density: 11.384,
    };

    let mut rng = thread_rng();

    let backside_image_size_px = 100;
    let mut backside_image_hist = ndhistogram!(
        UniformNoFlow::new(
            backside_image_size_px,
            material.geometry.min.y,
            material.geometry.max.y
        ),
        UniformNoFlow::new(
            backside_image_size_px,
            material.geometry.min.z,
            material.geometry.max.z
        )
    );

    let mut traces: Vec<Vec<(f32, f32, f32)>> = Vec::new();
    let mut exited_photons: usize = 0;
    for _idx in 0..total_photons {
        let mut photon = Photon {
            energy: initial_energy,
            r: Vector3 {
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
            direction: Vector3 {
                x: 1.0,
                y: 0.0,
                z: 0.0,
            },
        };
        let mut trace: Vec<Vector3> = Vec::new();
        trace.push(photon.r.clone());
        loop {
            let res = compton_mc_step(&mut rng, &mut photon, &material, &min_energy);
            trace.push(photon.r.clone());
            if let StepResult::Continue = res {
                continue;
            }
            if let StepResult::ParticleExited = res {
                exited_photons += 1;
                // adding exited photon to the histogram
                if let Some(backside_intersection) = find_ray_x_plane_intersection(
                    &photon.r,
                    &photon.direction,
                    &material.geometry.max.x,
                ) {
                    // println!("{} {}", backside_intersection.y, backside_intersection.z);
                    backside_image_hist.fill(&(backside_intersection.y, backside_intersection.z));
                }
            }
            break;
        }

        traces.push(trace.into_iter().map(|r| (r.x, r.y, r.z)).collect());
    }
    println!("Total particles exited from material: {}", exited_photons);

    plot_lines_3d(
        traces,
        &format!(
            "Compton-scattered photons (E={:.0} KeV, N={}, Emin={:.1}, absorbed {:.2}%)",
            initial_energy,
            total_photons,
            min_energy,
            100.0 * (1.0 - (exited_photons as f32) / (total_photons as f32))
        ),
        (
            AxLim::Range(
                material.geometry.min.x as f64,
                material.geometry.max.x as f64,
            ),
            AxLim::Range(
                material.geometry.min.y as f64,
                material.geometry.max.y as f64,
            ),
            AxLim::Range(
                material.geometry.min.z as f64,
                material.geometry.max.z as f64,
            ),
        ),
        &format!(
            "out/compton/e={:.2}-emin={:.2}-3d-traces.png",
            initial_energy, min_energy
        ),
    )
    .expect("Failed to plot photon traces");

    let max_bin_count = backside_image_hist
        .iter()
        .map(|item| item.value)
        .fold(f64::MIN, |acc, v| v.max(acc)) as f32;
    let mut img = Image::new(backside_image_size_px as u32, backside_image_size_px as u32);
    let y_px_size =
        (material.geometry.max.y - material.geometry.min.y) / (backside_image_size_px as f32);
    let z_px_size =
        (material.geometry.max.z - material.geometry.min.z) / (backside_image_size_px as f32);
    for (y_px, z_px) in img.coordinates() {
        let pixel_count = backside_image_hist
            .value(&(
                (material.geometry.min.y + (y_px as f32) * y_px_size),
                (material.geometry.min.z + (z_px as f32) * z_px_size),
            ))
            .unwrap()
            .clone() as f32;
        let pixel_intensity = (1.0 + pixel_count).ln() / (1.0 + max_bin_count).ln();
        let pixel_value = (255.0 * pixel_intensity) as u32;
        img.set_pixel(y_px, z_px, px!(pixel_value, pixel_value, pixel_value));
    }
    let _ = img.save(&format!(
        "out/compton/e={:.2}-emin={:.2}-backside-image.png",
        initial_energy, min_energy
    ));
}

fn find_ray_x_plane_intersection(
    ray_start: &Vector3,
    ray_dir: &Vector3,
    plane_x: &f32,
) -> Option<Vector3> {
    let normal_dir_dot = ray_dir.x; // normal = (1, 0, 0)
    if normal_dir_dot.abs() < 1e-15 {
        // = ray is parallel to the plane (with some tolerance)
        return None;
    }

    let t_intersection = (plane_x - ray_start.x) / normal_dir_dot;
    if !t_intersection.is_finite() {
        // something's not right!
        return None;
    }
    if t_intersection < 0.0 {
        // = ray goes away from the plane
        return None;
    }

    Some(Vector3 {
        x: ray_start.x + t_intersection * ray_dir.x,
        y: ray_start.y + t_intersection * ray_dir.y,
        z: ray_start.z + t_intersection * ray_dir.z,
    })
}
