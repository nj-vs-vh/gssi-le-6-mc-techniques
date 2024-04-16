use rand::prelude::*;

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

    fn add_inplace(&mut self, other: &Vector3) {
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
    }
}

struct Cuboid {
    pub min: Vector3,
    pub max: Vector3,
}

impl Cuboid {
    fn is_inside(&self, point: &Vector3) -> bool {
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

fn step(rng: &mut ThreadRng, photon: &mut Photon, material: &Material) -> bool {
    // moving to the next interaction point
    let lambda = mean_free_path(photon, material);
    let step_size = sample_exponential(rng, 1.0 / lambda);
    photon.direction.normalize();
    photon.r.add_inplace(&photon.direction.mult(step_size));
    if !material.geometry.is_inside(&photon.r) {
        return false;
    }
    // simulating scattering
    let scattered_photon = sample_compton_scattering(rng, &photon.k());
    photon.energy *= scattered_photon.energy_ratio;
    let local_theta = scattered_photon.theta;
    let local_phi = std::f32::consts::PI * rng.gen::<f32>();
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
    return true;
}

pub fn run_compton_mc() {
    let material = Material {
        geometry: Cuboid {
            min: Vector3 {
                x: 0.0,
                y: -5.0,
                z: -5.0,
            },
            max: Vector3 {
                x: 10.0,
                y: 5.0,
                z: 5.0,
            },
        },
        z: 82,
        a: 208,
        density: 11.384,
    };

    let mut rng = thread_rng();

    let mut traces: Vec<Vec<(f32, f32, f32)>> = Vec::new();
    for _ in 0..10 {
        let mut photon = Photon {
            energy: 3.0 * 511.0,
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
        let mut tracing = true;
        trace.push(photon.r.clone());
        while tracing {
            tracing = step(&mut rng, &mut photon, &material);
            trace.push(photon.r.clone());
        }
        println!("Steps made: {}", trace.len());
        traces.push(trace.into_iter().map(|r| (r.x, r.y, r.z)).collect());
    }

    plot_lines_3d(
        traces,
        "Compton-scattered photon trace",
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
        "out/compton/traces.png",
    )
    .expect("Failed to plot photon traces");
}
