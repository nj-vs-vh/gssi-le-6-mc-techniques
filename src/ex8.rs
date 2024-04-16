use itertools::Itertools;
use ndhistogram::{axis::UniformNoFlow, ndhistogram, Histogram};
use rand::prelude::*;

// use std::fs::File;
// use std::io::prelude::*;

use crate::plot::{plot_histogram, plot_lines, Line};

fn compton_energy_ratio(costheta: &f32, k: &f32) -> f32 {
    1.0 / (1.0 + k * (1.0 - costheta))
}

// d sigma / d Omega, normalized to unity at theta = 0
fn klein_nishina_cross_section(costheta: &f32, k: &f32) -> f32 {
    let e_ratio = compton_energy_ratio(costheta, k);
    0.5 * e_ratio.powi(2) * (e_ratio + (1.0 / e_ratio) - (1.0 - costheta.powi(2)))
}

struct ScatteredPhoton {
    pub theta: f32,
    pub energy_ratio: f32,
}

fn sample_compton_scattering(rng: &mut ThreadRng, k: &f32) -> ScatteredPhoton {
    let costheta = loop {
        // sampling in cos(theta) space, so no need for sin(theta) jacobian
        let costheta_try = 1.0 - 2.0 * rng.gen::<f32>();
        let kn_value: f32 = rng.gen();
        if kn_value < klein_nishina_cross_section(&costheta_try, k) {
            break costheta_try;
        }
    };
    ScatteredPhoton {
        theta: costheta.acos(),
        energy_ratio: compton_energy_ratio(&costheta, k),
    }
}

pub fn ex8() {
    let mut rng = thread_rng();
    for k in [0.01, 0.1, 0.5, 1.0, 10.0, 100.0, 1000.0] {
        plot_compton_scattering_sample(&mut rng, &k);
    }
}

fn plot_compton_scattering_sample(rng: &mut ThreadRng, k: &f32) {
    let primary_energy_kev = k * 511.0;
    let costheta = (0..1000)
        .map(|i| 1.0 - 2.0 * (i as f32) / 1000.0)
        .collect_vec();
    let kn = costheta
        .iter()
        .map(|th| klein_nishina_cross_section(th, &k))
        .collect_vec();
    plot_lines(
        vec![Line {
            data: costheta
                .iter()
                .zip(kn.iter())
                .map(|(rx, ry)| (*rx, *ry))
                .collect_vec(),
            label: "Klein-Nishina cross section".to_owned(),
        }],
        &format!("Klein-Nishina cross section for E = {:.1} m_e c^2", k),
        "cos(theta)",
        "normalized cross-section, r_e^2",
        &format!("out/ex8/k={:.2}-klein-nishina.png", k),
    )
    .expect("Failed to plot KN cross section plot");

    let sample_size = 1_000_000;
    let scattered_photons = (0..sample_size)
        .map(|_| sample_compton_scattering(rng, &k))
        .collect_vec();

    // plotting theta distribution
    let mut theta_hist = ndhistogram!(UniformNoFlow::new(100, 0.0, std::f64::consts::PI));
    for photon in scattered_photons.iter() {
        theta_hist.fill(&(photon.theta as f64));
    }
    plot_histogram(
        &theta_hist,
        &format!(
            "Compton-scattered photon angles distribution for E = {:.2} KeV",
            primary_energy_kev,
        ),
        "theta",
        &format!("out/ex8/k={:.2}-theta-dist.png", k),
        crate::plot::XLim::FromData,
        None,
    )
    .expect("Failed to plot theta histogram");

    // plotting E distribution
    let mut energy_hist = ndhistogram!(UniformNoFlow::new(
        100,
        compton_energy_ratio(&-1.0, &k) as f64,
        compton_energy_ratio(&1.0, &k) as f64,
    ));
    for photon in scattered_photons.iter() {
        energy_hist.fill(&(photon.energy_ratio as f64));
    }
    plot_histogram(
        &energy_hist,
        &format!(
            "Compton-scattered photon energy distribution for E = {:.2} KeV",
            primary_energy_kev
        ),
        "E_fin / E_in",
        &format!("out/ex8/k={:.2}-E-dist.png", k),
        crate::plot::XLim::FromData,
        None,
    )
    .expect("Failed to plot E histogram");
}
