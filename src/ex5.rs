use itertools::Itertools;
use ndhistogram::{axis::UniformNoFlow, ndhistogram, Histogram};
use rand::prelude::*;

use crate::{
    plot::{plot_histogram, XLim},
    utils::mean_std,
};

fn onedim_mc_integral(rng: &mut ThreadRng, pow: i32, sample_count: usize) -> f32 {
    let func = |x: f32| x.powi(pow);
    let n = sample_count as f32;
    let func_mean = (0..sample_count)
        .map(|_| func(rng.gen::<f32>()) / n)
        .sum::<f32>();
    func_mean // *(b - a), but in this case it's 1
}

fn plot_onedim_integral_mc_distributions(rng: &mut ThreadRng) {
    for pow in [1, 2, 3, 4, 5] {
        for sample_count in [1_000, 10_000, 100_000] {
            let estimations_count = 10_000;
            let estimations = (0..estimations_count)
                .map(|_| onedim_mc_integral(rng, pow, sample_count))
                .collect::<Vec<_>>();
            let (mu, sigma) = mean_std(&estimations);
            let k_est = sigma * (sample_count as f32).sqrt();

            let true_value = 1.0 / (1.0 + pow as f64);

            let lo = (mu - 5.0 * sigma) as f64;
            let hi = (mu + 5.0 * sigma) as f64;
            let mut hist = ndhistogram!(UniformNoFlow::new(100, lo, hi));
            for value in estimations {
                hist.fill(&(value as f64));
            }
            plot_histogram(
                &hist,
                &format!(
                    "One-dimensional MC integration with {} samples: {:.4}+/-{:.4}, k = {:.3}",
                    sample_count, mu, sigma, k_est
                ),
                &format!("\\int_0^1 x^{} dx", pow),
                &format!(
                    "out/ex5/1/x^{}-mc-integral-{}-samples.png",
                    pow, sample_count
                ),
                XLim::Range(lo, hi),
                Some(vec![true_value]),
            )
            .expect("Failed to plot histogram");
        }
    }
}

enum Integrand {
    SquaresSum,
    ExpProduct,
}

impl Integrand {
    pub fn name(&self) -> &str {
        match self {
            Integrand::SquaresSum => "squares-sum",
            Integrand::ExpProduct => "exp-product",
        }
    }

    pub fn true_integral(&self, ndim: usize) -> f64 {
        match self {
            Integrand::SquaresSum => (ndim as f64) / 3.0,
            Integrand::ExpProduct => (1.0 - 1.0 / std::f64::consts::E).powi(ndim as i32),
        }
    }

    pub fn random_value(&self, rng: &mut ThreadRng, ndim: usize) -> f32 {
        let random_point: Vec<f32> = (0..ndim).map(|_| rng.gen::<f32>()).collect();
        self.value_at(random_point.iter())
    }

    pub fn value_at<'a>(&self, point: impl Iterator<Item = &'a f32>) -> f32 {
        match self {
            Integrand::SquaresSum => point.map(|coord| (*coord).powi(2)).sum::<f32>(),
            Integrand::ExpProduct => point.map(|coord| (-*coord).exp()).product::<f32>(),
        }
    }
}

fn ndim_mc_integral(
    rng: &mut ThreadRng,
    ndim: usize,
    sample_count: usize,
    integrand: &Integrand,
) -> f32 {
    let n = sample_count as f32;
    let func_mean = (0..sample_count)
        .map(|_| integrand.random_value(rng, ndim) / n)
        .sum::<f32>();
    func_mean // times hypercube volume, but in this case it's 1
}

fn ndim_midpoint_integral(ndim: usize, cells_per_dim: usize, integrand: &Integrand) -> f32 {
    let step = 1.0 / (cells_per_dim as f32);
    let cell_volume = step.powi(ndim as i32);
    let onedim_cell_centers = (0..cells_per_dim)
        .map(|idx| step * (idx as f32 + 0.5))
        .collect::<Vec<_>>();
    let cell_centers = (0..ndim)
        .map(|_| onedim_cell_centers.iter())
        .multi_cartesian_product();
    cell_centers
        .map(|cell| cell_volume * integrand.value_at(cell.into_iter()))
        .sum()
}

fn integration_params(ndim: &usize) -> (usize, usize) {
    let cells_per_dim = match *ndim {
        1 => 65536,
        2 => 256,
        3 => 40,
        4 => 16,
        5 => 9,
        6 => 6,
        7 => 5,
        8 => 4,
        _ => panic!("Unexpected dimension number"),
    };
    (cells_per_dim, cells_per_dim.pow(*ndim as u32))
}

fn plot_ndim_integral_mc_vs_midpoint(rng: &mut ThreadRng, integrand: Integrand) {
    for ndim in 1..=8 {
        let (cells_per_dim, sample_count) = integration_params(&ndim);
        let mut timer = std::time::Instant::now();
        let int_midpoint = ndim_midpoint_integral(ndim, cells_per_dim, &integrand) as f64;
        println!(
            "Midpoint integration took {:.4} sec",
            timer.elapsed().as_secs_f32()
        );

        let evals: usize = 500;
        timer = std::time::Instant::now();
        let int_mc_evals: Vec<f32> = (0..evals)
            .map(|_| ndim_mc_integral(rng, ndim, sample_count, &integrand))
            .collect();

        let time = timer.elapsed().as_secs_f32();
        println!(
            "{} MC integrations took {:.2} sec, on average {:.4} sec",
            evals,
            time,
            time / (evals as f32)
        );
        let (mu, sigma) = mean_std(&int_mc_evals);
        let lo = (mu - 5.0 * sigma) as f64;
        let hi = (mu + 5.0 * sigma) as f64;
        let mut hist = ndhistogram!(UniformNoFlow::new(30, lo, hi));
        for value in int_mc_evals {
            hist.fill(&(value as f64));
        }
        let true_value = integrand.true_integral(ndim);

        plot_histogram(
            &hist,
            &format!(
                "MC integration vs midpoint integration for {} dimensions",
                ndim
            ),
            "integral value",
            &format!(
                "out/ex5/2/{}-mc-vs-midpoint-{}-dim.png",
                integrand.name(),
                ndim
            ),
            XLim::enlarged_range(
                lo.min(true_value.min(int_midpoint)),
                hi.max(true_value.max(int_midpoint)),
                0.05,
            ),
            Some(vec![true_value, int_midpoint]),
        )
        .expect("Failed to plot histogram");
    }
}

pub fn ex5() {
    let mut rng = rand::thread_rng();

    // plot_onedim_integral_mc_distributions(&mut rng);
    plot_ndim_integral_mc_vs_midpoint(&mut rng, Integrand::SquaresSum);
    plot_ndim_integral_mc_vs_midpoint(&mut rng, Integrand::ExpProduct);
}
