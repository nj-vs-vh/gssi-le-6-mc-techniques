use rand::prelude::*;

use ndhistogram::{axis::UniformNoFlow, ndhistogram, Histogram};

use crate::{
    plot::{plot_histogram, XLim},
    utils::mean_std,
};

fn onedim_mc_integration(rng: &mut ThreadRng, pow: i32, sample_count: usize) -> f32 {
    let func = |x: f32| x.powi(pow);
    let n = sample_count as f32;
    let func_mean = (0..sample_count)
        .map(|_| func(rng.gen::<f32>()) / n)
        .sum::<f32>();
    func_mean // *(b - a), but in this case it's 1
}

pub fn ex5() {
    let mut rng = rand::thread_rng();

    for pow in [1, 2, 3, 4, 5] {
        for sample_count in [1_000, 10_000, 50_000] {
            let estimations_count = 10_000;
            let estimations = (0..estimations_count)
                .map(|_| onedim_mc_integration(&mut rng, pow, sample_count))
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
                Some(true_value),
            )
            .expect("Failed to plot histogram");
        }
    }
}
