use crate::{plot::AxLim, utils::mean_std};
use ndhistogram::{axis::UniformNoFlow, ndhistogram, Histogram};
use rand::prelude::*;

use crate::plot::{plot_histogram, plot_lines, Line};

const PI32: f32 = std::f32::consts::PI;
const PI64: f64 = std::f64::consts::PI;

fn throw_into_circle(rng: &mut ThreadRng) -> bool {
    let x = 2.0 * rng.gen::<f32>() - 1.0;
    let y = 2.0 * rng.gen::<f32>() - 1.0;
    return (x.powi(2) + y.powi(2)) < 1.0;
}

fn estimate_pi(thrown_inside: usize, thrown_total: usize) -> f32 {
    4.0 * (thrown_inside as f32 / thrown_total as f32)
}

fn produce_pi_series(
    rng: &mut ThreadRng,
    total_steps: usize,
    record_each: Option<usize>,
    logging: bool,
) -> Vec<(f32, f32)> {
    let mut pi_estimates: Vec<(f32, f32)> = Vec::new();
    let mut thrown_inside: usize = 0;
    for step in 1..total_steps + 1 {
        if throw_into_circle(rng) {
            thrown_inside += 1;
        }
        if let Some(record_each) = record_each {
            if step % record_each == 0 {
                let logstep = (step as f32).log10();
                pi_estimates.push((logstep, estimate_pi(thrown_inside, step)))
            }
        }
    }

    if record_each.is_none() {
        pi_estimates.push((
            (total_steps as f32).log10(),
            estimate_pi(thrown_inside, total_steps),
        ));
    }
    if logging {
        println!(
            "Thrown {} points, recorded {} pi estimation(s)",
            total_steps,
            pi_estimates.len()
        );
    }
    pi_estimates
}

fn plot_pi_series(rng: &mut ThreadRng) {
    let steps_total: usize = 5_000_000;
    let record_each: usize = 10_000;
    let pi_samples: Vec<Vec<(f32, f32)>> = (0..3)
        .map(|_| produce_pi_series(rng, steps_total, Some(record_each), true))
        .collect();

    plot_lines(
        pi_samples
            .into_iter()
            .enumerate()
            .map(|(idx, pi_series)| Line {
                data: pi_series
                    .iter()
                    .map(|(logstep, pi_est)| (*logstep, (PI32 - pi_est).abs().log10()))
                    .collect(),
                label: format!("Run #{}", idx),
            })
            .collect(),
        "Pi value estimation convergence",
        "log_10(step)",
        "log_10(|error|)",
        "out/ex3/pi.png",
    )
    .expect("Failed to plot line");
}

fn estimate_pi_error(rng: &mut ThreadRng) {
    for total_throws in [50, 100, 500, 1000, 5000] {
        let sample_size = 1_000_000 * 100 / total_throws; // to keep time reasonable
        let pi_est_sample: Vec<f32> = (0..sample_size)
            .map(|_| {
                produce_pi_series(rng, total_throws, None, false)
                    .last()
                    .unwrap()
                    .1
            })
            .collect();

        let (mu, sigma) = mean_std(&pi_est_sample);
        let k_est = sigma * (total_throws as f32).sqrt();
        println!(
            "Sample moments: {:e} +/- {:e}, k = sigma*sqrt(N) = {:.5}",
            mu, sigma, k_est,
        );

        let sqrt_throws = (total_throws as f64).sqrt();
        // 10 / sqrt(N) scaling is ad hoc
        let low = PI64 - 10.0 / sqrt_throws;
        let high = PI64 + 10.0 / sqrt_throws;
        // it makes no sense to make more bins than max precision, limited by total throws
        let bins = ((high - low) * total_throws as f64 / 4.0) as usize;
        println!(
            "Histogram params for {} throws: range {:.4}..{:.4}, {} bins",
            total_throws, low, high, bins
        );
        let mut hist = ndhistogram!(UniformNoFlow::new(bins, low, high));
        for value in pi_est_sample {
            hist.fill(&(value as f64));
        }
        plot_histogram(
            &hist,
            &format!(
                "Pi estimations for {} throws, sample of {}: {:.3} +/- {:.3}, k = {:.3}",
                total_throws, sample_size, mu, sigma, k_est,
            ),
            "pi estimation",
            &format!("out/ex3/pi-est-distribution-{}-throws.png", total_throws),
            AxLim::FromData,
            Some(vec![PI64]),
        )
        .expect("Failed to plot histogram");
    }
}

pub fn ex3() {
    let mut rng = rand::thread_rng();
    plot_pi_series(&mut rng);
    estimate_pi_error(&mut rng);
}
