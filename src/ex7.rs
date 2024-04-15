use ndhistogram::{axis::UniformNoFlow, ndhistogram, Histogram};
use rand::prelude::*;

use crate::plot::plot_histogram;

fn sample_exponential(rng: &mut ThreadRng, mu: f32) -> f32 {
    -(rng.gen::<f32>()).ln() / mu
}

fn sample_min_exponential1(rng: &mut ThreadRng, mu1: f32, mu2: f32) -> f32 {
    sample_exponential(rng, mu1).min(sample_exponential(rng, mu2))
}

fn sample_min_exponential2(rng: &mut ThreadRng, mu1: f32, mu2: f32) -> f32 {
    sample_exponential(rng, mu1 + mu2)
}

pub fn ex7() {
    let mut rng = thread_rng();
    let mu1 = 1.0;
    let mu2 = 2.0;

    let sample_size = 1_000_000;

    for (idx, func) in [sample_min_exponential1, sample_min_exponential2]
        .into_iter()
        .enumerate()
    {
        let time = std::time::Instant::now();
        let mut hist = ndhistogram!(UniformNoFlow::new(100, 0.0, 3.0));
        for _ in 0..sample_size {
            hist.fill(&(func(&mut rng, mu1, mu2) as f64));
        }
        println!(
            "sampling with method #{} took {:.3} sec",
            idx,
            time.elapsed().as_secs_f32()
        );

        plot_histogram(
            &hist,
            &format!("Min of exponentials, method {}", idx + 1),
            "s = min(s1, s2)",
            &format!("out/ex7/sum-exp-sample-{}.png", idx + 1),
            crate::plot::XLim::FromData,
            None,
        )
        .expect("Failed to plot histogram");
    }

    let mu_process = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    let mut counts = Vec::<usize>::new();
    (0..mu_process.len()).for_each(|_| counts.push(0));
    let sample_size = 1_000_000;
    for _ in 0..sample_size {
        let (argmin, _) = mu_process
            .iter()
            .map(|mu_i| sample_exponential(&mut rng, *mu_i))
            .enumerate()
            .fold(
                (None, std::f32::INFINITY),
                |(argmin, min), (idx, current)| {
                    if current < min {
                        (Some(idx), current)
                    } else {
                        (argmin, min)
                    }
                },
            );
        if let Some(argmin) = argmin {
            counts[argmin] += 1;
        }
    }

    let mu = mu_process.iter().sum::<f32>();
    for (idx, (count, mu_i)) in counts.into_iter().zip(mu_process.iter()).enumerate() {
        println!(
            "process #{}: {:.3}% of samples, 100*mu_i/mu = {:3}",
            idx,
            100.0 * (count as f32) / (sample_size as f32),
            100.0 * mu_i / mu,
        )
    }
}
