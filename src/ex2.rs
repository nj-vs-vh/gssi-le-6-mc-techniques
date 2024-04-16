use ndhistogram::{axis::UniformNoFlow, ndhistogram, Histogram};
use rand::prelude::*;

use crate::plot::{plot_histogram, AxLim};

fn plot_cauchy_dist() {
    let mut rng = rand::thread_rng();
    let sample_size = 10000;
    let sample: Vec<f64> = (0..sample_size)
        // inverse CDF C^-1(t) = tan(pi*(t - 1/2))
        .map(|_| (std::f64::consts::PI * (rng.gen::<f64>() - 0.5)).tan())
        .collect();
    let mut hist = ndhistogram!(UniformNoFlow::new(100, -10.0, 10.0));
    sample.iter().map(|s| hist.fill(s)).count();
    plot_histogram(
        &hist,
        "Cauchy distribution sample",
        "x",
        "out/ex2/cauchy.png",
        AxLim::Range(-10.0, 10.0),
        Some(vec![0.0]),
    )
    .expect("Failed to plot histogram");
}

const TWO_PI: f32 = 2.0 * std::f32::consts::PI;

fn sample_circle_analytic(rng: &mut ThreadRng) -> (f32, f32) {
    let t1: f32 = rng.gen();
    let t2: f32 = rng.gen();
    let sqrt_t1 = t1.sqrt();
    return (sqrt_t1 * (TWO_PI * t2).cos(), sqrt_t1 * (TWO_PI * t2).sin());
}

fn sample_circle_rejection(rng: &mut ThreadRng) -> (f32, f32) {
    loop {
        let x = 2.0 * rng.gen::<f32>() - 1.0;
        let y = 2.0 * rng.gen::<f32>() - 1.0;
        if x.powi(2) + y.powi(2) <= 1.0 {
            return (x, y);
        }
    }
}

pub fn ex2() {
    println!("Plotting Cauchy distribution sample");
    plot_cauchy_dist();

    let sample_size = 100_000_000;
    let mut rng = rand::thread_rng();
    println!("\nSampling {} points inside the unit circle", sample_size);
    println!("Analytic method...");
    let mut timer = std::time::Instant::now();
    for _ in 0..sample_size {
        sample_circle_analytic(&mut rng);
    }
    println!("... done in {:.2} sec", timer.elapsed().as_secs_f32());
    println!("Rejection method...");
    timer = std::time::Instant::now();
    for _ in 0..sample_size {
        sample_circle_rejection(&mut rng);
    }
    println!("... done in {:.2} sec", timer.elapsed().as_secs_f32());
}
