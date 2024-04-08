use std::f64::consts::PI;

use ndhistogram::{axis::UniformNoFlow, ndhistogram, Histogram};
use rand::prelude::*;

use crate::utils::plot_histogram::{plot_histogram, XLim};

fn sample_cauchy_dist() {
    let mut rng = rand::thread_rng();
    let sample_size = 10000;
    let sample: Vec<f64> = (0..sample_size)
        // inverse CDF C^-1(t) = tan(pi*(t - 1/2))
        .map(|_| (PI * (rng.gen::<f64>() - 0.5)).tan())
        .collect();
    let mut hist = ndhistogram!(UniformNoFlow::new(100, -10.0, 10.0));
    sample.iter().map(|s| hist.fill(s)).count();
    plot_histogram(
        &hist,
        "Cauchy distribution sample",
        "x",
        "out/ex2/cauchy.png",
        XLim::Range(-10.0, 10.0),
    )
    .expect("Failed to plot histogram");
}

pub fn ex2() {
    sample_cauchy_dist();
}
