use crate::utils::plot_histogram::plot_histogram;
use ndhistogram::{
    axis::{Axis, UniformNoFlow},
    ndhistogram, AxesTuple, Histogram, VecHistogram,
};

struct MLCG {
    // Multiplicative Linear Congruential Generator
    init: u32,
    mult: u32,
    current: u32,
    step: usize,
}

impl MLCG {
    fn new(init: u32, mult: u32) -> MLCG {
        MLCG {
            init,
            mult,
            current: init,
            step: 0,
        }
    }

    fn reset(&mut self) -> () {
        self.current = self.init;
        self.step = 0;
    }
}

impl Iterator for MLCG {
    type Item = u32;

    fn next(&mut self) -> Option<u32> {
        self.step += 1;
        (self.current, _) = self.current.overflowing_mul(self.mult);
        if self.current != self.init {
            return Some(self.current);
        } else {
            return None;
        }
    }
}

pub fn uniform_random_sampling() {
    let mut timer = std::time::Instant::now();

    let mut mlcg = MLCG::new(987654321, 663608941);
    // for _ in &mut mlcg {} // exhausting the iterator

    println!(
        "Full iteration: {} steps (log2(steps) = {}), took {:.2} seconds",
        mlcg.step,
        (mlcg.step as f32).log2(),
        timer.elapsed().as_secs_f32()
    );

    for (sample_idx, sample_size) in [
        Some(2usize.pow(8)),
        Some(2usize.pow(10)),
        Some(2usize.pow(15)),
        Some(2usize.pow(20)),
        Some(2usize.pow(25)),
        None,
    ]
    .iter()
    .enumerate()
    {
        println!("\nHistograming with sample size = {:?}", sample_size);
        timer = std::time::Instant::now();
        let mut hist = ndhistogram!(UniformNoFlow::new(100, 0.0, u32::MAX as f64));
        mlcg.reset();
        for (idx, i) in (&mut mlcg).enumerate() {
            hist.fill(&(i as f64));
            if let Some(limit) = sample_size {
                if idx == *limit {
                    break;
                }
            }
        }
        println!(
            "Histogram filled in {:.2} sec",
            timer.elapsed().as_secs_f32()
        );
        let chi2_per_dof = uniform_chi2_per_dof(&hist, mlcg.step);
        println!("Histogram chi2 / d.o.f. = {}", chi2_per_dof);

        plot_histogram(
            &hist,
            &format!(
                "Simple MLCG, 2^{:.0} samples, chi2/dof = {:.2}",
                (mlcg.step as f64).log2(),
                chi2_per_dof
            ),
            "u32 value",
            &format!("pics/1-simple-mlcg-{}.png", sample_idx),
        )
        .expect("Failed to plot histogram");
    }
}

fn uniform_chi2_per_dof(
    hist: &VecHistogram<AxesTuple<(UniformNoFlow<f64>,)>, f64>,
    total_samples: usize,
) -> f64 {
    let ax = &hist.axes().as_tuple().0;
    let bin_count = ax.num_bins() as f64;
    let bin_count_mean = total_samples as f64 / bin_count;
    let mut chi2: f64 = 0.0;
    for value in hist.values() {
        // Pearson's chi squared for binned data
        chi2 += (bin_count_mean - value).powi(2) / bin_count_mean;
    }
    return chi2 / bin_count;
}
