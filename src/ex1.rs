use crate::plot::{plot_histogram, AxLim};
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
impl Iterator for MLCG {
    type Item = u32;

    fn next(&mut self) -> Option<u32> {
        self.step += 1;
        (self.current, _) = self.current.overflowing_mul(self.mult);
        if self.current != 0 && self.current != self.init {
            return Some(self.current);
        } else {
            return None;
        }
    }
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

    fn calculate_period(&mut self) -> usize {
        self.reset();
        for _ in self.into_iter() {} // exhausting the iterator
        let period = self.step.to_owned();
        self.reset();
        return period;
    }

    fn run_tests(mut self, name: &str) -> () {
        println!("\nTesting {}", name);
        self.reset();

        let mut timer = std::time::Instant::now();
        let period = self.calculate_period();
        println!(
            "Full iteration: {} steps (log2(steps) = {}), took {:.2} seconds",
            period,
            (period as f64).log2(),
            timer.elapsed().as_secs_f32()
        );

        for sample_size in [
            Some(2usize.pow(10)),
            Some(2usize.pow(15)),
            Some(2usize.pow(20)),
            Some(2usize.pow(25)),
            None,
        ] {
            println!("\nHistograming with sample size = {:?}", sample_size);
            timer = std::time::Instant::now();
            let mut hist = ndhistogram!(UniformNoFlow::new(100, 0.0, u32::MAX as f64));
            self.reset();
            for (idx, i) in (&mut self).enumerate() {
                hist.fill(&(i as f64));
                if let Some(limit) = sample_size {
                    if idx == limit {
                        break;
                    }
                }
            }
            println!(
                "Histogram filled in {:.2} sec",
                timer.elapsed().as_secs_f32()
            );
            let chi2_per_dof = uniform_chi2_per_dof(&hist, self.step);
            println!("Histogram chi2 / d.o.f. = {}", chi2_per_dof);

            let log2_size = (self.step as f64).log2();
            plot_histogram(
                &hist,
                &format!(
                    "{}, 2^{:.0} samples, chi2/dof = {:.2}",
                    name, log2_size, chi2_per_dof
                ),
                "u32 value",
                &format!("out/ex1/{}-2^{:.0}.png", name, log2_size),
                AxLim::FromData,
                None,
            )
            .expect("Failed to plot histogram");
        }
    }
}

pub fn ex1() {
    println!("\nEx 1: Uniform random sampling");
    MLCG::new(987654321, 663608941).run_tests("Simple MLCG");
    MLCG::new(10, 9).run_tests("Bad MLCG");

    println!("\nEx 1.1: MINSTD algorithm");
    for seed in [
        1,
        2,
        3,
        5,
        16,
        100,
        128,
        1025,
        7u32.pow(5),
        7u32.pow(5) - 1,
        29147598,
    ] {
        let mut minstd = MLCG::new(seed, 16807);
        let period = minstd.calculate_period();
        println!(
            "Seed {:<10} -> period {:>9} ({:.2}% of u32 range)",
            seed,
            period,
            100.0 * (period as f64) / (u32::MAX as f64)
        );
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
