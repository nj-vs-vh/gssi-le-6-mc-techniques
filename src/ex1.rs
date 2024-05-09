use crate::plot::{plot_histogram, AxLim};
use ndhistogram::{
    axis::{Axis, UniformNoFlow},
    ndhistogram, AxesTuple, Histogram, VecHistogram,
};

use num::{cast, traits::WrappingMul, Integer, NumCast};
use std::fmt::Display;

enum Modulus<T> {
    TypeMax,
    Value(T),
}

struct MCG<T: Integer + Display + WrappingMul + Copy> {
    // Multiplicative Congruential Generator
    init: T,
    mult: T,
    modulus: Modulus<T>,

    // generator state
    current: T,
    step: usize,
}
impl<T: Integer + Display + WrappingMul + Copy> Iterator for MCG<T> {
    type Item = T;

    fn next(&mut self) -> Option<T> {
        self.step += 1;
        self.current = match self.modulus {
            Modulus::TypeMax => self.current.wrapping_mul(&self.mult),
            Modulus::Value(m) => (self.current * self.mult).rem(m),
        };

        if !self.current.is_zero() && self.current != self.init {
            return Some(self.current);
        } else {
            return None;
        }
    }
}

impl<T: Integer + Display + WrappingMul + Copy> std::fmt::Display for MCG<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&format!(
            "MCG {{ init: {}, mult: {} }}",
            self.init, self.mult
        ))
    }
}

impl<T: Integer + Display + WrappingMul + NumCast + Copy> MCG<T> {
    fn new(init: T, mult: T, modulus: Modulus<T>) -> MCG<T> {
        MCG {
            init,
            mult,
            modulus,
            current: init,
            step: 0,
        }
    }
    fn new_modulo_type_max(init: T, mult: T) -> MCG<T> {
        MCG::new(init, mult, Modulus::TypeMax)
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
        println!("\nTesting {}: {}", name, &self);
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
            let mut hist = ndhistogram!(UniformNoFlow::new(100, 0.0, 1.0));
            let scale = 1.0 / (u32::MAX as f64 - 1.0);
            self.reset();
            for (idx, i) in (&mut self).enumerate() {
                hist.fill(&(scale * cast::<T, f64>(i).unwrap()));
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
                "standard uniform value",
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
    MCG::<u32>::new_modulo_type_max(987654321, 663608941).run_tests("from exercise");
    MCG::<u32>::new_modulo_type_max(10, 9).run_tests("poorly selected params");

    println!("\n\nEx 1.1: MINSTD algorithm");
    let m = 2u64.pow(31) - 1;
    let m_f64 = m as f64;
    for seed in [1, 7u64.pow(2), 2u64.pow(4)] {
        let mut minstd = MCG::<u64>::new(seed, 16807, Modulus::Value(m));
        let start = std::time::Instant::now();
        let period = minstd.calculate_period();
        println!(
            "Seed {:<10} -> period {:>9} ({:.2}% of range, calculated in {:.2} sec)",
            seed,
            period,
            100.0 * (period as f64) / (m_f64),
            start.elapsed().as_secs_f32(),
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
