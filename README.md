<!-- DO NOT EDIT, GENERATED AUTOMATICALLY -->

# LE-6 Monte Carlo techniques: Excercises

*Igor Vaiman*

## Introduction

The exercises are done in Rust programming language. It's not widely used in physics
community but hopefully can see wider adoption in the future. In combines access to
low-level performance and fine control with high-level features of modern programming
languages. Its emphasis on memory safety, "fearless concurrency" and compile-time error
detection might be helpful in parallelizing scientific tasks.

## Ex. 1: Uniform random sampling



<details>
<summary>Source code</summary>

```rust
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

```

</details>




<details>
<summary>Execution log</summary>

```

Ex 1: Uniform random sampling

Testing from exercise: MCG { init: 987654321, mult: 663608941 }
Full iteration: 1073741824 steps (log2(steps) = 30), took 1.03 seconds

Histograming with sample size = Some(1024)
Histogram filled in 0.00 sec
Histogram chi2 / d.o.f. = 0.7402439024390247
Histogram has been saved to out/ex1/from exercise-2^10.png

Histograming with sample size = Some(32768)
Histogram filled in 0.00 sec
Histogram chi2 / d.o.f. = 1.0199697885196377
Histogram has been saved to out/ex1/from exercise-2^15.png

Histograming with sample size = Some(1048576)
Histogram filled in 0.00 sec
Histogram chi2 / d.o.f. = 1.3205646414140306
Histogram has been saved to out/ex1/from exercise-2^20.png

Histograming with sample size = Some(33554432)
Histogram filled in 0.09 sec
Histogram chi2 / d.o.f. = 0.8355132125165103
Histogram has been saved to out/ex1/from exercise-2^25.png

Histograming with sample size = None
Histogram filled in 3.01 sec
Histogram chi2 / d.o.f. = 0.000000016503036022602618
Histogram has been saved to out/ex1/from exercise-2^30.png

Testing poorly selected params: MCG { init: 10, mult: 9 }
Full iteration: 268435456 steps (log2(steps) = 28), took 0.25 seconds

Histograming with sample size = Some(1024)
Histogram filled in 0.00 sec
Histogram chi2 / d.o.f. = 1.0192682926829275
Histogram has been saved to out/ex1/poorly selected params-2^10.png

Histograming with sample size = Some(32768)
Histogram filled in 0.00 sec
Histogram chi2 / d.o.f. = 1.2246144221672923
Histogram has been saved to out/ex1/poorly selected params-2^15.png

Histograming with sample size = Some(1048576)
Histogram filled in 0.00 sec
Histogram chi2 / d.o.f. = 0.833544613318812
Histogram has been saved to out/ex1/poorly selected params-2^20.png

Histograming with sample size = Some(33554432)
Histogram filled in 0.09 sec
Histogram chi2 / d.o.f. = 0.7983531150712634
Histogram has been saved to out/ex1/poorly selected params-2^25.png

Histograming with sample size = None
Histogram filled in 0.75 sec
Histogram chi2 / d.o.f. = 0.00000009968876838725717
Histogram has been saved to out/ex1/poorly selected params-2^28.png


Ex 1.1: MINSTD algorithm
Seed 1          -> period 2147483646 (100.00% of range, calculated in 8.80 sec)
Seed 49         -> period 2147483646 (100.00% of range, calculated in 8.80 sec)
Seed 16         -> period 2147483646 (100.00% of range, calculated in 8.80 sec)


Total runtime: 31.67 sec

```

</details>



The exercise is to write a Multiplicative Congruent Generator (MCG) that produces
a stream of numbers according to the rule $X_{i+1} = (M \cdot X_{i}) \text{mod} N$.
$N$ is fixed to $2^{32}$ by the fact that we're dealing with unsigned 32-bit integers
(`u32` Rust datatype) and use overflowing multiplication that wraps numbers on overflow.
The MCG then has two parameters: a multiplier $M$ and an initial value $X_0$.

Results for MCG with initial value: 987654321 and multiplier: 663608941:
```
Full iteration: 1073741824 steps (log2(steps) = 30), took 1.03 seconds
```

Samples from this MCG with different sizes and with $\chi^2 / \text{d.o.f.}$ test
of uniformity:

![2-to-the-10th](out/ex1/from%20exercise-2^10.png)
![2-to-the-15th](out/ex1/from%20exercise-2^15.png)
![2-to-the-20th](out/ex1/from%20exercise-2^20.png)
![2-to-the-25th](out/ex1/from%20exercise-2^25.png)
![2-to-the-30th](out/ex1/from%20exercise-2^30.png)


I decided to try other, rather small and not specifically selected parameters: 
initial value: 10 and multiplier: 9. The cycle size is 4 times shorter, $2^{28}$, but the
distributions look fairly flat anyway. Perhaps, the sample is much more correlated
in this case.

<details>
<summary>Samples</summary>

![1](out/ex1/poorly%20selected%20params-2^10.png)
![2](out/ex1/poorly%20selected%20params-2^15.png)
![3](out/ex1/poorly%20selected%20params-2^20.png)
![4](out/ex1/poorly%20selected%20params-2^25.png)
![5](out/ex1/poorly%20selected%20params-2^28.png)

</details>

### Ex. 1.1: MINSTD algorithm

Again, MINST is an MCG with multiplier $M = 7^5 = 16807$, but $N$ is not a Mersenne prime $2^{31} - 1$.
This requires one to use 64-bit integer type and perform modulus operation explicitly. But for that we
get a longer sequence, spanning the full range of $2^{31} - 1$, and this seems to happen irrespective
of the seed (initial value).
