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

## Ex. 2: Random sampling



<details>
<summary>Source code</summary>

```rust
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

```

</details>




<details>
<summary>Execution log</summary>

```
Plotting Cauchy distribution sample
Histogram has been saved to out/ex2/cauchy.png

Sampling 100000000 points inside the unit circle
Analytic method...
... done in 1.08 sec
Rejection method...
... done in 1.95 sec


Total runtime: 3.05 sec

```

</details>



### Ex. 2.1: Inversion

Inverted CDF for Cauchy distribution is: $C^{-1}(t) = \tan(\pi (t - \frac{1}{2}))$.

Using it, we obtain a sample:

![cauchy-hist](out/ex2/cauchy.png)

### Ex. 2.2: Inversion and rejection

In principle, the result could depend on compiler optimizations. In Rust, optimizations are
controlled by `--release` flag. Both with and without this flag the analytic method works
faster, at least on my machine, but optimization makes the lead more pronounced:

|                             	| Analytic 	| Rejection 	|
|-----------------------------	|----------	|-----------	|
| Debug build (non-optimized) 	| 52.78    	| 65.49     	|
| Release build (optimized)   	| 1.07     	| 1.94      	|


## Ex. 3: Numerical estimation of $\pi$



<details>
<summary>Source code</summary>

```rust
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

```

</details>




<details>
<summary>Execution log</summary>

```
Thrown 5000000 points, recorded 500 pi estimation(s)
Thrown 5000000 points, recorded 500 pi estimation(s)
Thrown 5000000 points, recorded 500 pi estimation(s)
Plot has been saved to out/ex3/pi.png
Sample moments: 3.139188e0 +/- 2.3990864e-1, k = sigma*sqrt(N) = 1.69641
Histogram params for 50 throws: range 1.7274..4.5558, 35 bins
Histogram has been saved to out/ex3/pi-est-distribution-50-throws.png
Sample moments: 3.140801e0 +/- 1.17763884e-1, k = sigma*sqrt(N) = 1.17764
Histogram params for 100 throws: range 2.1416..4.1416, 50 bins
Histogram has been saved to out/ex3/pi-est-distribution-100-throws.png
Sample moments: 3.1418881e0 +/- 4.7691856e-2, k = sigma*sqrt(N) = 1.06642
Histogram params for 500 throws: range 2.6944..3.5888, 111 bins
Histogram has been saved to out/ex3/pi-est-distribution-500-throws.png
Sample moments: 3.1414034e0 +/- 6.4998284e-2, k = sigma*sqrt(N) = 2.05543
Histogram params for 1000 throws: range 2.8254..3.4578, 158 bins
Histogram has been saved to out/ex3/pi-est-distribution-1000-throws.png
Sample moments: 3.1415582e0 +/- 2.1306079e-2, k = sigma*sqrt(N) = 1.50657
Histogram params for 5000 throws: range 3.0002..3.2830, 353 bins
Histogram has been saved to out/ex3/pi-est-distribution-5000-throws.png


Total runtime: 5.93 sec

```

</details>



To check the convergence, I plot the logarithm of absolute estimation error
vs. the number of throws. All traces converge (error decreases), but rather
slowly and with large fluctuations.

![pi-throws](out/ex3/pi.png)


### Ex. 3.1: Uncertainty evaluation

To estimate the error, we can fix the number of throws $N$ and repeat the procedure to get the
distribution of the $\pi$ estimate. From it, we obtain the Gaussian-looking sample of $\pi$ and
calculate the scaling factor $k \equiv \sigma \sqrt{N}$, which turns out to be around 1-2,
although it fluctuates significantly and seems to depend on the number of throws.
So, using the conservative value of $k=2$ we conclude, that to get the precision
$\sigma = 10^{-4}$, we need to throw at least
$N = \left( \frac{k}{\sigma} \right)^2 = 4 \cdot 10^8$ points.

![pi-est-1](out/ex3/pi-est-distribution-50-throws.png)
![pi-est-2](out/ex3/pi-est-distribution-100-throws.png)
![pi-est-3](out/ex3/pi-est-distribution-500-throws.png)
![pi-est-4](out/ex3/pi-est-distribution-1000-throws.png)
![pi-est-5](out/ex3/pi-est-distribution-5000-throws.png)


## Ex. 5: Monte-Carlo integration



<details>
<summary>Source code</summary>

```rust
use itertools::Itertools;
use ndhistogram::{axis::UniformNoFlow, ndhistogram, Histogram};
use rand::prelude::*;

use crate::{
    plot::{plot_histogram, AxLim},
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
        println!("\nIntegrating x^{}", pow);
        let sample_count = 1_000_000;
        let estimations_count = 1_000;
        let estimations = (0..estimations_count)
            .map(|_| onedim_mc_integral(rng, pow, sample_count))
            .collect::<Vec<_>>();
        let (mu, sigma) = mean_std(&estimations);
        println!("Estimate sigma: {}", sigma);
        let k_est = sigma * (sample_count as f32).sqrt();

        let true_value = 1.0 / (1.0 + pow as f64);

        let lo = (mu - 5.0 * sigma) as f64;
        let hi = (mu + 5.0 * sigma) as f64;
        let mut hist = ndhistogram!(UniformNoFlow::new(60, lo, hi));
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
            AxLim::Range(lo, hi),
            Some(vec![true_value]),
        )
        .expect("Failed to plot histogram");
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
            AxLim::enlarged_range(
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

    println!("\nEx. 5.1");
    plot_onedim_integral_mc_distributions(&mut rng);
    println!("\nEx. 5.2");
    plot_ndim_integral_mc_vs_midpoint(&mut rng, Integrand::SquaresSum);
    println!("\nEx. 5.3");
    plot_ndim_integral_mc_vs_midpoint(&mut rng, Integrand::ExpProduct);
}

```

</details>




<details>
<summary>Execution log</summary>

```

Ex. 5.1

Integrating x^1
Estimate sigma: 0.00055939646
Histogram has been saved to out/ex5/1/x^1-mc-integral-1000000-samples.png

Integrating x^2
Estimate sigma: 0.00028628026
Histogram has been saved to out/ex5/1/x^2-mc-integral-1000000-samples.png

Integrating x^3
Estimate sigma: 0.00033982936
Histogram has been saved to out/ex5/1/x^3-mc-integral-1000000-samples.png

Integrating x^4
Estimate sigma: 0.00012207031
Histogram has been saved to out/ex5/1/x^4-mc-integral-1000000-samples.png

Integrating x^5
Estimate sigma: 0.00033707765
Histogram has been saved to out/ex5/1/x^5-mc-integral-1000000-samples.png

Ex. 5.2
Midpoint integration took 0.0016 sec
500 MC integrations took 0.92 sec, on average 0.0018 sec
Histogram has been saved to out/ex5/2/squares-sum-mc-vs-midpoint-1-dim.png
Midpoint integration took 0.0015 sec
500 MC integrations took 1.11 sec, on average 0.0022 sec
Histogram has been saved to out/ex5/2/squares-sum-mc-vs-midpoint-2-dim.png
Midpoint integration took 0.0016 sec
500 MC integrations took 1.31 sec, on average 0.0026 sec
Histogram has been saved to out/ex5/2/squares-sum-mc-vs-midpoint-3-dim.png
Midpoint integration took 0.0017 sec
500 MC integrations took 1.56 sec, on average 0.0031 sec
Histogram has been saved to out/ex5/2/squares-sum-mc-vs-midpoint-4-dim.png
Midpoint integration took 0.0016 sec
500 MC integrations took 1.61 sec, on average 0.0032 sec
Histogram has been saved to out/ex5/2/squares-sum-mc-vs-midpoint-5-dim.png
Midpoint integration took 0.0013 sec
500 MC integrations took 1.42 sec, on average 0.0028 sec
Histogram has been saved to out/ex5/2/squares-sum-mc-vs-midpoint-6-dim.png
Midpoint integration took 0.0022 sec
500 MC integrations took 2.63 sec, on average 0.0053 sec
Histogram has been saved to out/ex5/2/squares-sum-mc-vs-midpoint-7-dim.png
Midpoint integration took 0.0020 sec
500 MC integrations took 2.43 sec, on average 0.0049 sec
Histogram has been saved to out/ex5/2/squares-sum-mc-vs-midpoint-8-dim.png

Ex. 5.3
Midpoint integration took 0.0017 sec
500 MC integrations took 1.00 sec, on average 0.0020 sec
Histogram has been saved to out/ex5/2/exp-product-mc-vs-midpoint-1-dim.png
Midpoint integration took 0.0018 sec
500 MC integrations took 1.29 sec, on average 0.0026 sec
Histogram has been saved to out/ex5/2/exp-product-mc-vs-midpoint-2-dim.png
Midpoint integration took 0.0020 sec
500 MC integrations took 1.56 sec, on average 0.0031 sec
Histogram has been saved to out/ex5/2/exp-product-mc-vs-midpoint-3-dim.png
Midpoint integration took 0.0023 sec
500 MC integrations took 1.89 sec, on average 0.0038 sec
Histogram has been saved to out/ex5/2/exp-product-mc-vs-midpoint-4-dim.png
Midpoint integration took 0.0026 sec
500 MC integrations took 2.06 sec, on average 0.0041 sec
Histogram has been saved to out/ex5/2/exp-product-mc-vs-midpoint-5-dim.png
Midpoint integration took 0.0020 sec
500 MC integrations took 1.83 sec, on average 0.0037 sec
Histogram has been saved to out/ex5/2/exp-product-mc-vs-midpoint-6-dim.png
Midpoint integration took 0.0034 sec
500 MC integrations took 3.41 sec, on average 0.0068 sec
Histogram has been saved to out/ex5/2/exp-product-mc-vs-midpoint-7-dim.png
Midpoint integration took 0.0031 sec
500 MC integrations took 3.12 sec, on average 0.0062 sec
Histogram has been saved to out/ex5/2/exp-product-mc-vs-midpoint-8-dim.png


Total runtime: 63.04 sec

```

</details>



### Ex. 5.1: Unidimensional integration

To estimate the number of samples required for a given precision, we can use the same
dependence of estimation's variance $\sigma \propto N^{-1/2}$. In the previous exercise
the coefficient $k$ in this dependence was shown to be of order $1$. Using this rough
value, we obtain the number of samples required for precision $\sigma = 0.001$:
$N = \sigma^{-2} = 10^6$.

From the samples we obtain the actual value $k \approx 0.2$, so the actual error is
factor of $5$ smaller than the required value.

![x-1-mc-integral](out/ex5/1/x^1-mc-integral-1000000-samples.png)
![x-2-mc-integral](out/ex5/2/x^1-mc-integral-1000000-samples.png)
![x-3-mc-integral](out/ex5/3/x^1-mc-integral-1000000-samples.png)
![x-4-mc-integral](out/ex5/4/x^1-mc-integral-1000000-samples.png)
![x-5-mc-integral](out/ex5/5/x^1-mc-integral-1000000-samples.png)

### Ex. 5.2: Multidimensional integration

For each dimension I produce 500 MC integrations and compare them with the midpoint
formula and the true integral value. Since the number of samples for MC and cells
for midpoint formula is chosen to be equal, evaluation time of the methods is within
an order of magnitude from each other, with MC integration taking $\approx 2$ times more
time for large numbers of dimensions.

On the plots below, blue line is the true integral value, orange line is midpoint approximation
result, yellow histogram is MC integration results distribution.

For low number of dimensions the MC integration introduces additional error:

![mc-vs-midpoint-low-dim](out/ex5/2/squares-sum-mc-vs-midpoint-3-dim.png)

For higher dimensionality, the midpoint approximation introduces significant bias,
while MC integration stays consistent.

![mc-vs-midpoint-high-dim](out/ex5/2/squares-sum-mc-vs-midpoint-7-dim.png)

<details>

<summary>Plots for all dimensions</summary>


![mc-vs-midpoint-1-dim](out/ex5/2/squares-sum-mc-vs-midpoint-1-dim.png)
![mc-vs-midpoint-2-dim](out/ex5/2/squares-sum-mc-vs-midpoint-2-dim.png)
![mc-vs-midpoint-3-dim](out/ex5/2/squares-sum-mc-vs-midpoint-3-dim.png)
![mc-vs-midpoint-4-dim](out/ex5/2/squares-sum-mc-vs-midpoint-4-dim.png)
![mc-vs-midpoint-5-dim](out/ex5/2/squares-sum-mc-vs-midpoint-5-dim.png)
![mc-vs-midpoint-6-dim](out/ex5/2/squares-sum-mc-vs-midpoint-6-dim.png)
![mc-vs-midpoint-7-dim](out/ex5/2/squares-sum-mc-vs-midpoint-7-dim.png)
![mc-vs-midpoint-8-dim](out/ex5/2/squares-sum-mc-vs-midpoint-8-dim.png)

</details>


### Ex. 5.3: Multidimensional integration for product of exponents

Qualitatively the results are the same as for previous exercise, but the "critical"
number of dimensions, at which midpoint integration becomes worse than MC, is higher.
For example, the previous function in 5 dimensions would be better integrated by MC,
while for this one the errors of two methods are comparable.

The explaination is as follows, considering one dimension for clarity:
on $x \in [0, 1]$ the previous function $x^2$ varies in $[0, 1]$, while the exponent
varies only in $[0, e^{-1} \approx 0.37]$. So, the previous function is "cuspier",
making midpoint formula more biased as it doesn't probe the peak and underestimates
the integral.

<details>

<summary>Plots for all dimensions</summary>


![mc-vs-midpoint-1-dim-exp](out/ex5/2/exp-product-mc-vs-midpoint-1-dim.png)
![mc-vs-midpoint-2-dim-exp](out/ex5/2/exp-product-mc-vs-midpoint-2-dim.png)
![mc-vs-midpoint-3-dim-exp](out/ex5/2/exp-product-mc-vs-midpoint-3-dim.png)
![mc-vs-midpoint-4-dim-exp](out/ex5/2/exp-product-mc-vs-midpoint-4-dim.png)
![mc-vs-midpoint-5-dim-exp](out/ex5/2/exp-product-mc-vs-midpoint-5-dim.png)
![mc-vs-midpoint-6-dim-exp](out/ex5/2/exp-product-mc-vs-midpoint-6-dim.png)
![mc-vs-midpoint-7-dim-exp](out/ex5/2/exp-product-mc-vs-midpoint-7-dim.png)
![mc-vs-midpoint-8-dim-exp](out/ex5/2/exp-product-mc-vs-midpoint-8-dim.png)

</details>


