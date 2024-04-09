use num::{cast, Float};

fn sum_inverse<F: Float>(n: usize) -> F {
    let n_float: F = cast(n).unwrap();
    let one: F = cast(1.0).unwrap();
    (0..n)
        .map(|_| one.div(n_float))
        .fold(cast(0.0).unwrap(), |acc, v| acc + v)
}

pub fn ex6() {
    println!(
        "{:>8} | {:<12} | {:<12} | {:<12}",
        "n", "half (f16)", "single (f32)", "double (f64)",
    );
    for n in [100, 1000, 10_000, 100_000, 1_000_000] {
        println!(
            "{:>8} | {:<12} | {:<12} | {:<12}",
            n,
            sum_inverse::<half::f16>(n),
            sum_inverse::<f32>(n),
            sum_inverse::<f64>(n),
        );
    }
}
