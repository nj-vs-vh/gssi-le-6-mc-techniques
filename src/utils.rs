pub fn mean_std(series: &[f32]) -> (f32, f32) {
    let count = series.len() as f32;
    let (sum, square_sum) = series
        .iter()
        .map(|v| (v, v.powi(2)))
        .fold((0.0, 0.0), |(acc, acc_sq), (v, v_sq)| {
            (acc + v, acc_sq + v_sq)
        });
    let mean = sum / count;
    let var = square_sum / count - (sum / count).powi(2);
    (mean, var.sqrt())
}
