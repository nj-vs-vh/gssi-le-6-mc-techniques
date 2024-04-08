use ndhistogram::{axis::UniformNoFlow, AxesTuple, Histogram as Hist};
use plotters::prelude::*;

pub fn plot_histogram(
    hist: &dyn Hist<AxesTuple<(UniformNoFlow,)>, f64>,
    title: &str,
    x_caption: &str,
    filename: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let hist_axis = hist.axes().as_tuple().0.to_owned();
    let hist_max_count = hist.values().fold(f64::NEG_INFINITY, |a, &b| a.max(b));

    let fig = BitMapBackend::new(filename, (640, 480)).into_drawing_area();

    fig.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&fig)
        .margin(5)
        .y_label_area_size(60)
        .x_label_area_size(40)
        .caption(title, ("sans-serif", 20))
        .build_cartesian_2d(
            hist_axis.low().to_owned()..hist_axis.high().to_owned(),
            0.0..(hist_max_count * 1.05),
        )?;

    chart
        .configure_mesh()
        // .bold_line_style(WHITE.mix(0.3))
        .y_desc("Values per bucket, #")
        .x_desc(x_caption)
        .axis_desc_style(("sans-serif", 14))
        .draw()?;

    chart
        .draw_series(AreaSeries::new(
            hist.iter()
                .filter(|item| item.bin.start().is_some() && item.bin.end().is_some())
                .flat_map(|item| {
                    let value = item.value.to_owned();
                    return [
                        (item.bin.start().unwrap(), value),
                        (item.bin.end().unwrap(), value),
                    ];
                }),
            0.0,
            RED.mix(0.6),
        ))
        .unwrap();

    fig.present()?;
    println!("Histogram has been saved to {}", filename);
    Ok(())
}
