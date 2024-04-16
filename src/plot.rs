use ndhistogram::{axis::UniformNoFlow, AxesTuple, Histogram as Hist};
use plotters::prelude::*;

pub enum AxLim {
    FromData,
    Range(f64, f64),
}

impl AxLim {
    pub fn enlarged_range(min: f64, max: f64, enlarge: f64) -> AxLim {
        let range = max - min;
        AxLim::Range(min - range * enlarge, max + range * enlarge)
    }
}

pub fn plot_histogram(
    hist: &dyn Hist<AxesTuple<(UniformNoFlow,)>, f64>,
    title: &str,
    x_caption: &str,
    filename: &str,
    xlim: AxLim,
    vertical_lines: Option<Vec<f64>>,
) -> Result<(), Box<dyn std::error::Error>> {
    let hist_nonempty_bins = hist
        .iter()
        .filter(|item| *item.value > 0.0)
        .map(|item| item.bin)
        .collect::<Vec<_>>();
    if hist_nonempty_bins.is_empty() {
        return Err("Hist is empty".into());
    }
    let hist_max_count = hist.values().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    let y_min = 0.0;
    let y_max = hist_max_count * 1.05;

    let fig = BitMapBackend::new(filename, (640, 480)).into_drawing_area();

    fig.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&fig)
        .margin(5)
        .y_label_area_size(60)
        .x_label_area_size(40)
        .caption(title, ("sans-serif", 20))
        .build_cartesian_2d(
            match xlim {
                AxLim::FromData => {
                    hist_nonempty_bins.first().unwrap().start().unwrap()
                        ..hist_nonempty_bins.last().unwrap().start().unwrap()
                }
                AxLim::Range(low, high) => low..high,
            },
            y_min..y_max,
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
            Palette9999::pick(0).mix(0.5),
        ))
        .unwrap();
    if let Some(values) = vertical_lines {
        for (idx, value) in values.into_iter().enumerate() {
            chart
                .draw_series(LineSeries::new(
                    [(value, y_min), (value, y_max)],
                    Palette9999::pick(idx + 1).stroke_width(2),
                ))
                .unwrap();
        }
    }

    fig.present()?;
    println!("Histogram has been saved to {}", filename);
    Ok(())
}

type Bounds = (f32, f32);

const DEFAULT_BOUNDS: Bounds = (f32::INFINITY, f32::NEG_INFINITY);

fn find_extent<'a>(
    coords: impl Iterator<Item = &'a f32>,
    enlarge_left: f32,
    enlarge_right: f32,
) -> Option<(f32, f32)> {
    let (min, max) = coords.fold(DEFAULT_BOUNDS, |acc, &item| {
        (item.min(acc.0), item.max(acc.1))
    });
    if min.is_finite() && max.is_finite() {
        let range = max - min;
        Some((min - range * enlarge_left, max + range * enlarge_right))
    } else {
        None
    }
}

pub struct Line {
    pub data: Vec<(f32, f32)>,
    pub label: String,
}

pub fn plot_lines(
    lines: Vec<Line>,
    title: &str,
    x_caption: &str,
    y_caption: &str,
    filename: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let fig = BitMapBackend::new(filename, (640, 480)).into_drawing_area();

    let xlim = find_extent(
        lines
            .iter()
            .flat_map(|l| l.data.iter().map(|c| &c.0))
            .filter(|v| v.is_finite()),
        0.1,
        0.1,
    )
    .ok_or("Failed to compute x limits")?;
    let ylim = find_extent(
        lines
            .iter()
            .flat_map(|l| l.data.iter().map(|c| &c.1))
            .filter(|v| v.is_finite()),
        0.1,
        0.1,
    )
    .ok_or("Failed to compute y limits")?;

    fig.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&fig)
        .margin(5)
        .y_label_area_size(60)
        .x_label_area_size(40)
        .caption(title, ("sans-serif", 20))
        .build_cartesian_2d(xlim.0..xlim.1, ylim.0..ylim.1)?;

    chart
        .configure_mesh()
        // .bold_line_style(WHITE.mix(0.3))
        .y_desc(y_caption)
        .x_desc(x_caption)
        .axis_desc_style(("sans-serif", 14))
        .draw()?;

    for (idx, line) in lines.into_iter().enumerate() {
        chart
            .draw_series(LineSeries::new(line.data, Palette99::pick(idx)))
            .unwrap()
            .label(line.label);
    }

    fig.present()?;
    println!("Plot has been saved to {}", filename);
    Ok(())
}

type Line3D = Vec<(f32, f32, f32)>;

pub fn plot_lines_3d(
    lines: Vec<Line3D>,
    title: &str,
    limits: (AxLim, AxLim, AxLim),
    filename: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let area = BitMapBackend::new(filename, (1024, 760)).into_drawing_area();

    area.fill(&WHITE)?;

    let xlim = match limits.0 {
        AxLim::FromData => find_extent(
            lines
                .iter()
                .flat_map(|l| l.iter().map(|c| &c.0))
                .filter(|v| v.is_finite()),
            0.1,
            0.1,
        )
        .ok_or("Failed to compute x limits")?,
        AxLim::Range(min, max) => (min as f32, max as f32),
    };
    let ylim = match limits.1 {
        AxLim::FromData => find_extent(
            lines
                .iter()
                .flat_map(|l| l.iter().map(|c| &c.1))
                .filter(|v| v.is_finite()),
            0.1,
            0.1,
        )
        .ok_or("Failed to compute y limits")?,
        AxLim::Range(min, max) => (min as f32, max as f32),
    };
    let zlim = match limits.2 {
        AxLim::FromData => find_extent(
            lines
                .iter()
                .flat_map(|l| l.iter().map(|c| &c.2))
                .filter(|v| v.is_finite()),
            0.1,
            0.1,
        )
        .ok_or("Failed to compute z limits")?,
        AxLim::Range(min, max) => (min as f32, max as f32),
    };

    let mut chart = ChartBuilder::on(&area)
        .caption(title.to_string(), ("sans", 20))
        .build_cartesian_3d(xlim.0..xlim.1, ylim.0..ylim.1, zlim.0..zlim.1)?;

    chart.with_projection(|mut pb| {
        pb.yaw = 0.5;
        pb.pitch = 0.3;
        pb.scale = 0.9;
        pb.into_matrix()
    });

    chart
        .configure_axes()
        .light_grid_style(BLACK.mix(0.15))
        .max_light_lines(3)
        .draw()?;

    for (idx, line) in lines.iter().enumerate() {
        chart.draw_series(LineSeries::new(
            line.iter().map(|coords| coords.to_owned()),
            &Palette9999::pick(idx),
        ))?;
    }

    // To avoid the IO failure being ignored silently, we manually call the present function
    area.present()
        .expect("Unable to write result to file, please make sure output dir exists");
    println!("3D line plot has been saved to {}", filename);
    Ok(())
}
