use clap::Parser;

mod ex1;
mod ex2;
mod ex3;
mod plot;

#[derive(Parser, Debug)]
#[command(name = "AoC 2022 solutions")]
#[command(author = "Igor V. <gosha.vaiman@gmail.com>")]
#[command(version = "1.3.1.2")]
struct CliArgs {
    exercise: String,
}

fn main() {
    let args = CliArgs::parse();
    match args.exercise.as_str() {
        "1" => ex1::ex1(),
        "2" => ex2::ex2(),
        "3" => ex3::ex3(),
        _ => {
            println!("Exercise {:?} not found", args.exercise);
        }
    }
}
