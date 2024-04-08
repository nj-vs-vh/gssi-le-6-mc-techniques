use clap::Parser;

mod ex1;

#[derive(Parser, Debug)]
#[command(name = "AoC 2022 solutions")]
#[command(author = "Igor V. <gosha.vaiman@gmail.com>")]
#[command(version = "1.3.1.2")]
struct CliArgs {
    exercise: String,
}

fn main() {
    let args = CliArgs::parse();
    println!("Running ex. {}\n", args.exercise);

    match args.exercise.as_str() {
        "1" => ex1::uniform_random_sampling(),
        _ => {
            println!("Exercise not found");
        }
    }
}
