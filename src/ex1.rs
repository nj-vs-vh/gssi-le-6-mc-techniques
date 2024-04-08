struct MLCG {
    // Multiplicative Linear Congruential Generator
    init: u32,
    mult: u32,
    current: u32,
    step: u32,
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
}

impl Iterator for MLCG {
    type Item = u32;

    fn next(&mut self) -> Option<u32> {
        self.step += 1;
        (self.current, _) = self.current.overflowing_mul(self.mult);
        if self.current != self.init {
            return Some(self.current);
        } else {
            return None;
        }
    }
}

pub fn uniform_random_sampling() {
    let now = std::time::Instant::now();

    let mut mlcg = MLCG::new(987654321, 663608941);
    for _ in &mut mlcg {} // exhausting the iterator

    println!(
        "Full iteration was {} steps (log_2(steps) = {}) and took {:.2} seconds",
        mlcg.step,
        (mlcg.step as f32).log2(),
        now.elapsed().as_secs_f32()
    );
}
