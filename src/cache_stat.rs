use log::info;

#[derive(Default)]
pub struct CacheStat {
    pub hit: usize,
    pub miss: usize,
}

impl CacheStat {
    pub fn new() -> Self { Self { hit: 0, miss: 0 } }

    pub fn hit(&mut self) { self.hit += 1; }

    pub fn miss(&mut self) { self.miss += 1; }

    pub fn show(&self, tag: &str) {
        info!(
            "[Cache] tag = {tag}, hit = {}, miss = {}, miss ratio = {:.6e}",
            self.hit,
            self.miss,
            self.miss as f64 / (self.miss + self.hit) as f64,
        );
    }
}
