pub use F64x2::test_utils::*;

pub fn init_logger() { let _ = env_logger::builder().is_test(true).try_init(); }
