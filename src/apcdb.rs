use crate::traits::MyReal;
use byteorder::{BigEndian, ReadBytesExt, WriteBytesExt};
use log::{debug, info};
use std::{
    io::{BufReader, BufWriter},
    path::Path,
};
use F64x2::f64x2;

#[allow(non_snake_case)]
pub fn write_APCDB(
    roots: &[f64x2],
    file_path: impl AsRef<Path>,
    n0: usize,
    n1: usize,
) -> Result<(), std::io::Error> {
    let file_path = file_path.as_ref();
    info!("Writing zeta zeros to {}", file_path.display());
    let file = std::fs::File::create(file_path)?;
    let mut buf_writer = BufWriter::new(file);
    buf_writer.write_u64::<BigEndian>(n0 as u64)?;
    buf_writer.write_u64::<BigEndian>(n1 as u64)?;
    assert!(roots.len() == n1 - n0);
    for &root in roots {
        buf_writer.write_f64::<BigEndian>(root.hi)?;
        buf_writer.write_f64::<BigEndian>(root.lo)?;
    }

    Ok(())
}

#[allow(non_snake_case)]
pub fn read_APCDB<T: MyReal>(
    data_path: impl AsRef<Path>,
    limit: f64,
    mut f: impl FnMut(T),
) -> Result<(), std::io::Error> {
    let file_path = data_path.as_ref();
    info!("Writing zeta zeros to {}", file_path.display());
    let file = std::fs::File::open(file_path)?;
    let mut buf_reader = BufReader::new(file);
    let n0 = buf_reader.read_u64::<BigEndian>()?;
    let n1 = buf_reader.read_u64::<BigEndian>()?;
    for _ in n0..n1 {
        let hi = buf_reader.read_f64::<BigEndian>()?;
        let lo = buf_reader.read_f64::<BigEndian>()?;
        let z = T::zero() + hi + lo;
        if z.fp() > limit {
            break;
        }
        f(z);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{read_APCDB, write_APCDB};
    use crate::lmfdb::read_LMFDB;
    use crate::zeta_zeros::{try_isolate, HybridPrecHardyZ};

    use log::debug;
    use F64x2::{f64x2, test_utils::assert_close};

    fn all_zeta_zeros_upto(t: f64) -> Vec<f64x2> {
        // well, technically g_{-1} = 0 but I was too lazy to do that...
        let mut roots = vec![f64x2::new(14.134725141734695, -8.407109157053214e-16)];

        const PI: f64 = std::f64::consts::PI;
        let mut hardy_z = HybridPrecHardyZ::<f64x2>::new(t, 10, 1e-18);
        let n0 = 0;
        let n1 = t / 2.0 / PI * (t / 2.0 / PI).ln()
            - (0.112 * t.ln() + 0.278 * t.ln().ln() + 3.385 + 0.2 / t);
        let n1 = n1 as usize + 200;
        let stats = try_isolate(&mut hardy_z, n0, n1, 1e-18, 1e-30);
        assert!(stats.height >= t);
        roots.extend(stats.roots.iter().copied().filter(|&x| x.fp() <= t));
        let n_calls_separate = stats.count_separate;
        let n_calls_locate = stats.count_locate;
        let n_zeros = stats.roots.len();
        debug!("{} zeros located. The next call should starts with n0 = {}", n_zeros, n0 + n_zeros);
        debug!(
            "{:.3} calls to separate, {:.3} calls to locate, total = {:.3}",
            n_calls_separate as f64 / n_zeros as f64,
            n_calls_locate as f64 / n_zeros as f64,
            (n_calls_separate + n_calls_locate) as f64 / n_zeros as f64,
        );

        let r1 = stats.roots[0];
        let rn = stats.roots[n_zeros - 1];
        debug!("the first zero is {}, and the last zero is {}", r1, rn);

        roots
    }

    #[test]
    fn check_write_apcdb() {
        crate::init();
        let roots = all_zeta_zeros_upto(1e4);
        let _ = write_APCDB(&roots, "data/apcdb/height_100000.dat", 0, roots.len());
    }

    #[test]
    fn check_read_apcdb() {
        crate::init();
        let mut apcdb_zeros = vec![];
        let _ = read_APCDB("data/apcdb/height_100000.dat", 1e4, |x: f64x2| apcdb_zeros.push(x));

        let mut lmfdb_zeros = vec![];
        let _ = read_LMFDB("data/lmfdb", 1e4, |x: f64x2| lmfdb_zeros.push(x));

        let n1 = apcdb_zeros.len();
        let n2 = lmfdb_zeros.len();
        debug!("{n1} zeros from APCDB and {n2} zeros from LMFDB");
        for i in 0..n1.min(n2) {
            let (r1, r2) = (apcdb_zeros[i], lmfdb_zeros[i]);
            assert_close(r1, r2, 1e-18, 1e-18);
        }
    }
}
