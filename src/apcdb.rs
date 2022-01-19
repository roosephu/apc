use crate::traits::MyReal;
use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use log::{debug, info};
use std::fs::File;
use std::{
    io::{BufRead, BufReader, BufWriter},
    marker::PhantomData,
    path::{Path, PathBuf},
};
use F64x2::f64x2;

#[allow(non_snake_case)]
pub fn write_APCDB(
    roots: &[f64x2],
    file_path: impl AsRef<Path>,
    n0: usize,
    n1: usize,
    atol: f64,
    rtol: f64,
) -> Result<(), std::io::Error> {
    let file_path = file_path.as_ref();
    info!("Writing zeta zeros to {}", file_path.display());
    let file = std::fs::File::create(file_path)?;
    let mut buf_writer = BufWriter::new(file);
    buf_writer.write_u64::<LittleEndian>(n0 as u64)?;
    buf_writer.write_u64::<LittleEndian>(n1 as u64)?;
    buf_writer.write_f64::<LittleEndian>(atol)?;
    buf_writer.write_f64::<LittleEndian>(rtol)?;
    assert!(roots.len() == n1 - n0);
    for &root in roots {
        buf_writer.write_f64::<LittleEndian>(root.hi)?;
        buf_writer.write_f64::<LittleEndian>(root.lo)?;
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
    info!("Reading zeta zeros from {}", file_path.display());
    let file = std::fs::File::open(file_path)?;
    let mut buf_reader = BufReader::new(file);
    let n0 = buf_reader.read_u64::<LittleEndian>()?;
    let n1 = buf_reader.read_u64::<LittleEndian>()?;
    let _atol = buf_reader.read_f64::<LittleEndian>()?;
    let _rtol = buf_reader.read_f64::<LittleEndian>()?;
    for _ in n0..n1 {
        let hi = buf_reader.read_f64::<LittleEndian>()?;
        let lo = buf_reader.read_f64::<LittleEndian>()?;
        let z = T::zero() + hi + lo;
        if z.fp() > limit {
            break;
        }
        f(z);
    }
    Ok(())
}

pub struct APCDB<T> {
    reader: BufReader<File>,
    n_zeros: u64,
    atol: f64,
    rtol: f64,
    _marker: PhantomData<T>,
}

impl<T: MyReal> APCDB<T> {
    pub fn new(path: impl AsRef<Path>) -> std::io::Result<Self> {
        let path = path.as_ref();
        let file = std::fs::File::open(&path)?;
        let mut reader = BufReader::new(file);
        debug!("[APCDB] loading {}", path.display());

        let n0 = reader.read_u64::<LittleEndian>()?;
        let n1 = reader.read_u64::<LittleEndian>()?;
        let atol = reader.read_f64::<LittleEndian>()?;
        let rtol = reader.read_f64::<LittleEndian>()?;

        Ok(Self { reader, n_zeros: n1 - n0, atol, rtol, _marker: PhantomData })
    }

    pub fn directory(path: impl AsRef<Path>) -> impl Iterator<Item = (T, f64)> {
        let path = path.as_ref().to_path_buf();
        let file = std::fs::File::open(path.join("md5.txt")).unwrap();
        let reader = std::io::BufReader::new(file);
        let mut indices = vec![];
        for line in reader.lines() {
            let line = line.unwrap();
            let index = &line[41..line.len() - 4];
            let index = index.parse::<i64>().unwrap();
            indices.push(index);
        }

        indices
            .into_iter()
            .flat_map(move |index| Self::new(path.join(format!("height_{index}.dat"))).unwrap())
    }
}

impl<T: MyReal> Iterator for APCDB<T> {
    type Item = (T, f64);

    fn next(&mut self) -> Option<Self::Item> {
        let reader = &mut self.reader;
        let mut inner = || {
            if self.n_zeros == 0 {
                return None;
            }
            let hi = reader.read_f64::<LittleEndian>().ok()?;
            let lo = reader.read_f64::<LittleEndian>().ok()?;
            let z = T::zero() + hi + lo;

            self.n_zeros -= 1;
            Some((z, z.fp() * self.rtol + self.atol))
        };
        inner()
    }
}

#[cfg(test)]
mod tests {
    use super::{read_APCDB, write_APCDB, APCDB};
    use crate::lmfdb::read_LMFDB;
    use crate::zeta_zeros::{try_isolate, HybridPrecHardyZ};

    use log::debug;
    use F64x2::{f64x2, test_utils::assert_close};

    #[test]
    fn check_write_apcdb() {
        crate::init();
        let atol = 1e-18;
        let rtol = 1e-18;
        let t = 1e4;

        let mut hardy_z = HybridPrecHardyZ::<f64x2>::new(t, 10, atol);
        let stats = try_isolate(&mut hardy_z, -1, t, atol, rtol, true);
        let roots = stats.roots.iter().copied().filter(|&x| x.fp() <= t).collect::<Vec<_>>();
        let _ = write_APCDB(&roots, "data/apcdb/height_100000.dat", 0, roots.len(), atol, rtol);
    }

    #[test]
    fn check_rosser_rule() {
        crate::init();

        let n = 13999505 - 100;
        let atol = 1e-18;
        let rtol = 1e-18;
        let t = 6820053.0;

        let mut hardy_z = HybridPrecHardyZ::<f64x2>::new(t, 10, atol);
        let stats = try_isolate(&mut hardy_z, n, t, atol, rtol, true);
        let _roots = stats.roots.iter().copied().filter(|&x| x.fp() <= t).collect::<Vec<_>>();
    }

    #[test]
    fn check_read_apcdb() {
        crate::init();
        let limit = 1e4;

        let mut apcdb_zeros = vec![];
        let _ = read_APCDB("data/apcdb/height_100000.dat", limit, |x: f64x2| apcdb_zeros.push(x));

        let mut lmfdb_zeros = vec![];
        let _ = read_LMFDB("data/lmfdb", limit, |x: f64x2| lmfdb_zeros.push(x));

        let n1 = apcdb_zeros.len();
        let n2 = lmfdb_zeros.len();
        debug!("{n1} zeros from APCDB and {n2} zeros from LMFDB");
        for i in 0..n1.min(n2) {
            let (r1, r2) = (apcdb_zeros[i], lmfdb_zeros[i]);
            assert_close(r1, r2, 1e-18, 1e-18);
        }

        let apcdb2 = APCDB::<f64x2>::directory("data/apcdb")
            .map_while(|(x, _)| if x.fp() <= limit { Some(x) } else { None })
            .collect::<Vec<_>>();
        assert!(apcdb_zeros == apcdb2);
    }
}
