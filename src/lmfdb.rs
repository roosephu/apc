use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::{io::BufRead, path::PathBuf};

use log::{debug, info};

use crate::traits::MyReal;
use byteorder::{LittleEndian, ReadBytesExt};

/// Wish if this can be turned into a Generator and then an Iterator.
pub(crate) fn read_LMFDB<T: MyReal>(
    data_path: impl AsRef<Path>,
    limit: f64,
    mut f: impl FnMut(T),
) -> Result<(), std::io::Error> {
    info!("Loading zeta zeros up to {}", limit);

    let data_path = data_path.as_ref();
    let file = std::fs::File::open(data_path.join("md5.txt")).expect("Missing md5.txt");
    let reader = std::io::BufReader::new(file);
    let mut indices = vec![];
    for line in reader.lines() {
        let line = line?;
        let index = &line[40..line.len() - 4];
        let index = index.parse::<i64>().unwrap();
        indices.push(index);
    }
    indices.sort_unstable();

    let eps = T::mp(2.0).powi(-101);
    let limit = T::mp(limit);

    let mut n_roots = 0usize;
    for &index in indices.iter() {
        let path = data_path.join(format!("zeros_{}.dat", index));
        let file = std::fs::File::open(&path)?;
        let mut reader = std::io::BufReader::new(file);
        debug!("[LMFDB] loading {}", path.display());

        let n_blocks = reader.read_u64::<LittleEndian>()?;

        for _ in 0..n_blocks {
            let t0 = reader.read_f64::<LittleEndian>()?;
            let _ = reader.read_f64::<LittleEndian>()?;
            let n0 = reader.read_u64::<LittleEndian>()?;
            let n1 = reader.read_u64::<LittleEndian>()?;
            debug!("n0 = {n0}, n1 = {n1}");

            let t0 = T::mp(t0);
            let mut z = 0u128;

            for _ in n0..n1 {
                let z1 = reader.read_u64::<LittleEndian>()? as u128;
                let z2 = reader.read_u32::<LittleEndian>()? as u128;
                let z3 = reader.read_u8()? as u128;
                z = z + z1 + (z2 << 64) + (z3 << 96);

                let zz = t0 + T::from_u128(z).unwrap() * eps;

                if zz > limit {
                    info!("# zeros = {}", n_roots);
                    return Ok(());
                }
                n_roots += 1;
                f(zz);
            }
        }
    }
    panic!("Insufficient zeta zeros data")
}

#[allow(clippy::upper_case_acronyms)]
pub struct LMFDB<T> {
    reader: BufReader<File>,
    n_blocks: u64,
    n_zeros: u64,
    z: u128,
    t0: T,
    eps: T,
}

impl<T: MyReal> LMFDB<T> {
    pub fn new(path: impl AsRef<Path>) -> std::io::Result<Self> {
        let path = path.as_ref();
        let file = std::fs::File::open(&path)?;
        let mut reader = BufReader::new(file);
        debug!("[LMFDB] loading {}", path.display());

        let n_blocks = reader.read_u64::<LittleEndian>()?;
        let eps = T::mp(2.0).powi(-101);

        Ok(Self { eps, n_blocks, n_zeros: 0, reader, z: 0, t0: T::mp(0.0) })
    }

    pub fn directory(path: impl AsRef<Path>) -> impl Iterator<Item = (T, f64)> {
        let path = path.as_ref().to_path_buf();
        let file = std::fs::File::open(path.join("md5.txt")).unwrap();
        let reader = std::io::BufReader::new(file);
        let mut indices = vec![];
        for line in reader.lines() {
            let line = line.unwrap();
            let index = &line[40..line.len() - 4];
            let index = index.parse::<i64>().unwrap();
            indices.push(index);
        }
        indices.sort_unstable();

        indices
            .into_iter()
            .flat_map(move |index| Self::new(path.join(format!("zeros_{}.dat", index))).unwrap())
    }
}

impl<T: MyReal> Iterator for LMFDB<T> {
    type Item = (T, f64);

    fn next(&mut self) -> Option<Self::Item> {
        let reader = &mut self.reader;
        let mut inner = || loop {
            if self.n_zeros != 0 {
                let z1 = reader.read_u64::<LittleEndian>().ok()? as u128;
                let z2 = reader.read_u32::<LittleEndian>().ok()? as u128;
                let z3 = reader.read_u8().ok()? as u128;
                self.z = self.z + z1 + (z2 << 64) + (z3 << 96);

                let zz = self.t0 + T::from_u128(self.z)? * self.eps;
                self.n_zeros -= 1;

                break Some((zz, self.eps.fp()));
            } else if self.n_blocks != 0 {
                let t0 = reader.read_f64::<LittleEndian>().ok()?;
                let _ = reader.read_f64::<LittleEndian>().ok()?;
                let n0 = reader.read_u64::<LittleEndian>().ok()?;
                let n1 = reader.read_u64::<LittleEndian>().ok()?;

                self.t0 = T::mp(t0);
                self.z = 0u128;
                self.n_zeros = n1 - n0;
                self.n_blocks -= 1;
            } else {
                break None;
            }
        };
        inner()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lmfdb_iter() {
        crate::init();

        let limit = 1e7;
        let path = "./data/lmfdb";

        let mut lmfdb = vec![];
        let _ = read_LMFDB::<f64>(path, limit, |x| lmfdb.push(x));

        let db = LMFDB::<f64>::directory(path);
        let mut lmfdb2 = vec![];
        for (x, _) in db {
            if x > limit {
                break;
            }
            lmfdb2.push(x)
        }
        println!("len = {}, {}", lmfdb.len(), lmfdb2.len());

        assert!(lmfdb == lmfdb2);
    }
}
