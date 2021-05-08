use std::io::BufRead;

use log::{debug, info};

use crate::traits::MyReal;
use byteorder::{LittleEndian, ReadBytesExt};

const LMFDB_DATA_PATH: &str = "./data/zeros";

fn LMFDB_read_ckpt_list() -> Result<Vec<i64>, std::io::Error> {
    let file = std::fs::File::open(format!("{}/md5.txt", LMFDB_DATA_PATH))?;
    let reader = std::io::BufReader::new(file);
    let mut indices = vec![];
    for line in reader.lines() {
        let line = line?;
        let length = line.len();
        let index = &line[40..line.len() - 4];
        let index = index.parse::<i64>().unwrap();
        indices.push(index);
    }
    indices.sort_unstable();
    Ok(indices)
}

// Actually, I think it's better not to take a function as input, but return an Iterator.
// However, I'm too lazy to convert it to an iterator. Wish Generator was stablized.
pub(crate) fn LMFDB_reader<T: MyReal, F: FnMut(T)>(
    limit: f64,
    mut f: F,
) -> Result<(), std::io::Error> {
    info!("Loading zeta zeros up to {}", limit);
    let ckpts = LMFDB_read_ckpt_list()?;

    let eps = T::from_f64(2.0).unwrap().powi(-101);
    let limit = T::from_f64(limit).unwrap();

    let mut n_roots = 0usize;
    for &ckpt in ckpts.iter() {
        let path = format!("{}/zeros_{}.dat", LMFDB_DATA_PATH, ckpt);
        let file = std::fs::File::open(&path)?;
        let mut reader = std::io::BufReader::new(file);
        debug!("[LMFDB] loading {}", path);

        let n_blocks = reader.read_u64::<LittleEndian>()?;

        for b in 0..n_blocks {
            let t0 = reader.read_f64::<LittleEndian>()?;
            let t1 = reader.read_f64::<LittleEndian>()?;
            let n0 = reader.read_u64::<LittleEndian>()?;
            let n1 = reader.read_u64::<LittleEndian>()?;

            let t0 = T::from_f64(t0).unwrap();
            let mut z = 0u128;

            for i in n0..n1 {
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
