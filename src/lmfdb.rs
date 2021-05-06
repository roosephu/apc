use std::io::BufRead;

use log::{debug, info};

use crate::traits::MyReal;
use byteorder::{LittleEndian, ReadBytesExt};

/// https://stackoverflow.com/questions/40107797/can-array-lengths-be-inferred-in-rust
macro_rules! arr {
    ($id: ident $name: ident: [$ty: ty; _] = $value: expr) => {
        $id $name: [$ty; $value.len()] = $value;
    }
}

const LMFDB_DATA_PATH: &str = "./data/zeros";
arr!(const LMFDB_CKPTS: [i64; _] = [14, 5000, 26000, 236000, 446000, 2546000, 4646000, 6746000, 8846000]);

pub fn LMFDB_read_ckpt_list() -> Result<Vec<i64>, std::io::Error> {
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

pub(crate) fn LMFDB_reader<T: MyReal>(limit: f64) -> Result<Vec<T>, std::io::Error> {
    info!("Loading zeta zeros up to {}", limit);
    let ckpts = LMFDB_read_ckpt_list()?;

    let eps = T::from_f64(2.0).unwrap().powi(-101);
    let limit = T::from_f64(limit).unwrap();

    let mut roots = vec![];
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
                    info!("largest zeta roots {:?}, # zeros = {}", roots.last(), roots.len());
                    return Ok(roots);
                }
                // debug!("read zero: {}", z);
                roots.push(zz);
            }
        }
    }
    panic!("Insufficient zeta zeros data")
}
