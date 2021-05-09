use std::io::BufRead;
use std::path::Path;

use log::{debug, info};

use crate::traits::MyReal;
use byteorder::{LittleEndian, ReadBytesExt};

// Actually, I think it's better not to take a function as input, but return an Iterator.
// However, I'm too lazy to convert it to an iterator. Wish Generator was stablized.
pub(crate) fn LMFDB_reader<T: MyReal, F: FnMut(T)>(
    data_path: &Path,
    limit: f64,
    mut f: F,
) -> Result<(), std::io::Error> {
    info!("Loading zeta zeros up to {}", limit);

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

    let eps = T::from_f64(2.0).unwrap().powi(-101);
    let limit = T::from_f64(limit).unwrap();

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

            let t0 = T::from_f64(t0).unwrap();
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
