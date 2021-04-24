use log::{debug, info};

use crate::traits::MyReal;
use byteorder::{LittleEndian, ReadBytesExt};

pub(crate) fn LMFDB_reader<T: MyReal>(limit: T) -> Result<Vec<T>, std::io::Error> {
    let data_files = [
        "./data/zeros/zeros_14.dat",
        "./data/zeros/zeros_5000.dat",
        "./data/zeros/zeros_26000.dat",
        "./data/zeros/zeros_236000.dat",
        "./data/zeros/zeros_446000.dat",
        "./data/zeros/zeros_2546000.dat",
        "./data/zeros/zeros_4646000.dat",
    ];
    info!("Loading zeta zeros up to {}", limit);

    let eps = 2.0.unchecked_cast::<T>().powi(-101);

    let mut ret = vec![];
    for &file_name in data_files.iter() {
        let mut file = std::fs::File::open(file_name)?;
        let n_blocks = file.read_u64::<LittleEndian>()?;

        for b in 0..n_blocks {
            let t0 = file.read_f64::<LittleEndian>()?;
            let t1 = file.read_f64::<LittleEndian>()?;
            let n0 = file.read_u64::<LittleEndian>()?;
            let n1 = file.read_u64::<LittleEndian>()?;
            debug!(
                "[LMFDB] loading {} block {}, from N({}) = {} to N({}) = {}",
                file_name, b, t0, n0, t1, n1
            );

            let t0 = T::from_f64(t0).unwrap();
            let mut z = 0u128;

            for i in n0..n1 {
                let z1 = file.read_u64::<LittleEndian>()? as u128;
                let z2 = file.read_u32::<LittleEndian>()? as u128;
                let z3 = file.read_u8()? as u128;
                z = z + z1 + (z2 << 64) + (z3 << 96);

                let zz = t0 + T::from_u128(z).unwrap() * eps;

                if zz > limit {
                    return Ok(ret);
                }
                // debug!("read zero: {}", z);
                ret.push(zz);
            }
        }
    }
    panic!("Insufficient zeta zeros data")
}
