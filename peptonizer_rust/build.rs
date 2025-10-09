use std::fs::File;
use std::io::Write;
use std::path::Path;

const N: usize = 1024;
const MIN_X: f64 = 1e-10;
const MAX_X: f64 = 1.0;

fn main() {
    let step = (MAX_X - MIN_X) / (N as f64 - 1.0);

    let out_dir = std::env::var("OUT_DIR").unwrap();
    let path = Path::new(&out_dir).join("log_table.rs");
    let mut file = File::create(&path).unwrap();

    write!(file, "pub const LOG_TABLE: [f64; {}] = [\n", N).unwrap();
    for i in 0..N {
        let x = MIN_X + i as f64 * step;
        writeln!(file, "    {:.16},", x.ln()).unwrap();
    }
    writeln!(file, "];").unwrap();

    // Tell Cargo to rerun build.rs if this file changes (optional)
    println!("cargo:rerun-if-changed=build.rs");
}