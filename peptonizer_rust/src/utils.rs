/// Cross-platform logging utilities for Rust/WASM.
///
/// These functions provide a unified interface for logging messages to the console,
/// both when running natively (non-WASM) and in a WebAssembly (WASM) environment.

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
extern "C" {
    /// Logs a message to the JavaScript console (`console.log`) when running in WASM.
    #[wasm_bindgen(js_namespace = console)]
    pub fn log(s: &str);
}

/// Logs a message to the console when running natively (non-WASM).
#[cfg(not(target_arch = "wasm32"))]
pub fn log(s: &str) {
    println!("{s}");
}

