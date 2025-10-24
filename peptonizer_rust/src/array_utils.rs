
include!(concat!(env!("OUT_DIR"), "/log_table.rs"));

const N: usize = 1024;
const MIN_X: f64 = 1e-10;
const MAX_X: f64 = 1.0;
const STEP: f64 = (MAX_X - MIN_X) / (N as f64 - 1.0);


/// Approximates the natural logarithm of `x` using a precomputed lookup table.
///
/// # Arguments
/// * `x` - Input value between `MIN_X` and `MAX_X`.
///
/// # Returns
/// * Approximate value of `ln(x)` retrieved from a static lookup table.
#[inline(always)]
pub fn ln_from_table(x: f64) -> f64 {
    let idx = ((x - MIN_X) / STEP) as usize;
    LOG_TABLE[idx]
}


/// Computes the sum of logarithms for pairs of values in batches to improve performance.
///
/// Each element in `rows` is a `[f64; 2]` pair. The function multiplies elements in
/// small batches (to minimize floating-point underflow) and applies logarithms to
/// the products using `ln_from_table`.
///
/// # Arguments
/// * `rows` - A vector of `[f64; 2]` pairs representing numeric data.
///
/// # Returns
/// * `[f64; 2]` where each component is the batched sum of logarithms across the corresponding column.
pub fn sum_logs_batched(rows: &Vec<[f64; 2]>) -> [f64; 2] {

    let mut acc0 = 0.0f64;
    let mut acc1 = 0.0f64;

    let block_size = 10;
    let blocks = rows.len() / block_size;

    for i in 0..blocks {
        let start = i * block_size;
        let mut prod0 = 1.0f64;
        let mut prod1 = 1.0f64;
        for &row in &rows[start..(start+block_size)] {
            prod0 *= row[0];
            prod1 *= row[1];
        }

        acc0 += ln_from_table(prod0);
        acc1 += ln_from_table(prod1);
    }

    let rem_start = blocks * block_size;
    let mut prod0 = 1.0;
    let mut prod1 = 1.0;
    for &row in &rows[rem_start..] {
        prod0 *= row[0];
        prod1 *= row[1];
    }

    [acc0 + prod0.ln(), acc1 + prod1.ln()]
}


/// Returns a copy of a vector excluding the element at the specified index.
///
/// # Arguments
/// * `v` - Reference to the input vector.
/// * `idx` - Index of the element to remove.
///
/// # Returns
/// * A new vector containing all elements of `v` except the one at position `idx`.
pub fn copy_without_index<T: Clone>(v: &Vec<T>, idx: usize) -> Vec<T> {
    let mut out = Vec::with_capacity(v.len() - 1);
    out.extend_from_slice(&v[..idx]);
    out.extend_from_slice(&v[idx + 1..]);
    out
}


/// Normalizes a vector of floating-point values so that the sum of all elements equals 1.
/// 
/// Mathematically: `x_i = x_i / Σx_j` for all elements `x_i` in the array.
/// 
/// # Arguments
/// * `array` - A mutable reference to a vector of `f64` values.
pub fn normalize(array: &mut Vec<f64>) {
    let sum: f64 = array.iter().sum();
    for val in array.iter_mut() {
        *val /= sum;
    }
}

pub fn normalize_arr(array: &mut [f64; 2]) {
    let sum: f64 = array[0] + array[1];
    array[0] /= sum;
    array[1] /= sum;
}


/// Normalizes a vector of 2D points so that the sum of all components equals 1.
/// 
/// Mathematically: For each `[a, b]`, `a = a / Σ(a_i + b_i)` and `b = b / Σ(a_i + b_i)`.
/// 
/// # Arguments
/// * `array` - A mutable reference to a vector of `[f64; 2]` values.
pub fn normalize_2d(array: &mut Vec<[f64; 2]>) {
    let sum: f64 = array.iter().map(|x| x[0] + x[1]).sum();
    for val in array.iter_mut() {
        val[0] /= sum;
        val[1] /= sum;
    }
}


/// Applies log-normalization to a [f64; 2] using the log-sum-exp trick for stability.
/// 
/// Mathematically: `x_i = exp(x_i - log(Σ exp(x_j)))`
/// 
/// # Arguments
/// * `array` - A mutable reference to a vector of `f64` values.
/// 
/// # Notes
/// - Subtracts `max(x)` to prevent overflow.
pub fn log_normalize(array: &mut [f64; 2]) {
    let max_val = array[0].max(array[1]);
    let log_sum_exp = ((array[0] - max_val).exp() + (array[1] - max_val).exp()).ln();
    array[0] = (array[0] - max_val - log_sum_exp).exp();
    array[1] = (array[1] - max_val - log_sum_exp).exp();
}


/// Applies log-normalization to a vector of 2D points using the log-sum-exp trick.
/// 
/// Mathematically: For each `[a, b]`:
/// `a = exp(a - log(Σ exp(exp(a_i) + exp(b_i))))`
/// `b = exp(b - log(Σ exp(exp(a_i) + exp(b_i))))`
/// 
/// # Arguments
/// * `array` - A mutable reference to a vector of `[f64; 2]` values.
/// 
/// # Notes
/// - Subtracts `max(a, b)` across the entire array for numerical stability.
pub fn log_normalize_2d(array: &mut Vec<[f64;2]>) {
    let max_val = array.iter().cloned().fold(f64::NEG_INFINITY, |acc, [a, b]| f64::max(f64::max(a, b), acc)); // Find the max value to prevent overflow
    let log_sum_exp = array.iter().map(|&[a, b]| (a - max_val).exp() + (b - max_val).exp()).sum::<f64>().ln(); // Calculate logsumexp

    array.iter_mut()
        .for_each(|x| {
            x[0] = (x[0] - max_val - log_sum_exp).exp();
            x[1] = (x[1] - max_val - log_sum_exp).exp();
        }); // Log-normalize and apply exp to each element
}


/// Prevents numerical underflow by setting a minimum threshold for all elements.
/// 
/// Mathematically: `x_i = max(x_i, 1e-30)`
/// 
/// # Arguments
/// * `array` - A mutable reference to a vector of `f64` values.
pub fn avoid_underflow(array: &mut Vec<f64>) {
    array.iter_mut().for_each(|x| if *x < 1e-30 { *x = 1e-30 });
}


/// Prevents numerical underflow in a fixed-size `[f64; 2]` array by setting a minimum threshold.
///
/// # Arguments
/// * `array` - Mutable reference to a `[f64; 2]` array.
pub fn avoid_underflow_arr(array: &mut [f64; 2]) {
    for i in 0..2 {
        if array[i] < 1e-30 {
            array[i] = 1e-30;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_normalize_basic() {
        let mut values = vec![1.0, 1.0, 2.0];
        normalize(&mut values);
        let sum: f64 = values.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_normalize_2d_basic() {
        let mut values = vec![[1.0, 1.0], [2.0, 2.0]];
        normalize_2d(&mut values);
        let total: f64 = values.iter().map(|x| x[0] + x[1]).sum();
        assert!((total - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_log_normalize_basic() {
        let mut values = [0.0, 0.0];
        log_normalize(&mut values);
        let sum: f64 = values[0] + values[1];
        assert!((sum - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_log_normalize_2d_basic() {
        let mut values = vec![[0.0, 0.0], [1.0, 1.0]];
        log_normalize_2d(&mut values);
        let total: f64 = values.iter().map(|x| x[0] + x[1]).sum();
        assert!((total - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_avoid_underflow_replaces_small_values() {
        let mut values = vec![1e-40, 1e-20];
        avoid_underflow(&mut values);
        assert!(values[0] >= 1e-30);
        assert!(values[1] >= 1e-30);
    }
}