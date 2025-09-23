
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


/// Applies log-normalization to a vector using the log-sum-exp trick for stability.
/// 
/// Mathematically: `x_i = exp(x_i - log(Σ exp(x_j)))`
/// 
/// # Arguments
/// * `array` - A mutable reference to a vector of `f64` values.
/// 
/// # Notes
/// - Subtracts `max(x)` to prevent overflow.
pub fn log_normalize(array: &mut Vec<f64>) {
    let max_val = array.iter().cloned().fold(f64::NEG_INFINITY, f64::max); // Find the max value to prevent overflow
    let log_sum_exp = array.iter().map(|&x| (x - max_val).exp()).sum::<f64>().ln(); // Calculate logsumexp

    array.iter_mut()
        .for_each(|x| *x = (*x - max_val - log_sum_exp).exp()); // Log-normalize and apply exp to each element
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