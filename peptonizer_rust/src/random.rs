#[cfg(target_arch = "wasm32")]
use js_sys::Math;

use std::collections::HashSet;

/// Selects `n` unique random sample indices from a list of weights.
///
/// The selection is based on weighted probabilities, where higher weights increase
/// the likelihood of being chosen. Sampling is without replacement (indices are unique).
///
/// # Arguments
/// * `weights` - Vector of non-negative weights for each element.
/// * `n` - Number of unique indices to select.
///
/// # Returns
/// A `HashSet<usize>` containing the selected indices.
///
/// # Errors
/// This function will panic if:
/// * `weights` is empty.
/// * `n` is larger than `weights.len()`.
#[cfg(target_arch = "wasm32")]
pub fn select_random_samples_with_weights(
    weights: Vec<f64>,
    n: usize,
) -> Result<HashSet<usize>, Box<dyn std::error::Error>> {

    let cumulative_weights: Vec<f64> = weights
        .iter()
        .scan(0.0, |acc, &weight| {
            *acc += weight;
            Some(*acc)
        })
        .collect();

    let mut samples: HashSet<usize> = HashSet::with_capacity(n);

    while samples.len() < n {
        let r = Math::random() as f64; 
        if r >= cumulative_weights[cumulative_weights.len()-1] {
            continue;
        }

        let chosen_idx: usize = cumulative_weights.binary_search_by(|&w| w.partial_cmp(&r).expect("Partial compare returned None")).unwrap_or_else(|x| x);

        samples.insert(chosen_idx);
    }

    Ok(samples)
}

/// Selects `n` unique random sample indices from a list of weights.
///
/// The selection is based on weighted probabilities, where higher weights increase
/// the likelihood of being chosen. Sampling is without replacement (indices are unique).
///
/// # Arguments
/// * `weights` - Vector of non-negative weights for each element.
/// * `n` - Number of unique indices to select.
///
/// # Returns
/// A `HashSet<usize>` containing the selected indices.
///
/// # Errors
/// This function will panic if:
/// * `weights` is empty.
/// * `n` is larger than `weights.len()`.
#[cfg(not(target_arch = "wasm32"))]
pub fn select_random_samples_with_weights(
    weights: Vec<f64>,
    n: usize,
) -> Result<HashSet<usize>, Box<dyn std::error::Error>> {
    use rand::prelude::*;
    use rand::distributions::WeightedIndex;

    // Create a WeightedIndex for sampling
    let dist = WeightedIndex::new(&weights)?;

    let mut rng = thread_rng();
    let mut selected_indices = HashSet::new();

    while selected_indices.len() < n {
        let idx = dist.sample(&mut rng);
        selected_indices.insert(idx);
    }

    Ok(selected_indices)
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn test_basic_sampling() {
        let weights = vec![1.0, 2.0, 3.0];
        let samples = select_random_samples_with_weights(weights.clone(), 2);

        assert_eq!(samples.len(), 2);
        for &idx in &samples {
            assert!(idx < weights.len());
        }
    }

    #[test]
    fn test_full_sampling() {
        let weights = vec![0.5, 1.5, 2.0];
        let samples = select_random_samples_with_weights(weights.clone(), 3);

        // Must return all indices
        assert_eq!(samples, HashSet::from([0, 1, 2]));
    }

    #[test]
    fn test_heavy_weight_bias() {
        let weights = vec![1000.0, 0.0001, 0.0001];
        let mut counts = vec![0; 3];

        for _ in 0..100 {
            let s = select_random_samples_with_weights(weights.clone(), 1);
            let idx = *s.iter().next().unwrap();
            counts[idx] += 1;
        }

        // Index 0 should dominate
        assert!(counts[0] > 90);
    }
}
