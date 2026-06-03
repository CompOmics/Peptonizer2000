use std::collections::HashMap;
use nori::zero_lookahead_bp_from_graph_bytes;


/// Runs belief propagation on a factor graph provided as a GraphML string.
///
/// This function constructs the factor graph, fills in factor tables and priors,
/// splits the graph into connected components, and performs loopy belief propagation
/// on each component. The result is returned as a CSV string.
///
/// # Arguments
///
/// * `graph` - GraphML representation of the factor graph.
/// * `alpha` - Noisy-OR factor alpha parameter.
/// * `beta` - Noisy-OR factor beta parameter.
/// * `regularized` - Whether to regularize factor tables to penalize large numbers of parents.
/// * `prior` - Prior belief for effect nodes.
/// * `max_iter` - Maximum number of belief propagation iterations.
/// * `tol` - Tolerance threshold for message convergence.
///
/// # Returns
///
/// CSV string with one row per node containing columns:
/// `[node_name, posterior_probability_1, node_category]`
pub fn run_belief_propagation(
    graph_bytes: &[u8],
    alpha: f32,
    beta: f32,
    regularized: bool,
    prior: f32,
    max_iter: Option<u32>,
    tol: Option<f32>
) -> Result<String, Box<dyn std::error::Error>> {
    let results = zero_lookahead_bp_from_graph_bytes(graph_bytes, alpha, beta, regularized, prior, max_iter, tol).unwrap();

    let effect_score_dict: HashMap<String, f32> = results
        .into_iter()
        .filter_map(|(key, values)| {
            values.get(1).map(|&v| (key, v))
        })
        .collect();

    Ok(serde_json::to_string(&effect_score_dict)?)
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn test_parse_effect_scores_basic() {
        let csv_content = "123,0.8,effect\n456,0.5,effect\n789,0.9,peptide\n".to_string();
        let json = parse_effect_scores(csv_content);
        assert!(json.is_ok());
        let json = json.unwrap();

        let parsed: HashMap<usize, f64> = serde_json::from_str(&json).unwrap();
        assert_eq!(parsed.get(&456), Some(&0.5));
        assert_eq!(parsed.get(&123), Some(&0.8));
        assert!(parsed.get(&789).is_none());
    }
}
