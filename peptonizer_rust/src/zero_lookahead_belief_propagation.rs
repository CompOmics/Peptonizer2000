use crate::factor_graph::CTFactorGraph;
use std::collections::HashMap;
use crate::messages::Messages;
use csv::Writer;
use serde_json;
use std::io::Cursor;
use csv::ReaderBuilder;


/// Calibrates multiple subgraphs (connected components) of a factor graph using loopy belief propagation.
///
/// # Arguments
///
/// * `ct_factor_graphs` - Vector of `CTFactorGraph` objects representing connected subgraphs of the full factor graph.
/// * `max_iterations` - Maximum number of iterations for message passing in case of non-convergence.
/// * `tolerance` - Convergence criterion; the maximum allowable change in messages between iterations.
///
/// # Returns
///
/// Tuple `(node_names, node_categories, results)`:
/// * `node_names` - Vector of node names in the same order as the belief results.
/// * `node_categories` - Vector of node types (categories) corresponding to the nodes.
/// * `results` - Vector of belief distributions for each node; each element is a vector `[P(0), P(1)]`.
fn calibrate_all_subgraphs(
    ct_factor_graphs: Vec<CTFactorGraph>,
    max_iterations: u32,
    tolerance: f64
) -> Result<(Vec<String>, Vec<String>, Vec<Vec<f64>>), Box<dyn std::error::Error>>{
    let mut results: Vec<Vec<f64>> = Vec::new();
    let mut node_categories: Vec<String> = Vec::new();
    let mut node_names: Vec<String> = Vec::new();

    for subgraph in ct_factor_graphs {
        if subgraph.node_count() > 2 {

            subgraph.add_node_names_categories(&mut node_names, &mut node_categories);

            let mut messages = Messages::new(subgraph);
            let beliefs: Vec<Vec<f64>> = messages.zero_lookahead_bp(
                max_iterations,
                tolerance
            )?;

            results.extend(beliefs);
        }
    }

    Ok((node_names, node_categories, results))
}


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
/// * `prior` - Prior belief for taxon nodes.
/// * `max_iter` - Maximum number of belief propagation iterations.
/// * `tol` - Tolerance threshold for message convergence.
///
/// # Returns
///
/// CSV string with one row per node containing columns:
/// `[node_name, posterior_probability_1, node_category]`
pub fn run_belief_propagation(
    graph: String,
    alpha: f64,
    beta: f64,
    regularized: bool,
    prior: f64,
    max_iter: u32,
    tol: f64
) -> Result<String, Box<dyn std::error::Error>> {
    let mut ct_factor_graph = CTFactorGraph::from_graphml(&graph)?;
    ct_factor_graph.fill_in_factors(alpha, beta, regularized);
    ct_factor_graph.fill_in_priors(prior);
    ct_factor_graph.add_ct_nodes();
    let ct_factor_graphs: Vec<CTFactorGraph> = ct_factor_graph.connected_components();

    let (node_names, node_types, results) = calibrate_all_subgraphs(
        ct_factor_graphs,
        max_iter,
        tol
    )?;

    Ok(generate_csv(node_names, node_types, results)?)
}


/// Generates a CSV string from node names, types, and belief results.
///
/// # Arguments
///
/// * `node_names` - Vector of node names.
/// * `node_types` - Vector of node types (categories) corresponding to `node_names`.
/// * `results` - Vector of belief distributions for each node; each element is a vector `[P(0), P(1)]`.
///
/// # Returns
///
/// CSV string with columns `[node_name, posterior_probability_1, node_category]`.
fn generate_csv(node_names: Vec<String>, node_types: Vec<String>, results: Vec<Vec<f64>>) -> Result<String, Box<dyn std::error::Error>> {

    let mut wtr = Writer::from_writer(vec![]);

    for i in 0..node_names.len() {
        let _ = wtr.write_record(&[
            node_names[i].clone(),
            results[i][1].to_string(),
            node_types[i].clone()
        ])?;
    }

    let csv: String = String::from_utf8(wtr.into_inner()?)?;

    Ok(csv)
}


/// Parses a CSV string of belief propagation results and extracts taxon scores.
///
/// Only rows with type "taxon" are included. The results are sorted by score in ascending order.
///
/// # Arguments
///
/// * `csv_content` - CSV string with columns `[id, score, type]`.
///
/// # Returns
///
/// JSON string mapping taxon IDs (`usize`) to their posterior probabilities (`f64`), sorted by score.
pub fn parse_taxon_scores(csv_content: String) -> Result<String, Box<dyn std::error::Error>> {
    let mut rdr = ReaderBuilder::new()
        .has_headers(false)
        .from_reader(Cursor::new(csv_content));

    let mut taxon_score_dict = HashMap::new();
    let mut records = Vec::new();

    for result in rdr.records() {
        let record = result?;
        
        let record_type = record.get(2).ok_or("Index 2 not in record")?;
        
        // Filter rows where "type" == "taxon"
        if record_type == "taxon" {
            let id: usize = record.get(0).ok_or("Index 0 not in record")?.parse()?;
            let score: f64 = record.get(1).ok_or("Index 1 not in record")?.parse()?;
            records.push((id, score));
        }
    }

    // Sort by score in ascending order
    records.sort_by(|a, b| a.1.partial_cmp(&b.1).expect("Partial compare returned None"));

    // Populate the HashMap with sorted values
    for (id, score) in records {
        taxon_score_dict.insert(id, score);
    }

    Ok(serde_json::to_string(&taxon_score_dict)?)
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn test_generate_csv_basic() {
        let node_names = vec!["n1".to_string(), "n2".to_string()];
        let node_types = vec!["taxon".to_string(), "peptide".to_string()];
        let results = vec![vec![0.3, 0.7], vec![0.6, 0.4]];

        let csv = generate_csv(node_names, node_types, results);

        assert!(csv.contains("n1"));
        assert!(csv.contains("0.7"));
        assert!(csv.contains("taxon"));
    }

    #[test]
    fn test_parse_taxon_scores_basic() {
        let csv_content = "123,0.8,taxon\n456,0.5,taxon\n789,0.9,peptide\n".to_string();
        let json = parse_taxon_scores(csv_content);

        let parsed: HashMap<usize, f64> = serde_json::from_str(&json).unwrap();
        assert_eq!(parsed.get(&456), Some(&0.5));
        assert_eq!(parsed.get(&123), Some(&0.8));
        assert!(parsed.get(&789).is_none());
    }

    #[test]
    fn test_calibrate_all_subgraphs_empty() {
        let (names, cats, results) = calibrate_all_subgraphs(vec![], 10, 1e-6);
        assert!(names.is_empty());
        assert!(cats.is_empty());
        assert!(results.is_empty());
    }

    #[test]
    fn test_run_belief_propagation_does_not_crash() {
        let minimal_graph = r#"<?xml version='1.0' encoding='utf-8'?>
        <graphml xmlns="http://graphml.graphdrawing.org/xmlns">
            <key id="d1" for="node" attr.name="InitialBelief_1" attr.type="double"/>
            <key id="d0" for="node" attr.name="InitialBelief_0" attr.type="double"/>
            <graph edgedefault="undirected">
                <node id="n0">
                    <data key="d0">0.0010000000000000009</data>
                    <data key="d1">0.999</data>
                    <data key="d2">peptide</data>
                </node>
                <node id="n1">
                    <data key="d0">0.0010000000000000009</data>
                    <data key="d1">0.999</data>
                    <data key="d2">peptide</data>
                </node>
                <node id="n2">
                    <data key="d2">taxon</data>
                </node>
                <edge source="n0" target="n1"/>
                <edge source="n1" target="n2"/>
            </graph>
        </graphml>
        "#.to_string();

        let csv = run_belief_propagation(
            minimal_graph,
            0.5,   // alpha
            0.5,   // beta
            true,  // regularized
            0.1,   // prior
            10,    // max_iter
            1e-6   // tol
        );

        assert!(csv.contains("n0"));
    }
}
