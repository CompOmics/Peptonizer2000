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
    max_iterations: i32,
    tolerance: f64
) -> (Vec<String>, Vec<String>, Vec<Vec<f64>>){
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
            );

            results.extend(beliefs);
        }
    }

    (node_names, node_categories, results)
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
    max_iter: i32,
    tol: f64
) -> String {
    let mut ct_factor_graph = CTFactorGraph::from_graphml(&graph).unwrap();
    ct_factor_graph.fill_in_factors(alpha, beta, regularized);
    ct_factor_graph.fill_in_priors(prior);
    ct_factor_graph.add_ct_nodes();
    let ct_factor_graphs: Vec<CTFactorGraph> = ct_factor_graph.connected_components();

    let (node_names, node_types, results) = calibrate_all_subgraphs(
        ct_factor_graphs,
        max_iter,
        tol
    );

    generate_csv(node_names, node_types, results)
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
fn generate_csv(node_names: Vec<String>, node_types: Vec<String>, results: Vec<Vec<f64>>) -> String {

    let mut wtr = Writer::from_writer(vec![]);

    for i in 0..node_names.len() {
        let _ = wtr.write_record(&[
            node_names[i].clone(),
            results[i][1].to_string(),
            node_types[i].clone()
        ]).unwrap();
    }

    let csv: String = String::from_utf8(wtr.into_inner().unwrap()).unwrap();

    csv
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
/// JSON string mapping taxon IDs (`i32`) to their posterior probabilities (`f64`), sorted by score.
pub fn parse_taxon_scores(csv_content: String) -> String {
    let mut rdr = ReaderBuilder::new()
        .has_headers(false)
        .from_reader(Cursor::new(csv_content));

    let mut taxon_score_dict = HashMap::new();
    let mut records = Vec::new();

    for result in rdr.records() {
        let record = result.unwrap();
        
        let record_type = record.get(2).unwrap();
        
        // Filter rows where "type" == "taxon"
        if record_type == "taxon" {
            let id: i32 = record.get(0).unwrap().parse().unwrap();
            let score: f64 = record.get(1).unwrap().parse().unwrap();
            records.push((id, score));
        }
    }

    // Sort by score in ascending order
    records.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    // Populate the HashMap with sorted values
    for (id, score) in records {
        taxon_score_dict.insert(id, score);
    }

    serde_json::to_string(&taxon_score_dict).unwrap()
}