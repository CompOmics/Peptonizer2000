extern crate serde_json;
extern crate serde;

mod utils;
mod http_client;
mod random;
mod weight_taxa;
pub mod zero_lookahead_belief_propagation;
mod node;
mod factor_graph;
mod messages;
mod convolution_tree;
mod array_utils;
mod fetch_unipept_taxa;
mod unipept_communicator;
mod taxa_clustering;
mod analyse_grid_search;
mod input_parser;
mod clean_csv;

#[cfg(target_arch = "wasm32")]
pub use wasm::*;

#[cfg(not(target_arch = "wasm32"))]
pub use pyo3::*;

#[cfg(target_arch = "wasm32")]
mod wasm {
    use wasm_bindgen::prelude::*;
    use crate::fetch_unipept_taxa::fetch_peptides_and_filter_taxa;
    use crate::weight_taxa::perform_taxa_weighing;
    use crate::zero_lookahead_belief_propagation::{run_belief_propagation, parse_taxon_scores};
    use crate::factor_graph::generate_graph;
    use crate::taxa_clustering::cluster_taxa;
    use crate::analyse_grid_search::compute_goodness;

    extern crate wasm_bindgen;
    extern crate web_sys;
    extern crate wasm_bindgen_futures;
    extern crate js_sys;
    extern crate console_error_panic_hook;

    /// Fetches taxa for peptides and filters them by rank and taxon query.
    ///
    /// # Arguments
    /// * `peptides` - JSON string of peptide sequences.
    /// * `rank` - Taxonomic rank used for filtering (e.g. "species").
    /// * `taxon_query` - JSON string of taxon IDs to filter against.
    ///
    /// # Returns
    /// JSON string mapping peptides to filtered taxon IDs.
    ///
    /// # Panics
    /// Panics if input JSON cannot be parsed or if result cannot be serialized.
    #[wasm_bindgen]
    pub fn fetch_unipept_taxa_wasm(
        peptides: String,
        rank: String,
        taxon_query: String
    ) -> String {
        fetch_peptides_and_filter_taxa(peptides, rank, taxon_query).unwrap()
    }

    /// Represents the main pipeline for weighting taxa based on peptide evidence.
    ///
    /// # Arguments
    ///
    /// * `pep_taxa` - JSON string mapping peptide sequences to lists of taxon IDs.
    /// * `pep_scores` - JSON string mapping peptide sequences to their scores (float).
    /// * `pep_psm_counts` - JSON string mapping peptide sequences to their PSM counts (int).
    /// * `max_taxa` - Maximum number of taxa to include in output.
    /// * `taxa_rank` - The taxonomic rank to normalize taxa to (e.g., "species").
    ///
    /// # Returns
    ///
    /// Tuple `(sequence_csv, taxa_weights_csv)`:
    /// * `sequence_csv` - CSV string of peptide sequences and their weights.
    /// * `taxa_weights_csv` - CSV string of taxa weights and uniqueness.
    #[wasm_bindgen]
    pub fn perform_taxa_weighing_wasm(
        pep_taxa: String,
        pep_scores: String,
        pep_psm_counts: String,
        max_taxa: usize,
        taxa_rank: String
    ) -> Box<[JsValue]> {
        console_error_panic_hook::set_once(); // Enable panic logging
        let (sequence_csv, taxa_weights_csv): (String, String) = perform_taxa_weighing(pep_taxa, pep_scores, pep_psm_counts, max_taxa, taxa_rank).unwrap();
        Box::new([JsValue::from(sequence_csv), JsValue::from(taxa_weights_csv)])
    }

    /// Generates a GraphML representation of a factor graph from a CSV string of taxon weights.
    ///
    /// # Arguments
    /// * `taxa_weights_csv` - A string containing CSV data for taxon weights.
    ///
    /// # Returns
    /// Returns a `Result` containing a GraphML string representation of the factor graph.
    ///
    /// # Errors
    /// Returns an error if CSV parsing fails or if any error occurs during graph construction.
    #[wasm_bindgen]
    pub fn generate_pepgm_graph_wasm(taxa_weights_csv: String) -> String {
        generate_graph(taxa_weights_csv).unwrap()
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
    #[wasm_bindgen]
    pub fn execute_pepgm_wasm(
        graph: String,
        alpha: f64,
        beta: f64,
        regularized: bool,
        prior: f64,
        max_iter: Option<i32>,
        tol: Option<f64>
    ) -> String {
        console_error_panic_hook::set_once(); // Enable panic logging
        let max_iter: i32 = max_iter.unwrap_or(10000);
        let tol: f64 = tol.unwrap_or(0.006);
        
        let csv: String = run_belief_propagation(graph, alpha, beta, regularized, prior, max_iter, tol).unwrap();

        parse_taxon_scores(csv).unwrap()
    }

    /// Clusters taxa based on peptidome similarity and returns a CSV.
    ///
    /// # Arguments
    /// * `graph_xml` - GraphML as string.
    /// * `taxa_weights_csv` - Taxa weights as CSV string.
    /// * `similarity_threshold` - Threshold for clustering.
    ///
    /// # Returns
    /// CSV string with taxa and their clusters.
    ///
    /// # Errors
    /// Returns an error if parsing, graph building, or clustering fails.
    #[wasm_bindgen]
    pub fn cluster_taxa_wasm(
        graph: String,
        taxa_weights_csv: String,
        similarity_threshold: f32
    ) -> String {
        cluster_taxa(graph, taxa_weights_csv, similarity_threshold).unwrap()
    }

    /// Computes a "goodness" score for clustering results by combining
    /// ranking similarity (via rank-biased overlap) and diversity (via entropy).
    /// 
    /// # Arguments
    /// * `clustered_taxa_weights_csv` - CSV string file containing clustered taxa weights.
    /// * `peptonizer_results` - JSON string containing taxa scores produced by Peptonizer.
    /// 
    /// # Returns
    /// A `Result<f64, Box<dyn std::error::Error>>` containing the computed goodness score,
    /// or an error if parsing fails.
    /// 
    /// # Errors
    /// This function may return an error if the input CSV or JSON cannot be parsed.
    #[wasm_bindgen]
    pub fn compute_goodness_wasm(
        clustered_taxa_weights_csv: String, 
        peptonizer_results: String
    ) -> f64 {
        compute_goodness(clustered_taxa_weights_csv, peptonizer_results).unwrap()
    }

}

#[cfg(not(target_arch = "wasm32"))]
mod pyo3 {
    use pyo3::prelude::*;
    use crate::fetch_unipept_taxa::fetch_peptides_and_filter_taxa;
    use crate::weight_taxa::perform_taxa_weighing;
    use crate::zero_lookahead_belief_propagation::{run_belief_propagation, parse_taxon_scores};
    use crate::factor_graph::generate_graph;
    use crate::taxa_clustering::cluster_taxa;
    use crate::analyse_grid_search::compute_goodness;
    use crate::input_parser::{parse_input_peptides, parse_unique_peptides};
    use crate::clean_csv::clean_csv;
    use crate::unipept_communicator::get_names_for_taxa;
    use serde_json;

    extern crate console_error_panic_hook;

    /// Parses peptides from a TSV string and returns JSON representations
    /// of scores and counts.
    ///
    /// # Arguments
    /// * `tsv_content` - Input TSV string with peptide data.
    ///
    /// # Returns
    /// A tuple containing:
    /// * `String` - JSON of peptide → max score mapping.
    /// * `String` - JSON of peptide → occurrence count mapping.
    ///
    /// # Errors
    /// Returns an error if parsing fails or if JSON serialization fails.
    #[pyfunction]
    pub fn parse_input_peptides_py(tsv_content: String) -> (String, String) {
        parse_input_peptides(tsv_content).unwrap()
    }

    /// Extracts unique peptides from a TSV string and returns them as JSON.
    ///
    /// # Arguments
    /// * `tsv_content` - Input TSV string with peptide data.
    ///
    /// # Returns
    /// JSON string containing the list of unique peptides.
    ///
    /// # Errors
    /// Returns an error if parsing fails or if JSON serialization fails.
    #[pyfunction]
    pub fn parse_unique_peptides_py(tsv_content: String) -> String {
        parse_unique_peptides(tsv_content).unwrap()
    }

    /// Fetches taxa for peptides and filters them by rank and taxon query.
    ///
    /// # Arguments
    /// * `peptides` - JSON string of peptide sequences.
    /// * `rank` - Taxonomic rank used for filtering (e.g. "species").
    /// * `taxon_query` - JSON string of taxon IDs to filter against.
    ///
    /// # Returns
    /// JSON string mapping peptides to filtered taxon IDs.
    ///
    /// # Panics
    /// Panics if input JSON cannot be parsed or if result cannot be serialized.
    #[pyfunction]
    pub fn fetch_unipept_taxa_py(
        peptides: String,
        rank: String,
        taxon_query: String
    ) -> String {
        fetch_peptides_and_filter_taxa(peptides, rank, taxon_query).unwrap()
    }

    /// Represents the main pipeline for weighting taxa based on peptide evidence.
    ///
    /// # Arguments
    ///
    /// * `pep_taxa` - JSON string mapping peptide sequences to lists of taxon IDs.
    /// * `pep_scores` - JSON string mapping peptide sequences to their scores (float).
    /// * `pep_psm_counts` - JSON string mapping peptide sequences to their PSM counts (int).
    /// * `max_taxa` - Maximum number of taxa to include in output.
    /// * `taxa_rank` - The taxonomic rank to normalize taxa to (e.g., "species").
    ///
    /// # Returns
    ///
    /// Tuple `(sequence_csv, taxa_weights_csv)`:
    /// * `sequence_csv` - CSV string of peptide sequences and their weights.
    /// * `taxa_weights_csv` - CSV string of taxa weights and uniqueness.
    #[pyfunction]
    fn perform_taxa_weighing_py(
        unipept_responses: String,
        pep_scores: String,
        pep_psm_counts: String,
        max_taxa: usize,
        taxa_rank: String
    ) -> (String, String) {
        perform_taxa_weighing(unipept_responses, pep_scores, pep_psm_counts, max_taxa, taxa_rank).unwrap()
    }

    /// Generates a GraphML representation of a factor graph from a CSV string of taxon weights.
    ///
    /// # Arguments
    /// * `taxa_weights_csv` - A string containing CSV data for taxon weights.
    ///
    /// # Returns
    /// Returns a `Result` containing a GraphML string representation of the factor graph.
    ///
    /// # Errors
    /// Returns an error if CSV parsing fails or if any error occurs during graph construction.
    #[pyfunction]
    pub fn generate_pepgm_graph_py(taxa_weights_csv: String) -> String {
        generate_graph(taxa_weights_csv).unwrap()
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
    #[pyfunction]
    pub fn execute_pepgm_py(
        graph: String,
        alpha: f64,
        beta: f64,
        regularized: bool,
        prior: f64,
        max_iter: Option<i32>,
        tol: Option<f64>
    ) -> String {
        console_error_panic_hook::set_once(); // Enable panic logging
        let max_iter: i32 = max_iter.unwrap_or(10000);
        let tol: f64 = tol.unwrap_or(0.006);
        
        run_belief_propagation(graph, alpha, beta, regularized, prior, max_iter, tol).unwrap()
    }

    /// Clusters taxa based on peptidome similarity and returns a CSV.
    ///
    /// # Arguments
    /// * `graph_xml` - GraphML as string.
    /// * `taxa_weights_csv` - Taxa weights as CSV string.
    /// * `similarity_threshold` - Threshold for clustering.
    ///
    /// # Returns
    /// CSV string with taxa and their clusters.
    ///
    /// # Errors
    /// Returns an error if parsing, graph building, or clustering fails.
    #[pyfunction]
    pub fn cluster_taxa_py(
        graph: String,
        taxa_weights_csv: String,
        similarity_threshold: f32
    ) -> String {
        cluster_taxa(graph, taxa_weights_csv, similarity_threshold).unwrap()
    }

    /// Computes a "goodness" score for clustering results by combining
    /// ranking similarity (via rank-biased overlap) and diversity (via entropy).
    /// 
    /// # Arguments
    /// * `clustered_taxa_weights_csv` - CSV string file containing clustered taxa weights.
    /// * `peptonizer_results_csv` - CSV string containing taxa scores produced by Peptonizer.
    /// 
    /// # Returns
    /// A `Result<f64, Box<dyn std::error::Error>>` containing the computed goodness score,
    /// or an error if parsing fails.
    /// 
    /// # Errors
    /// This function may return an error if the input CSV or JSON cannot be parsed.
    #[pyfunction]
    pub fn compute_goodness_py(
        clustered_taxa_weights_csv: String, 
        peptonizer_results_csv: String
    ) -> f64 {
        let taxon_scores = parse_taxon_scores(peptonizer_results_csv.clone()).unwrap();
        compute_goodness(clustered_taxa_weights_csv, taxon_scores).unwrap()
    }

    /// Returns a mapping from taxon ID to taxon name for all taxa provided.
    ///
    /// # Arguments
    /// * `target_taxa` - A list of taxon IDs for which all corresponding taxon names should be retrieved.
    ///
    /// # Errors
    /// Returns an error if the Unipept API server responds with a non-success status code
    /// or if something goes wrong with the network or JSON parsing.
    ///
    /// # Returns
    /// A JSON string mapping taxon IDs to their corresponding taxon names.
    #[pyfunction]
    pub fn get_names_for_taxa_py(target_taxa: Vec<i32>) -> String {
        let names = get_names_for_taxa(&target_taxa).unwrap();
        serde_json::to_string(&names).unwrap()
    }

    /// Read a CSV-file that was produced by the PepGM algorithm and use it to
    /// produce a new CSV-file that only contains the taxon-related information
    /// and scores. The string produced by this function can be written directly
    /// to a valid CSV-file and contains three columns: `taxon_name`, `taxon_id`,
    /// and `score`.
    ///
    /// # Arguments
    /// * `csv_content` - A CSV-file (as a string) that has been generated by running the PepGM algorithm.
    ///
    /// # Returns
    /// A `String` containing CSV rows with the columns: `taxon_name,taxon_id,score`.
    #[pyfunction]
    pub fn clean_csv_py(csv_content: String) -> String {
        clean_csv(csv_content).unwrap()
    }

    #[cfg(not(target_arch = "wasm32"))]
    #[pymodule]
    fn peptonizer_rust(_py: Python, m: &PyModule) -> PyResult<()> {
        m.add_function(wrap_pyfunction!(parse_input_peptides_py, m)?)?;
        m.add_function(wrap_pyfunction!(parse_unique_peptides_py, m)?)?;
        m.add_function(wrap_pyfunction!(fetch_unipept_taxa_py, m)?)?;
        m.add_function(wrap_pyfunction!(perform_taxa_weighing_py, m)?)?;
        m.add_function(wrap_pyfunction!(generate_pepgm_graph_py, m)?)?;
        m.add_function(wrap_pyfunction!(execute_pepgm_py, m)?)?;
        m.add_function(wrap_pyfunction!(cluster_taxa_py, m)?)?;
        m.add_function(wrap_pyfunction!(compute_goodness_py, m)?)?;
        m.add_function(wrap_pyfunction!(clean_csv_py, m)?)?;
        m.add_function(wrap_pyfunction!(get_names_for_taxa_py, m)?)?;
        Ok(())
    }
}
