
extern crate serde_json;
extern crate serde;

mod utils;
mod http_client;
mod random;
mod weight_taxa;
mod zero_lookahead_belief_propagation;
mod node;
mod factor_graph;
mod messages;
mod convolution_tree;
mod array_utils;
mod fetch_unipept_taxa;
mod unipept_communicator;

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

    extern crate wasm_bindgen;
    extern crate web_sys;
    extern crate wasm_bindgen_futures;
    extern crate js_sys;
    extern crate console_error_panic_hook;

    #[wasm_bindgen]
    pub fn fetch_unipept_taxa_wasm(
        peptides: String,
        rank: String,
        taxon_query: String
    ) -> String {
        fetch_peptides_and_filter_taxa(peptides, rank, taxon_query)
    }

    #[wasm_bindgen]
    pub fn perform_taxa_weighing_wasm(
        pep_taxa: String,
        pep_scores: String,
        pep_psm_counts: String,
        max_taxa: usize,
        taxa_rank: String
    ) -> Box<[JsValue]> {
        console_error_panic_hook::set_once(); // Enable panic logging
        let (sequence_csv, taxa_weights_csv): (String, String) = perform_taxa_weighing(pep_taxa, pep_scores, pep_psm_counts, max_taxa, taxa_rank);
        Box::new([JsValue::from(sequence_csv), JsValue::from(taxa_weights_csv)])
    }

    #[wasm_bindgen]
    pub fn generate_pepgm_graph_wasm(taxa_weights_csv: String) -> String {
        generate_graph(taxa_weights_csv).unwrap()
    }

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
        
        let csv: String = run_belief_propagation(graph, alpha, beta, regularized, prior, max_iter, tol);

        parse_taxon_scores(csv)
    }

}

#[cfg(not(target_arch = "wasm32"))]
mod pyo3 {
    use pyo3::prelude::*;
    use crate::weight_taxa::perform_taxa_weighing;

    #[pyfunction]
    fn perform_taxa_weighing_py(
        unipept_responses: String,
        pep_scores: String,
        pep_psm_counts: String,
        max_taxa: usize,
        taxa_rank: String
    ) -> (String, String) {
        perform_taxa_weighing(unipept_responses, pep_scores, pep_psm_counts, max_taxa, taxa_rank)
    }

    #[cfg(not(target_arch = "wasm32"))]
    #[pymodule]
    fn peptonizer_rust(_py: Python, m: &PyModule) -> PyResult<()> {
        m.add_function(wrap_pyfunction!(perform_taxa_weighing_py, m)?)?;
        Ok(())
    }
}
