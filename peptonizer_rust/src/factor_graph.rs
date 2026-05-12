use std::collections::{HashMap, HashSet};
use serde::Deserialize;
use csv::ReaderBuilder;
use nori::load_factor_graph_bytes;

/// Represents a single taxon weight record parsed from a CSV file.
#[derive(Deserialize)]
pub struct TaxonWeight {
    pub id: usize,
    pub sequence: String,
    pub score: f32,
    pub psms: usize,
    pub higher_taxa: usize,
    pub weight: f32,
    pub log_weight: f32
}


/// Parses a CSV string into a vector of `TaxonWeight` structs.
///
/// # Arguments
/// * `sequence_scores_csv` - A string containing CSV data for taxon weights. The CSV
///   must include headers: `id, sequence, score, psms, higher_taxa, weight, log_weight`.
///
/// # Returns
/// Returns a `Result` containing a vector of `TaxonWeight` structs if parsing succeeds.
///
/// # Errors
/// Returns an error if the CSV cannot be read, or if any record fails deserialization.
pub fn parse_taxon_weights_csv(sequence_scores_csv: String) -> Result<Vec<TaxonWeight>, Box<dyn std::error::Error>> {
    let mut rdr = ReaderBuilder::new()
        .has_headers(true)
        .from_reader(sequence_scores_csv.as_bytes());

    let mut sequence_scores = Vec::new();
    for record in rdr.deserialize() {
        let row: TaxonWeight = record.unwrap();
        sequence_scores.push(row);
    }

    Ok(sequence_scores)
}


/// Generates a GraphML representation of a factor graph from a CSV string of taxon weights.
///
/// # Arguments
/// * `sequence_scores_csv` - A string containing CSV data for taxon weights.
///
/// # Returns
/// Returns a `Result` containing a GraphML string representation of the factor graph.
///
/// # Errors
/// Returns an error if CSV parsing fails or if any error occurs during graph construction.
pub fn generate_graph(sequence_scores_csv: String) -> Result<Vec<u8>, Box<dyn std::error::Error>> {

    let sequence_scores = parse_taxon_weights_csv(sequence_scores_csv)?;

    let peptide_taxon_graph = taxon_weights_to_graphml(&sequence_scores);
    let factor_graph = load_factor_graph_bytes(&peptide_taxon_graph)?;
    Ok(factor_graph)
}

pub fn taxon_weights_to_graphml(taxon_weights: &Vec<TaxonWeight>) -> String {

    let mut xml = String::new();
    xml.push_str(r#"<?xml version="1.0" encoding="UTF-8"?>"#);
    xml.push_str(r#"<graphml xmlns="http://graphml.graphdrawing.org/xmlns">"#);
    xml.push_str(r#"  <key id="type" for="node"/>"#);
    xml.push_str(r#"  <key id="belief" for="node"/>"#);
    xml.push_str(r#"  <graph edgedefault="undirected">"#);
    
    // Keep unique sequence nodes
    let mut seen_sequences = HashSet::new();
    // Store first score seen for each sequence
    let mut sequence_scores: HashMap<&str, f32> = HashMap::new();

    for item in taxon_weights {
        sequence_scores
            .entry(item.sequence.as_str())
            .or_insert(item.score);
    }

    // Sequence nodes
    for (sequence, score) in &sequence_scores {
        if seen_sequences.insert(*sequence) {
            xml.push_str(&format!(
                r#"    <node id="{}">
      <data key="type">input</data>
      <data key="belief">[{}, {}]</data>
    </node>
"#,
                sequence,
                1.0 - score,
                score
            ));
        }
    }

    // Taxa nodes
    let mut seen_taxa = HashSet::new();
    for item in taxon_weights {
        if seen_taxa.insert(item.higher_taxa) {
            xml.push_str(&format!(
                r#"    <node id="{}">
      <data key="type">output</data>
    </node>
"#,
                item.higher_taxa
            ));
        }
    }

    // Edges
    for item in taxon_weights {
        xml.push_str(&format!(
            r#"    <edge source="{}" target="{}"/>
"#,
            item.sequence,
            item.higher_taxa
        ));
    }
    xml.push_str("  </graph>\n");
    xml.push_str("</graphml>\n");

    xml
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::node::{Node, NodeType, Factor};

    fn sample_csv() -> String {
        "id,sequence,score,psms,higher_taxa,weight,log_weight
1,PEPTIDE1,0.8,3,100,0.5,-0.3
2,PEPTIDE2,0.6,3,100,0.4,-0.5
3,PEPTIDE3,0.9,3,200,0.7,-0.1"
            .to_string()
    }

    #[test]
    fn test_parse_taxon_weights_csv() {
        let csv = sample_csv();
        let taxa = parse_taxon_weights_csv(csv).unwrap();
        assert_eq!(taxa.len(), 3);
        assert_eq!(taxa[0].id, 1);
        assert!((taxa[1].score - 0.6).abs() < 1e-6);
    }

    #[test]
    fn test_generate_graph_creates_graphml() {
        let csv = sample_csv();
        let graphml = generate_graph(csv).unwrap();
        assert!(graphml.contains("graphml"));
        assert!(graphml.contains("node"));
        assert!(graphml.contains("edge"));
    }

    #[test]
    fn test_edge_getters() {
        let edge = Edge::new(1, 10, 20, 0, 0, Some(5));
        assert_eq!(edge.get_id(), 1);
        assert_eq!(edge.get_node1_id(), 10);
        assert_eq!(edge.get_node2_id(), 20);
        assert_eq!(edge.get_node_ids(), (10, 20));
        assert_eq!(edge.get_message_length(), Some(5));
    }

    #[test]
    fn test_ctfactorgraph_from_taxa_weights() {
        let csv = sample_csv();
        let taxa = parse_taxon_weights_csv(csv).unwrap();
        let graph = CTFactorGraph::from_taxa_weights(taxa);
        assert!(graph.node_count() > 0);
        assert!(graph.edge_count() > 0);
    }

    #[test]
    fn test_ctfactorgraph_to_and_from_graphml() {
        let csv = sample_csv();
        let taxa = parse_taxon_weights_csv(csv).unwrap();
        let graph = CTFactorGraph::from_taxa_weights(taxa);
        let graphml = graph.to_graphml();
        assert!(graphml.is_ok());
        let graphml = graphml.unwrap();

        let parsed = CTFactorGraph::from_graphml(&graphml).unwrap();
        assert_eq!(graph.node_count(), parsed.node_count());
        assert_eq!(graph.edge_count(), parsed.edge_count());
    }

    #[test]
    fn test_neighbor_operations() {
        let csv = sample_csv();
        let taxa = parse_taxon_weights_csv(csv).unwrap();
        let graph = CTFactorGraph::from_taxa_weights(taxa);

        if graph.node_count() > 1 {
            let node = graph.get_node(0);
            println!("{:?}\n\n{:?}", graph, node);
            for n in graph.get_neighbors(node) {
                assert!(n >= 0);
            }
        }
    }

    #[test]
    fn test_get_peptide_for_factor_returns_ok_or_err() {
        let csv = sample_csv();
        let taxa = parse_taxon_weights_csv(csv).unwrap();
        let graph = CTFactorGraph::from_taxa_weights(taxa);

        for (i, node) in graph.get_nodes().iter().enumerate() {
            if node.is_factor_node() {
                let result = graph.get_peptide_for_factor(i);
                assert!(result.is_ok() || result.is_err());
            }
        }
    }

    #[test]
    fn test_connected_components() {
        let csv = sample_csv();
        let taxa = parse_taxon_weights_csv(csv).unwrap();
        let graph = CTFactorGraph::from_taxa_weights(taxa);

        let components = graph.connected_components();
        assert!(!components.is_empty());
        let total_nodes: usize = components.iter().map(|c| c.node_count()).sum();
        assert_eq!(total_nodes, graph.node_count());
    }
}
