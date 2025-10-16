use crate::factor_graph::CTFactorGraph;
use std::collections::{HashMap, HashSet};
use serde::{Serialize, Deserialize, Deserializer};
use csv::{ReaderBuilder, WriterBuilder};



/// Represents a taxonomic unit with attributes used for clustering.
#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct Taxon {
    /// Unique identifier of the taxon.
    pub id: i32,
     /// Identifier of the higher-level taxon it belongs to.
    pub higher_taxa: i32,
    /// Weight of the taxon, scaled for comparison.
    pub scaled_weight: f32,
    /// Whether the taxon is unique in the dataset.
    pub unique: bool,

    /// IDs of taxa belonging to the same cluster.
    /// Serialized as a string, deserialized back into a vector.
    #[serde(default, serialize_with = "vec_to_string", deserialize_with = "string_to_vec")]
    pub cluster_members: Vec<i32>
}


/// Converts a `Vec<i32>` to a serialized string representation.
///
/// # Arguments
/// * `vec` - Reference to the vector to be serialized.
/// * `serializer` - Serializer provided by Serde.
///
/// # Returns
/// Serialized string wrapped in the serializer's `Ok` type.
///
/// # Errors
/// Returns an error if serialization fails.
fn vec_to_string<S>(vec: &Vec<i32>, serializer: S) -> Result<S::Ok, S::Error>
where
    S: serde::Serializer,
{
    let joined = &format!(
        "[{}]", 
        vec.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", ")
    );
    serializer.serialize_str(&joined)
}


/// Converts a serialized string back into a `Vec<i32>`.
///
/// # Arguments
/// * `deserializer` - Deserializer provided by Serde.
///
/// # Returns
/// Vector of integers parsed from the string.
///
/// # Errors
/// Returns an error if the string cannot be deserialized or parsed into integers.
pub fn string_to_vec<'de, D>(deserializer: D) -> Result<Vec<i32>, D::Error>
where
    D: Deserializer<'de>,
{
    let s: &str = Deserialize::deserialize(deserializer)?;
    let vec = s[1..(s.len()-1)].split(',')
               .filter_map(|item| item.trim().parse::<i32>().ok())
               .collect();
    Ok(vec)
}


/// Parses a CSV string into a list of `Taxon`.
///
/// # Arguments
/// * `taxa_weights_csv` - CSV content as string.
///
/// # Returns
/// Vector of `Taxon`.
///
/// # Errors
/// Returns an error if CSV parsing fails.
pub fn parse_taxon_csv(taxa_weights_csv: String) -> Result<Vec<Taxon>, Box<dyn std::error::Error>> {
    let mut rdr = ReaderBuilder::new()
        .has_headers(true)
        .from_reader(taxa_weights_csv.as_bytes());
    
    let mut taxa_weights = Vec::new();
    for record in rdr.deserialize() {
        let row: Taxon = record?;
        taxa_weights.push(row);
    }

    Ok(taxa_weights)
}


/// Serializes a list of `Taxon` into CSV format.
///
/// # Arguments
/// * `taxa` - List of taxa to serialize.
///
/// # Returns
/// CSV as a string.
///
/// # Errors
/// Returns an error if serialization fails.
pub fn generate_taxa_cluster_csv(taxa: Vec<Taxon>) -> Result<String, Box<dyn std::error::Error>> {
    let mut wtr = WriterBuilder::new().from_writer(vec![]);

    for taxon in &taxa {
        wtr.serialize(taxon)?;
    }

    let data = String::from_utf8(wtr.into_inner()?)?;

    Ok(data)
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
pub fn cluster_taxa(graph_xml: String, taxa_weights_csv: String, similarity_threshold: f32) -> Result<String, Box<dyn std::error::Error>> {

    let graph = CTFactorGraph::from_graphml(&graph_xml)?;
    let taxa_weights = parse_taxon_csv(taxa_weights_csv)?;

    let peptidome_dict = get_peptides_per_taxon(&graph)?;
    let (similarities, taxon_index) = compute_detected_peptidome_similarity(peptidome_dict);

    let taxa_weights: Vec<Taxon> = taxa_weights
        .into_iter()
        .filter(|tw| taxon_index.contains_key(&tw.higher_taxa))
        .collect();

    let higher_taxa: Vec<i32> = taxa_weights.iter().map(|tw| tw.higher_taxa).collect();
    let mut weight_sorted_taxa: Vec<i32> = higher_taxa.clone();
    let mut taxa_clusters: Vec<Vec<i32>> = Vec::new();

    let mut cluster_heads: Vec<i32> = Vec::new();

    while ! weight_sorted_taxa.is_empty() {
        let taxon1 = weight_sorted_taxa[0];
        let mut cluster_list: Vec<i32> = Vec::new();
        cluster_heads.push(taxon1);

        for &taxon2 in &higher_taxa {
            if similarities[taxon_index[&taxon2] as usize][taxon_index[&taxon1] as usize] > similarity_threshold {
                cluster_list.push(taxon2);
                if weight_sorted_taxa.contains(&taxon2) {
                    weight_sorted_taxa.retain(|&taxon| taxon != taxon2);
                }
            }
        }

        taxa_clusters.push(cluster_list);
    }

    let mut cluster_weight_sorted_taxa: Vec<Taxon> = taxa_weights.into_iter()
        .filter(|tw| cluster_heads.contains(&tw.higher_taxa))
        .collect();
    // TODO: should we also add rows of taxa_weights where higher_taxa appear in taxa_clusters? 
    // This doesn't seem to do anything in the python code (Bug probably)
    for (taxon, cluster_members) in cluster_weight_sorted_taxa.iter_mut().zip(taxa_clusters.iter()) {
        taxon.cluster_members = cluster_members.clone();
    }

    Ok(generate_taxa_cluster_csv(cluster_weight_sorted_taxa)?)
}


/// Builds a dictionary of peptides per taxon.
///
/// # Arguments
/// * `graph` - Factor graph reference.
///
/// # Returns
/// Map from taxon ID to peptide set.
///
/// # Errors
/// Returns an error if node parsing fails.
fn get_peptides_per_taxon(graph: &CTFactorGraph) -> Result<HashMap<i32, HashSet<i32>>, Box<dyn std::error::Error>> {
    let mut peptidome_dict = HashMap::new();

    for node in graph.get_nodes() {
        if node.is_taxon_node() {
            let node_id: i32 = String::from(node.get_name()).parse()?;
            let neighbors: HashSet<i32> = graph.get_neighbors(node)
                .map(|factor_id| graph.get_peptide_for_factor(factor_id))
                .collect::<Result<_, _>>()?;
            peptidome_dict.insert(node_id, neighbors);
        }
    }

    Ok(peptidome_dict)
}


/// Computes similarity matrix and taxon index map.
///
/// # Arguments
/// * `peptidome_dict` - Map of taxon to peptides.
///
/// # Returns
/// Tuple of (similarity matrix, taxon index map).
fn compute_detected_peptidome_similarity(peptidome_dict: HashMap<i32, HashSet<i32>>) -> (Vec<Vec<f32>>, HashMap<i32, u32>) {
    let mut sim_matrix = Vec::new();
    let mut taxon_index: HashMap<i32, u32> = HashMap::new();

    let peptidome_keys = peptidome_dict.keys();
    for (index, taxon1) in peptidome_keys.clone().enumerate() {
        taxon_index.insert(*taxon1, index as u32);
        let set1 = &peptidome_dict[taxon1];
        let mut sim_row = Vec::new();
        for taxon2 in peptidome_keys.clone() {
            let set2 = &peptidome_dict[taxon2];
            let shared = set1.intersection(set2).count();
            let sim: f32 = if set2.len() == 0 {
                0.0
            } else {
                shared as f32 / set2.len() as f32
            };

            sim_row.push(sim);
        }
        sim_matrix.push(sim_row);
    }

    (sim_matrix, taxon_index)
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::{HashMap, HashSet};

    #[test]
    fn test_vec_to_string_and_string_to_vec() {
        let values = vec![1, 2, 3];
        let serialized = serde_json::to_string(&values).unwrap();
        assert_eq!(serialized, "[1,2,3]");

        let deserialized: Vec<i32> = serde_json::from_str(&serialized).unwrap();
        assert_eq!(deserialized, values);
    }

    #[test]
    fn test_parse_taxon_csv_and_generate_taxa_cluster_csv() {
        let csv_data = "\
id,higher_taxa,scaled_weight,unique
1,10,0.5,true
2,20,0.8,false
";

        let taxa = parse_taxon_csv(csv_data.to_string()).unwrap();
        assert_eq!(taxa.len(), 2);
        assert_eq!(taxa[0].id, 1);
        assert_eq!(taxa[1].higher_taxa, 20);

        let out_csv = generate_taxa_cluster_csv(taxa).unwrap();
        assert!(out_csv.contains("id,higher_taxa,scaled_weight,unique,cluster_members"));
        assert!(out_csv.contains("10"));
    }

    #[test]
    fn test_compute_detected_peptidome_similarity() {
        let mut peptidome_dict: HashMap<i32, HashSet<i32>> = HashMap::new();
        peptidome_dict.insert(1, HashSet::from([1, 2]));
        peptidome_dict.insert(2, HashSet::from([2, 3]));

        let (sim_matrix, taxon_index) = compute_detected_peptidome_similarity(peptidome_dict);
        println!("{:?}", sim_matrix);

        assert_eq!(sim_matrix.len(), 2);
        assert_eq!(sim_matrix[0].len(), 2);
        assert!(taxon_index.contains_key(&1));
        assert!(taxon_index.contains_key(&2));
        assert!(sim_matrix[0][1] > 0.0);
    }

    #[test]
    fn test_generate_taxa_cluster_csv_roundtrip() {
        let taxa = vec![
            Taxon { id: 1, higher_taxa: 10, scaled_weight: 0.5, unique: true, cluster_members: vec![10, 11] },
            Taxon { id: 2, higher_taxa: 20, scaled_weight: 0.8, unique: false, cluster_members: vec![20] },
        ];

        let csv_string = generate_taxa_cluster_csv(taxa.clone()).unwrap();
        let parsed = parse_taxon_csv(csv_string).unwrap();

        assert_eq!(parsed.len(), 2);
        assert_eq!(parsed[0].cluster_members, vec![10, 11]);
        assert_eq!(parsed[1].cluster_members, vec![20]);
    }
}
