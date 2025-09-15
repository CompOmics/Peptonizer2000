use crate::factor_graph::CTFactorGraph;
use std::collections::{HashMap, HashSet};
use serde::{Serialize, Deserialize, Deserializer};
use csv::{ReaderBuilder, WriterBuilder};


#[derive(Deserialize, Serialize, Debug)]
pub struct Taxon {
    pub id: i32,
    pub higher_taxa: i32,
    pub scaled_weight: f32,
    pub unique: bool,

    #[serde(default, serialize_with = "vec_to_string", deserialize_with = "string_to_vec")]
    pub cluster_members: Vec<i32>
}


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

pub fn string_to_vec<'de, D>(deserializer: D) -> Result<Vec<i32>, D::Error>
where
    D: Deserializer<'de>,
{
    let s: &str = Deserialize::deserialize(deserializer)?;
    let vec = s.split(',')
               .filter_map(|item| item.trim().parse::<i32>().ok())
               .collect();
    Ok(vec)
}


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

pub fn generate_taxa_cluster_csv(taxa: Vec<Taxon>) -> Result<String, Box<dyn std::error::Error>> {
    let mut wtr = WriterBuilder::new().from_writer(vec![]);

    for taxon in &taxa {
        wtr.serialize(taxon)?;
    }

    let data = String::from_utf8(wtr.into_inner()?)?;

    Ok(data)
}


pub fn cluster_taxa(graph_xml: String, taxa_weights_csv: String, similarity_threshold: f32) -> Result<String, Box<dyn std::error::Error>> {

    let graph = CTFactorGraph::from_graphml(&graph_xml)?;
    let taxa_weights = parse_taxon_csv(taxa_weights_csv)?;

    let peptidome_dict = get_peptides_per_taxon(&graph);
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


fn get_peptides_per_taxon(graph: &CTFactorGraph) -> HashMap<i32, HashSet<i32>> {
    let mut peptidome_dict = HashMap::new();

    for node in graph.get_nodes() {
        if node.is_taxon_node() {
            let node_id: i32 = String::from(node.get_name()).parse().expect("Taxon node name is no number");
            let neighbors = graph.get_neighbors(node);
            peptidome_dict.insert(node_id, neighbors[..neighbors.len()-4].iter().cloned().collect());
        }
    }

    peptidome_dict
}

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
                (shared / set2.len()) as f32
            };

            sim_row.push(sim);
        }
        sim_matrix.push(sim_row);
    }

    (sim_matrix, taxon_index)
}