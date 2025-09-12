use crate::unipept_communicator::{get_taxa_for_peptides, get_descendants_for_taxa};
use std::collections::{HashMap, HashSet};


pub fn fetch_peptides_and_filter_taxa(
    peptides: String,
    rank: String,
    taxon_query: String
) -> String {
    // Parse arguments
    let peptides: Vec<String> = serde_json::from_str(&peptides).unwrap();
    let taxon_query_ids: Vec<i32> = serde_json::from_str(&taxon_query).unwrap();
    
    // First we retrieve all taxa associated with the given peptids
    let mut peptides_taxa: HashMap<String, Vec<i32>> = get_taxa_for_peptides(peptides);

    // Then, we make sure to filter the taxa and only keep those that are associated 
    // to the taxa of interest indicated by the user. Retrieve all (in)direct children
    // of the filter taxa provided by the user
    let taxa_filter: HashSet<i32> = get_descendants_for_taxa(taxon_query_ids, rank);

    // Compute the intersection of the taxa that should be retained and the original list of taxa
    for taxa_list in peptides_taxa.values_mut() {
        taxa_list.retain(|taxon| taxa_filter.contains(taxon));
    }

    serde_json::to_string(&peptides_taxa).unwrap()
}
