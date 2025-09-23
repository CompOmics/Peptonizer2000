use std::collections::{HashMap, HashSet};
use std::cmp::min;
use crate::http_client::*;
use serde::{Serialize, Deserialize};
use serde_json::{ Value };


/// Base URL for the UniPept API
const UNIPEPT_URL: &str = "https://api.unipept.ugent.be";
/// Endpoint for mapping peptides to filtered taxa
const UNIPEPT_PEPT2FILTERED_ENDPOINT: &str = "/api/v2/pept2taxa";
/// Endpoint for retrieving taxonomic lineages
const UNIPEPT_TAXONOMY_ENDPOINT: &str = "/api/v2/taxonomy";

/// Maximum number of peptides per request to the peptide-to-taxa endpoint
const UNIPEPT_PEPTIDES_BATCH_SIZE: usize = 2000;

/// Maximum number of taxa per request to the taxonomy endpoint
const TAXONOMY_ENDPOINT_BATCH_SIZE: usize = 100;


/// Standard NCBI taxonomy ranks for lineage retrieval
const NCBI_RANKS: &[&str] = &[
    "superkingdom",
    "kingdom",
    "subkingdom",
    "superphylum",
    "phylum",
    "subphylum",
    "superclass",
    "class",
    "subclass",
    "superorder",
    "order",
    "suborder",
    "infraorder",
    "superfamily",
    "family",
    "subfamily",
    "tribe",
    "subtribe",
    "genus",
    "subgenus",
    "species_group",
    "species_subgroup",
    "species",
    "subspecies",
    "strain",
    "varietas",
    "forma"
];


/// Payload structure for requesting taxonomy lineages from UniPept
#[derive(Serialize, Deserialize, Debug)]
pub struct HTTPTaxonomyPayload {
    input: Vec<i32>,
    extra: bool
}

/// Payload structure for mapping peptides to taxa
#[derive(Serialize, Deserialize, Debug)]
pub struct HTTPPept2TaxaPayload { 
    input: Vec<String>,
    compact: bool,
    tryptic: bool
}

/// Response structure for peptide-to-taxa mapping
#[derive(Serialize, Deserialize, Debug)]
pub struct HTTPPept2TaxaResponse {
    peptide: String,
    taxa: Vec<i32>
}

/// Payload structure for retrieving descendants of taxa at specified ranks
#[derive(Serialize, Deserialize, Debug)]
pub struct HTTPTaxonomyDescendantsPayload {
    input: Vec<i32>,
    descendants: bool,
    descendants_ranks: Vec<String>
}

/// Response structure for retrieving descendants of a taxon
#[derive(Serialize, Deserialize, Debug)]
pub struct HTTPTaxonomyDescendantsResponse {
    taxon_id: i32,
    taxon_name: String,
    taxon_rank: String,
    descendants: Vec<i32>
}

/// Represents a response from the UniPept taxonomy API
#[derive(Serialize, Deserialize, Debug)]
struct TaxonomyResponse {
    taxon_id: i32,
    taxon_name: String,
}


/// Parses a JSON string returned by the UniPept API into a vector of key-value maps
///
/// # Arguments
/// * `http_response` - JSON string from UniPept API
///
/// # Returns
/// Vector of hash maps, where each key maps to an `Option<i32>` value. Only numeric or null values are retained.
fn parse_response_json_string(http_response: &str) -> Vec<HashMap<String, Option<i32>>> {
    let http_response_map: Vec<HashMap<String, Option<i32>>> = serde_json::from_str::<Vec<HashMap<String, Value>>>(http_response)
            .unwrap()
            .into_iter()
            .map(|mut obj: HashMap<String, Value>| {
                // Remove the key-value pair where the value is a string
                obj.retain(|_, v| v.is_null() || v.is_number());

                // Convert the remaining keys and values to `HashMap<String, Option<i32>>`
                obj.into_iter()
                    .map(|(key, value)| {
                        let value = if value.is_number() {
                            Some(value.as_i64().unwrap() as i32)
                        } else {
                            None
                        };
                        (key, value)
                    })
                    .collect()
            })
            .collect();
    
    http_response_map
}

/// Retrieves the unique lineage taxa IDs at a specified taxonomic rank.
///
/// This function queries the UniPept taxonomy API for the given `target_taxa` and extracts
/// the taxon IDs at the specified `taxa_rank`. To minimize API requests, it uses a cache
/// (`lineage_cache`) to store previously fetched lineages. 
///
/// # Arguments
///
/// * `target_taxa` - A reference to a vector of taxon IDs for which the lineage is requested.
/// * `taxa_rank` - The target taxonomic rank (e.g., "species", "genus") at which the unique lineage is extracted.
/// * `lineage_cache` - A mutable reference to a hash map that stores previously fetched lineages.
///
/// # Returns
///
/// A vector of unique taxon IDs corresponding to the specified `taxa_rank`.
///
/// # Panics
///
/// The function will panic if:
/// - The `taxa_rank` does not exist in the predefined `NCBI_RANKS`.
pub fn get_unique_lineage_at_specified_rank(target_taxa: &Vec<i32>, taxa_rank: &str, lineage_cache: &mut HashMap<i32, Vec<Option<i32>>>) -> Vec<i32> {

    let url: String = [UNIPEPT_URL, UNIPEPT_TAXONOMY_ENDPOINT].concat();

    // Remove duplicates from input
    let target_taxa: HashSet<i32> = target_taxa.iter().cloned().collect();

    // Prepare a list of taxa that are not yet in the cache
    let taxa_to_request: Vec<i32> = target_taxa.iter().filter(| tax_id | ! lineage_cache.contains_key(tax_id)).cloned().collect();

    let http_client = &create_http_client();
    // Fetch lineages from the API for taxa not in the cache
    for i in (0..taxa_to_request.len()).step_by(TAXONOMY_ENDPOINT_BATCH_SIZE) {

        let batch_size: usize = std::cmp::min(TAXONOMY_ENDPOINT_BATCH_SIZE, taxa_to_request.len() - i);
        let batch: Vec<i32> = taxa_to_request[i..(i + batch_size)].to_vec();
        let payload = HTTPTaxonomyPayload { input: batch, extra: true };

        // Perform the HTTP POST request
        let http_response:  String = 
            http_client.perform_post_request(url.clone(), &payload)
            .map_err(|e| format!("Failed to retrieve taxonomy data for batch {}. Error message: {}", (i / TAXONOMY_ENDPOINT_BATCH_SIZE), e))
            .unwrap();

        let http_response = parse_response_json_string(&http_response);

        for lineage_json in &http_response {
            let lineage: Vec<Option<i32>> = NCBI_RANKS.iter()
                    .filter_map(|key| lineage_json.get(&format!("{}_id", key)).cloned())
                    .collect();
            let taxon_id: i32 = lineage_json.get("taxon_id").unwrap().unwrap();
            lineage_cache.insert(taxon_id, lineage);
        }
        
    }

    let rank_idx = NCBI_RANKS.iter().position(|&ncbi_rank| ncbi_rank == taxa_rank).unwrap();
    let lineage: HashSet<i32> = target_taxa.iter()
                                            .filter_map(|taxon| lineage_cache[&taxon][rank_idx].clone())
                                            .collect();
    let lineage: Vec<i32> = lineage.into_iter().collect();

    lineage
}

/// Queries Unipept and returns all the taxa that are associated with the given list of peptides.
/// 
/// For each peptide in the input, an entry in the output map is created, which points to the
/// taxon IDs associated with this peptide.
/// 
/// # Arguments
/// 
/// * `peptides` - A list of peptide sequences for which all associated taxa should be queried.
/// 
/// # Errors
/// 
/// Returns an error if the Unipept API server responds with an error or if a network issue occurs.
/// 
/// # Returns
/// 
/// A map from each peptide in the input list to its associated taxa IDs.
pub fn get_taxa_for_peptides(peptides: Vec<String>) -> HashMap<String, Vec<i32>> {
    
    let url = [UNIPEPT_URL, UNIPEPT_PEPT2FILTERED_ENDPOINT].concat();

    let mut output = HashMap::new();

    let http_client = &create_http_client();
    // Split the peptides into batches of a predefined size
    for i in (0..peptides.len()).step_by(UNIPEPT_PEPTIDES_BATCH_SIZE) {

        let end_batch = min(i+UNIPEPT_PEPTIDES_BATCH_SIZE, peptides.len());
        let batch = peptides[i..end_batch].to_vec();
        let payload = HTTPPept2TaxaPayload { input: batch, compact: true, tryptic: true };

        let http_response:  String = http_client.perform_post_request(url.clone(), &payload)
            .map_err(|e| format!("Failed to retrieve taxa data for batch {}. Error message: {}", (i / UNIPEPT_PEPTIDES_BATCH_SIZE), e))
            .unwrap();

        let http_response = serde_json::from_str::<Vec<HTTPPept2TaxaResponse>>(&http_response).unwrap();

        for peptide_data in &http_response {
            let original_taxa: Vec<i32> = peptide_data.taxa.clone();
            output.insert(peptide_data.peptide.clone(), original_taxa);
        }
    }

    output
}


/// Returns a list of all taxon IDs that are descendants of the given taxa in `target_taxa`.
///
/// # Arguments
///
/// * `target_taxa` - A list of taxon IDs for which all descendants at a specific NCBI rank (and lower) should be retrieved.
/// * `descendants_rank` - The maximum rank that each of the descendants should have in the NCBI taxonomy.
///   All descendants that are defined at this rank or deeper are reported.
///
/// # Errors
///
/// Returns an error if the Unipept API server responds with an error, or if something goes wrong
/// with the network.
///
/// # Returns
///
/// A list of taxon IDs that meet the given rank criteria.
pub fn get_descendants_for_taxa(target_taxa: Vec<i32>, descendant_rank: String) -> HashSet<i32> {
    
    let url = [UNIPEPT_URL, UNIPEPT_TAXONOMY_ENDPOINT].concat();
    let mut all_descendants = HashSet::new();

    // We need to get all children at the requested level, AND at lower levels. That's what we are using the ranks array for.
    let rank_idx = NCBI_RANKS.iter().position(|&ncbi_rank| ncbi_rank == descendant_rank).unwrap();
    let descentants_ranks: Vec<String> = NCBI_RANKS[rank_idx..].iter().map(|&s| s.to_string()).collect();

    let http_client = &create_http_client();
    // Split the target taxa into batches of 15
    for i in (0..target_taxa.len()).step_by(TAXONOMY_ENDPOINT_BATCH_SIZE) {

        let end_batch = min(i+TAXONOMY_ENDPOINT_BATCH_SIZE, target_taxa.len());
        let batch = target_taxa[i..end_batch].to_vec();
        let payload = HTTPTaxonomyDescendantsPayload { input: batch, descendants: true, descendants_ranks: descentants_ranks.clone() };

        // Perform the HTTP Post request
        let http_response = http_client.perform_post_request(url.clone(), &payload)
            .map_err(|e| format!("Failed to retrieve taxonomy data for batch {}. Error message: {}", (i / TAXONOMY_ENDPOINT_BATCH_SIZE), e))
            .unwrap();
        let http_response = serde_json::from_str::<Vec<HTTPTaxonomyDescendantsResponse>>(&http_response).unwrap();
        
        for response in http_response {
            all_descendants.extend(response.descendants);
        }
    }

    all_descendants
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
/// A `HashMap<i32, String>` mapping taxon IDs to their corresponding taxon names.
pub fn get_names_for_taxa(target_taxa: &Vec<i32>) -> Result<HashMap<i32, String>, String> {
    let url = format!("{}{}", UNIPEPT_URL, UNIPEPT_TAXONOMY_ENDPOINT);
    let mut output: HashMap<i32, String> = HashMap::new();

    let http_client = &create_http_client();

    for i in (0..target_taxa.len()).step_by(TAXONOMY_ENDPOINT_BATCH_SIZE) {
        let batch: Vec<i32> = target_taxa[i..std::cmp::min(i + TAXONOMY_ENDPOINT_BATCH_SIZE, target_taxa.len())]
            .to_vec();

        let payload = serde_json::json!({
            "input": batch
        });

        // Perform the HTTP POST request
        let http_response = http_client.perform_post_request(url.clone(), &payload)
            .map_err(|e| format!("Communication error: {}", e)).unwrap();

        let http_response = serde_json::from_str::<Vec<TaxonomyResponse>>(&http_response).unwrap();
        
        for response in http_response {
            output.insert(response.taxon_id, response.taxon_name);
        }
    }

    Ok(output)
}