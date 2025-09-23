use std::collections::{HashMap, HashSet};
use crate::utils::log;
use crate::random::select_random_samples_with_weights;
use csv::Writer;
use crate::unipept_communicator::get_unique_lineage_at_specified_rank;

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
pub fn perform_taxa_weighing(
    pep_taxa: String,
    pep_scores: String,
    pep_psm_counts: String,
    max_taxa: usize,
    taxa_rank: String
) -> (String, String) {
    log("Parsing Unipept responses from disk...");
    let pep_taxa: HashMap<String, Vec<i32>> = serde_json::from_str(&pep_taxa).unwrap();

    let sequences: Vec<String> = pep_taxa.iter().map(|(seq, _taxa)| seq.to_owned()).collect();

    let mut taxa: Vec<Vec<i32>> = pep_taxa.into_iter().map(|(_seq, taxa)| taxa).collect();
    
    log("Started mapping all taxon ids to the specified rank...");
    normalize_unipept_responses(&mut taxa, &taxa_rank);
    // TODO: enable sampling (disabled to get deterministic output)
    // let chosen_idx: HashSet<usize> = weighted_random_sample(&taxa, 10000); // TODO: what if < 10000 or close to 10000 + hardcoded not a good idea
    let chosen_idx: Vec<usize> = (0..taxa.len()).collect();

    log(&format!("Using {} sequences as input...", chosen_idx.len()));

    log("Normalizing peptides and converting to vectors...");

    let sequences: Vec<String> = chosen_idx.iter().map(|idx| sequences[*idx].to_owned()).collect();
    let taxa: Vec<Vec<i32>> = chosen_idx.iter().map(|idx| taxa[*idx].to_owned()).collect();

    // Parse scores from JSON string to hashmap, only keep the randomly selected samples.
    let pep_scores_map: HashMap<String, f32> = serde_json::from_str(&pep_scores).unwrap();
    let mut pep_scores: Vec<f32> = vec![0.0; sequences.len()];
    for i in 0..sequences.len() {
        pep_scores[i] = pep_scores_map[&sequences[i]];
    }

    // parse counts from JSON string to hashmap, only keep the randomly selected samples.
    let pep_psm_counts_map: HashMap<String, i32> = serde_json::from_str(&pep_psm_counts).unwrap();
    let mut pep_psm_counts: Vec<i32> = vec![0; sequences.len()];
    for i in 0..sequences.len() {
        pep_psm_counts[i] = pep_psm_counts_map[&sequences[i]];
    }

    /* Score the degeneracy of a taxa, i.e.,
       how conserved a peptide sequence is between taxa.
       map all taxids in the list in the taxa column back to their taxid at species level (or the rank specified by the user)
       Right now, HigherTaxa is simply a copy of taxa. This step still needs to be optimized.
       Move taxa to highertaxa because taxa is not used anymore.
    */
    let higher_taxa: Vec<Vec<i32>> = taxa; 

    // Divide the number of PSMs of a peptide by the number of taxa the peptide is associated with, exponentiated by 3
    log("Started dividing the number of PSMS of a peptide by the number the peptide is associated with...");
    let weights: Vec<f32> = pep_psm_counts.iter()
                                          .zip(higher_taxa.iter().map(|taxa| taxa.len().pow(3)))
                                          .map(|(&count, len_cube)| count as f32 / len_cube as f32)
                                          .collect();

    let unique_psm_taxa: HashSet<i32> = higher_taxa.iter()
                                                   .filter(|tax| tax.len() == 1)
                                                   .map(|tax| tax[0])
                                                   .collect();

    // Sum up the weights of a taxon and sort by weight
    log("Started summing the weights of a taxon and sorting them by weight...");
    let log_weights: Vec<f32> = weights.iter().map(|w| (w + 1.0).log10()).collect();

    //  Since large proteomes tend to have more detectable peptides,
    // we adjust the weight by dividing by the size of the proteome i.e.,
    // the number of proteins that are associated with a taxon
    let scaled_weight = log_weights.clone();

    let mut tax_id_weights: HashMap<i32, f32> = HashMap::new();
    for (ids, weight) in higher_taxa.clone().into_iter().zip(scaled_weight.clone().into_iter()) {
        for id in ids {
            *tax_id_weights.entry(id).or_insert(0.0) += weight;
        }
    }
    let mut sorted_tax_id_weights: Vec<(i32, f32)> = tax_id_weights.into_iter().collect();
    sorted_tax_id_weights.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    let (tax_ids, tax_id_weights): (Vec<i32>, Vec<f32>) = sorted_tax_id_weights.into_iter().unzip();

    // Retrieves the specified taxonomic rank taxid in the lineage of each of the species-level taxids returned by
    // Unipept for both the UnipeptFrame and the TaxIdWeightFrame
    let higher_unique_psm_taxids = unique_psm_taxa;

    // Group the duplicate entries of higher up taxa and sum their weights
    let higher_taxid_weights = tax_id_weights;

    let higher_taxid_unique: Vec<bool> = tax_ids.iter().map(|id| higher_unique_psm_taxids.contains(&id)).collect();

    // TODO: check if duplicate removal is necessary, example datasets do not contain doubles. Link to sampling. Lower the amount of sampled?
    let sequence_csv;
    if higher_taxid_weights.len() < 50 {
        sequence_csv = generate_sequence_csv(None, false, sequences, pep_scores, pep_psm_counts, higher_taxa, weights, log_weights);
    } else {
        let mut taxa_to_include: HashSet<i32> = tax_ids.iter().take(max_taxa).cloned().collect();
        taxa_to_include.extend(higher_unique_psm_taxids);

        sequence_csv = generate_sequence_csv(Some(taxa_to_include), true, sequences, pep_scores, pep_psm_counts, higher_taxa, weights, log_weights);
    }

    let taxa_weights_csv = generate_taxa_weights_csv(tax_ids, higher_taxid_weights, higher_taxid_unique);

    (sequence_csv, taxa_weights_csv)
}

/// Generates a CSV for sequences with associated taxonomic weights and scores.
///
/// # Arguments
///
/// * `taxa_to_include` - Optional set of taxa IDs to filter the output.
/// * `filter_taxa` - Whether to filter sequences based on `taxa_to_include`.
/// * `sequences` - List of peptide sequences.
/// * `scores` - List of peptide scores corresponding to `sequences`.
/// * `psms` - List of PSM counts corresponding to `sequences`.
/// * `higher_taxa` - List of lists of higher taxa IDs for each sequence.
/// * `weights` - Computed weights for each sequence.
/// * `log_weights` - Log-transformed weights for each sequence.
///
/// # Returns
///
/// CSV string containing one row per peptide-taxon pair with columns:
/// "id", "sequence", "score", "psms", "higher_taxa", "weight", "log_weight".
fn generate_sequence_csv(taxa_to_include: Option<HashSet<i32>>, filter_taxa: bool, sequences: Vec<String>, scores: Vec<f32>, psms: Vec<i32>, higher_taxa: Vec<Vec<i32>>, weights: Vec<f32>, log_weights: Vec<f32>) -> String {

    let mut wtr = Writer::from_writer(vec![]);

    let _ = wtr.write_record(&["id", "sequence", "score", "psms", "higher_taxa", "weight", "log_weight"]);

    let mut id = 0;
    for i in 0..sequences.len() {
        for taxon in &higher_taxa[i] {
            if (! filter_taxa) || taxa_to_include.as_ref().unwrap().contains(&taxon) {
                let _ = wtr.write_record(&[
                    id.to_string(),
                    sequences[i].clone(), 
                    scores[i].to_string(), 
                    psms[i].to_string(), 
                    taxon.to_string(), 
                    weights[i].to_string(), 
                    log_weights[i].to_string()
                ]).unwrap();
                id += 1;
            }
        }
    }

    let csv: String = String::from_utf8(wtr.into_inner().unwrap()).unwrap();

    csv
}

/// Generates a CSV of taxa weights.
///
/// # Arguments
///
/// * `higher_taxa` - List of taxon IDs.
/// * `higher_taxid_weights` - List of computed weights corresponding to `higher_taxa`.
/// * `higher_taxid_unique` - List indicating whether each taxon is uniquely associated with a peptide.
///
/// # Returns
///
/// CSV string with columns: "id", "higher_taxa", "scaled_weight", "unique".
fn generate_taxa_weights_csv(higher_taxa: Vec<i32>, higher_taxid_weights: Vec<f32>, higher_taxid_unique: Vec<bool>) -> String {
    let mut wtr = Writer::from_writer(vec![]);

    let _ = wtr.write_record(&["id", "higher_taxa", "scaled_weight", "unique"]);

    for i in 0..higher_taxa.len() {
        let _ = wtr.write_record(&[
            i.to_string(),
            higher_taxa[i].to_string(),
            higher_taxid_weights[i].to_string(),
            higher_taxid_unique[i].to_string()
        ]).unwrap();
    }

    let csv: String = String::from_utf8(wtr.into_inner().unwrap()).unwrap();

    csv
    
}

/// Maps taxa lists onto the taxonomic rank specified by the user.
///
/// # Arguments
///
/// * `taxa` - Mutable reference to a vector of vectors of taxon IDs.
/// * `taxa_rank` - The desired taxonomic rank to normalize to (e.g., "species").
fn normalize_unipept_responses(taxa: &mut Vec<Vec<i32>>, taxa_rank: &str) {
    
    let mut lineage_cache: HashMap<i32, Vec<Option<i32>>> = HashMap::new();

    // Map all taxa onto the rank specified by the user
    for i in 0..taxa.len() {
        taxa[i] = get_unique_lineage_at_specified_rank(&taxa[i], taxa_rank, &mut lineage_cache);
    }
}

/// Selects `n` random indices from `taxa` vectors, weighted by inverse degeneracy.
///
/// # Arguments
///
/// * `taxa` - Vector of vectors of taxon IDs for each peptide.
/// * `n` - Number of random samples to select.
///
/// # Returns
///
/// A `HashSet` of selected indices, chosen with probability proportional to
/// `1 / number_of_taxa_per_peptide`.
fn weighted_random_sample(taxa: &Vec<Vec<i32>>, n: usize) -> HashSet<usize> {
    
    // Calculate normalized weights based on the length of the taxa array
    let weights: Vec<f64> = taxa.iter().map(|taxon| if taxon.len() == 0 { 0.0 } else { 1.0 / taxon.len() as f64 }).collect();
    let total_weight: f64 = weights.iter().sum();
    let normalized_weights: Vec<f64> = weights.iter().map(|w| w / total_weight as f64).collect();

    let samples: HashSet<usize> = select_random_samples_with_weights(normalized_weights, n);

    samples
}
