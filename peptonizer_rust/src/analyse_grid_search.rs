use std::collections::{HashSet, HashMap};
use crate::taxa_clustering::{Taxon, parse_taxon_csv};
use crate::utils::log;

pub fn compute_goodness(
    clustered_taxa_weights_csv: String, 
    peptonizer_results: String
) -> Result<f64, Box<dyn std::error::Error>> {

    let taxid_weights: Vec<Taxon> = parse_taxon_csv(clustered_taxa_weights_csv)?;
    let higher_taxa: Vec<i32> = taxid_weights.iter().map(|t| t.higher_taxa).collect();

    let taxa_scores: HashMap<String, f64> = serde_json::from_str(&peptonizer_results)?;
    let mut taxa_scores: Vec<(&String, &f64)> = taxa_scores.iter().collect();

    taxa_scores.sort_by(|a, b| b.1.partial_cmp(a.1).unwrap()); // ascending order

    let sorted_ids: Vec<i32> = taxa_scores.iter()
        .map(|(k, _)| (*k).clone().parse::<i32>())
        .collect::<Result<Vec<_>, _>>()?;
    let sorted_scores: Vec<f64> = taxa_scores.iter().map(|(_, v)| **v).collect::<Vec<_>>();

    let entropy = entropy(&sorted_scores);
    let rbo = rbo(&higher_taxa, &sorted_ids);

    Ok(rbo / entropy.powi(2))
}


fn rbo(list1: &[i32], list2: &[i32]) -> f64 {
    let k = list1.len().min(list2.len());
    let mut sum = 0.0;

    let mut seen1 = HashSet::new();
    let mut seen2 = HashSet::new();

    for d in 1..=k {
        seen1.insert(list1[d - 1]);
        seen2.insert(list2[d - 1]);

        let overlap = seen1.intersection(&seen2).count() as f64;
        let agreement = overlap / d as f64;

        sum += agreement;
    }

    sum / k as f64
}


fn entropy(values: &[f64]) -> f64 {
    let sum: f64 = values.iter().sum();

    if sum == 0.0 {
        return 0.0;
    }

    values.iter()
        .map(|&v| {
            let p = v / sum;
            if p > 0.0 {
                -p * p.log2()
            } else {
                0.0
            }
        })
        .sum()
}
