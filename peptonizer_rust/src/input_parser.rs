use std::collections::HashMap;
use serde_json;

fn parse_peptides(tsv_content: String) -> Result<(HashMap<String, f64>, HashMap<String, u32>), Box<dyn std::error::Error>> {
    let mut peptides_scores: HashMap<String, f64> = HashMap::new();
    let mut peptides_counts: HashMap<String, u32> = HashMap::new();

    let mut lines = tsv_content.lines().map(|l| l.trim()).filter(|l| !l.is_empty());

    // Skip the header
    lines.next();

    for line in lines {
        let mut parts = line.split('\t');
        let peptide = parts.next()
            .ok_or_else(|| format!("Invalid line (missing peptide): {}", line))?
            .to_string();
        let score_str = parts.next()
            .ok_or_else(|| format!("Invalid line (missing score): {}", line))?;
        let score: f64 = score_str.parse().map_err(|_| {
            format!("Invalid line (score not a number): {} (score={})", line, score_str)
        })?;

        // Update counts
        let count = peptides_counts.entry(peptide.clone()).or_insert(0);
        *count += 1;

        // Update max score
        let entry = peptides_scores.entry(peptide).or_insert(f64::NEG_INFINITY);
        if score > *entry {
            *entry = score;
        }
    }

    Ok((peptides_scores, peptides_counts))
}

pub fn parse_input_peptides(tsv_content: String) -> Result<(String, String), Box<dyn std::error::Error>> {
    let (peptides_scores, peptides_counts) = parse_peptides(tsv_content)?;

    // Convert HashMaps to JSON strings
    let scores_json = serde_json::to_string(&peptides_scores)?;
    let counts_json = serde_json::to_string(&peptides_counts)?;

    Ok((scores_json, counts_json))
}

pub fn parse_unique_peptides(tsv_content: String) -> Result<String, Box<dyn std::error::Error>> {
    let (peptides_scores, peptides_counts) = parse_peptides(tsv_content)?;
    
    let peptides: Vec<String> = peptides_scores.keys().cloned().collect();
    let peptides_json = serde_json::to_string(&peptides)?;

    Ok(peptides_json)
}