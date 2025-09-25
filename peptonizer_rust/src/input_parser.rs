use std::collections::HashMap;
use serde_json;


/// Parses peptide data from a TSV string.
///
/// The TSV is expected to have a header row and contain at least two columns:
/// peptide sequence and score. For each peptide, the maximum score is stored,
/// along with the number of times the peptide appears.
///
/// # Arguments
/// * `tsv_content` - Input TSV string with peptide data.
///
/// # Returns
/// A tuple containing:
/// * `HashMap<String, f64>` - Maximum score per peptide.
/// * `HashMap<String, u32>` - Count of occurrences per peptide.
///
/// # Errors
/// Returns an error if lines are malformed, scores cannot be parsed,
/// or the input is otherwise invalid.
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
pub fn parse_input_peptides(tsv_content: String) -> Result<(String, String), Box<dyn std::error::Error>> {
    let (peptides_scores, peptides_counts) = parse_peptides(tsv_content)?;

    // Convert HashMaps to JSON strings
    let scores_json = serde_json::to_string(&peptides_scores)?;
    let counts_json = serde_json::to_string(&peptides_counts)?;

    Ok((scores_json, counts_json))
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
pub fn parse_unique_peptides(tsv_content: String) -> Result<String, Box<dyn std::error::Error>> {
    let (peptides_scores, peptides_counts) = parse_peptides(tsv_content)?;
    
    let peptides: Vec<String> = peptides_scores.keys().cloned().collect();
    let peptides_json = serde_json::to_string(&peptides)?;

    Ok(peptides_json)
}


#[cfg(test)]
mod tests {
    use super::*;

    fn sample_tsv() -> String {
        "peptide\tscore\nPEP1\t0.5\nPEP2\t1.0\nPEP1\t0.8\n".to_string()
    }

    #[test]
    fn test_parse_peptides_basic() {
        let (scores, counts) = parse_peptides(sample_tsv()).unwrap();

        assert_eq!(scores.len(), 2);
        assert_eq!(counts.len(), 2);

        // PEP1 appears twice, max score = 0.8
        assert_eq!(counts.get("PEP1"), Some(&2));
        assert!((scores.get("PEP1").unwrap() - 0.8).abs() < 1e-10);

        // PEP2 appears once, score = 1.0
        assert_eq!(counts.get("PEP2"), Some(&1));
        assert!((scores.get("PEP2").unwrap() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_parse_input_peptides() {
        let (scores_json, counts_json) = parse_input_peptides(sample_tsv()).unwrap();

        let scores: HashMap<String, f64> = serde_json::from_str(&scores_json).unwrap();
        let counts: HashMap<String, u32> = serde_json::from_str(&counts_json).unwrap();

        assert_eq!(scores["PEP1"], 0.8);
        assert_eq!(scores["PEP2"], 1.0);
        assert_eq!(counts["PEP1"], 2);
        assert_eq!(counts["PEP2"], 1);
    }

    #[test]
    fn test_parse_unique_peptides() {
        let peptides_json = parse_unique_peptides(sample_tsv()).unwrap();
        let peptides: Vec<String> = serde_json::from_str(&peptides_json).unwrap();

        assert_eq!(peptides.len(), 2);
        assert!(peptides.contains(&"PEP1".to_string()));
        assert!(peptides.contains(&"PEP2".to_string()));
    }
}
