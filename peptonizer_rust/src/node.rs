use std::collections::HashMap;
use minidom::Element;
use serde::Serialize;
use std::fmt::Write;

#[derive(Debug, Serialize, Clone)]
pub struct Factor {
    // represents noisy OR cpds, has dimension n(parensports)xn(peptide states(=2))
    pub array: Vec<[f64; 2]>,
    pub array_labels: Vec<String>
}

#[derive(Debug, Serialize, Clone)]
pub enum NodeType {
    PeptideNode { initial_belief_0: f64, initial_belief_1: f64 },
    FactorNode { parent_number: i32, initial_belief: Factor },
    TaxonNode { initial_belief_0: f64, initial_belief_1: f64 },
    ConvolutionTreeNode { number_of_parents: i32 }
}

#[derive(Debug, Clone)]
pub struct Node {
    id: i32,
    name: String,
    incident_edges: Vec<i32>,
    subtype: NodeType
}

impl Node {

    pub fn new(id: i32, name: String, subtype: NodeType) -> Self {
        let incident_edges: Vec<i32> = Vec::new();
     
        Self { id, name, incident_edges, subtype }
    }

    pub fn copy_with_id(&self, new_id: i32) -> Self {
        let mut copy: Node = self.clone();
        copy.id = new_id;
        copy
    }

    pub fn new_convolution_node(id: i32, name: String, number_of_parents: i32) -> Self {
        Self { id, name, incident_edges: Vec::new(), subtype: NodeType::ConvolutionTreeNode { number_of_parents } }
    }

    pub fn add_incident_edge(&mut self, edge: i32) {
        self.incident_edges.push(edge);
    }

    pub fn get_name(&self) -> &str { 
        &self.name
    }

    pub fn get_id(&self) -> i32 {
        self.id
    }

    pub fn get_subtype(&self) -> &NodeType {
        &self.subtype
    }

    pub fn set_subtype(&mut self, subtype: NodeType) {
        self.subtype = subtype;
    }

    pub fn neighbors_count(&self) -> usize {
        self.incident_edges.len()
    }

    pub fn get_incident_edge(&self, neighbor_id: i32) -> i32 {
        self.incident_edges[neighbor_id as usize]
    }

    pub fn get_incident_edges(&self) -> &Vec<i32> {
        &self.incident_edges
    }

    pub fn set_incident_edges(&mut self, new_incident_edges: Vec<i32>) {
        self.incident_edges = new_incident_edges;
    }

    pub fn category(&self) -> &str {
        match self.subtype {
            NodeType::PeptideNode { .. } => "peptide",
            NodeType::FactorNode { .. } => "factor",
            NodeType::TaxonNode { .. } => "taxon",
            NodeType::ConvolutionTreeNode { .. } => "convolution_tree"
        }
    }

    pub fn is_factor_node(&self) -> bool {
        matches!(self.subtype, NodeType::FactorNode { .. })
    }

    pub fn is_taxon_node(&self) -> bool {
        matches!(self.subtype, NodeType::TaxonNode { .. })
    }

    pub fn fill_in_prior(&mut self, prior: f64) {
        if matches!(self.subtype, NodeType::TaxonNode { .. }) {
            self.subtype = NodeType::TaxonNode { initial_belief_0: 1.0 - prior, initial_belief_1: prior };
        }
    }

    pub fn fill_in_factor(&mut self, alpha: f64, beta: f64, regularized: bool) {
        if let NodeType::FactorNode { parent_number, .. } = self.subtype {
            let degree: i32 = parent_number;

            let mut cpd_array: Vec<[f64; 2]> = Vec::with_capacity(degree as usize + 1);
            let mut cpd_array_regularized = cpd_array.clone();
            let exponent_array: Vec<i32> = (0..=degree).collect();
            let divide_array: Vec<f64> = std::iter::once(1i32).chain(1..=degree).map(|x| x as f64).collect();
            
            // regularize cpd priors to penalize higher number of parents
            // log domain to avoid underflow
            let mut cpd_sum: f64 = 0.0;
            let mut cpd_regularized_sum: f64 = 0.0;
            for (i, exp) in exponent_array.iter().enumerate() {
                let cpd_0 = (1.0 - alpha).powi(*exp) * (1.0 - beta);
                let cpd_1 = 1.0 - cpd_0;
                cpd_sum += cpd_0 + cpd_1;
                cpd_array.push([cpd_0, cpd_1]);

                let cpd_regularized_0 = (cpd_0.powi(*exp) * (1.0 - beta)) / divide_array[i];
                let cpd_regularized_1 = 1.0 - cpd_regularized_0;
                cpd_regularized_sum += cpd_regularized_0 + cpd_regularized_1;
                cpd_array_regularized.push([cpd_regularized_0, cpd_regularized_1]);
            }

            // Normalize arrays (assuming normalize and avoid_underflow are implemented)
            Self::normalize_cpd(&mut cpd_array, cpd_sum, false);
            Self::normalize_cpd(&mut cpd_array_regularized, cpd_regularized_sum, true);
            
            // Create factor
            let factor_to_add = if regularized {
                Factor {
                    array: cpd_array_regularized,
                    array_labels: vec![format!("placeholder"), format!("{}0", self.id), format!("{}1", self.id)],
                }
            } else {
                Factor {
                    array: cpd_array,
                    array_labels : vec![format!("placeholder"), format!("{}0", self.id), format!("{}1", self.id)],
                }
            };
            
            // Add factor to the node's attributes
            self.subtype = NodeType::FactorNode { parent_number: parent_number, initial_belief: factor_to_add };
        }
    }

    fn normalize_cpd(arr: &mut Vec<[f64; 2]>, sum: f64, avoid_underflow: bool) {
        for cpd in arr.iter_mut() {
            cpd[0] /= sum;
            cpd[1] /= sum;
    
            if avoid_underflow {
                if cpd[0] < 1e-30 {
                    cpd[0] = 1e-30;
                }
                if cpd[1] < 1e-30 {
                    cpd[1] = 1e-30;
                }
            }
        }
    }

    fn parse_data(data: &Element) -> (String, String) {
        let key = data.attr("key").unwrap().to_string();
        let value = data.text();
    
        (key, value)
    }

    pub fn parse_node(node: &Element, id: i32) -> Result<Self, String> {
        // Process a node
        let name = node.attr("id").unwrap().to_string();
    
        // Initialize data for this node
        let mut current_node_data = HashMap::new();
        for data in node.children().filter(|d| d.name() == "data") {
            let (data_key, data_val) = Self::parse_data(data);
            current_node_data.insert(data_key, data_val);
        }
    
        let subtype: NodeType = match current_node_data.get("d2").map(String::as_str) {
            Some("factor") => {
                let parent_number: i32 = current_node_data.get("d3").unwrap().parse().unwrap();
                let initial_belief: Factor = Factor { array: Vec::new(), array_labels: Vec::new() };
                NodeType::FactorNode { parent_number, initial_belief }
            }
            Some("peptide") => {
                let initial_belief_0: f64 = current_node_data.get("d0").unwrap().parse().unwrap();
                let initial_belief_1: f64 = current_node_data.get("d1").unwrap().parse().unwrap();
                NodeType::PeptideNode { initial_belief_0, initial_belief_1 }
            }
            Some("taxon") => {
                let initial_belief_0: f64 = 0.0;
                let initial_belief_1: f64 = 0.0;
                NodeType::TaxonNode { initial_belief_0, initial_belief_1 }
            }
            _ => {
                return Err("Node data has unknown type".to_string());
            }
        };

        Ok(Self { id, name, incident_edges: Vec::new(), subtype })
    }

    pub fn to_graphml(&self) -> String {

        let mut graphml = String::new();

        writeln!(&mut graphml, r#"<node id="{}">"#, self.name).unwrap();
        writeln!(&mut graphml, r#"<data key="d2">{}</data>"#, self.category()).unwrap();

        match &self.subtype {
            NodeType::PeptideNode { initial_belief_0, initial_belief_1 } => {
                writeln!(&mut graphml, r#"<data key="d0">{}</data>"#, initial_belief_0).unwrap();
                writeln!(&mut graphml, r#"<data key="d1">{}</data>"#, initial_belief_1).unwrap();
            },
            NodeType::FactorNode { parent_number, initial_belief: _ } => {
                writeln!(&mut graphml, r#"<data key="d3">{}</data>"#, parent_number).unwrap();
            },
            _ => {}
        }

        writeln!(&mut graphml, r#"</node>"#).unwrap();

        graphml
    }
}
