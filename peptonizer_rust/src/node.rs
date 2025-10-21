use std::collections::HashMap;
use minidom::Element;
use serde::Serialize;
use std::fmt::Write;

/// Represents a noisy-OR conditional probability distribution (CPD) for factor nodes.
#[derive(Debug, Serialize, Clone)]
pub struct Factor {
    /// CPD array with dimensions n(parents) × 2 (inactive/active).
    pub array: Vec<[f64; 2]>,
    /// Labels for the CPD entries.
    pub array_labels: Vec<String>
}


/// Defines the type of node in the factor graph with its initial beliefs.
#[derive(Debug, Serialize, Clone)]
pub enum NodeType {
    /// Peptide node with prior probabilities.
    PeptideNode { initial_belief_0: f64, initial_belief_1: f64 },
    /// Factor node with parent count and CPD.
    FactorNode { parent_number: u32, initial_belief: Factor },
    /// Taxon node with prior probabilities.
    TaxonNode { initial_belief_0: f64, initial_belief_1: f64 },
    /// Convolution tree node with a number of parents.
    ConvolutionTreeNode { number_of_parents: u32 }
}


/// Represents a node in the factor graph with its attributes and connections.
#[derive(Debug, Clone)]
pub struct Node {
    id: u32,
    name: String,
    incident_edges: Vec<u32>,
    subtype: NodeType
}


impl Node {

    /// Creates a new node with given ID, name, and subtype.
    ///
    /// # Arguments
    /// * `id` - Unique node identifier.
    /// * `name` - Node label.
    /// * `subtype` - Type of the node with initial state.
    ///
    /// # Returns
    /// A new `Node`.
    pub fn new(id: usize, name: String, subtype: NodeType) -> Self {
        let incident_edges: Vec<u32> = Vec::new();
     
        Self { id: id as u32, name, incident_edges, subtype }
    }

    /// Creates a copy of the node with a new ID.
    ///
    /// # Arguments
    /// * `new_id` - Replacement node ID.
    ///
    /// # Returns
    /// A cloned node with updated ID.
    pub fn copy_with_id(&self, new_id: usize) -> Self {
        let mut copy: Node = self.clone();
        copy.id = new_id as u32;
        copy
    }

    /// Creates a new convolution tree node.
    ///
    /// # Arguments
    /// * `id` - Node ID.
    /// * `name` - Node label.
    /// * `number_of_parents` - Number of parents in convolution tree.
    ///
    /// # Returns
    /// A new convolution tree node.
    pub fn new_convolution_node(id: usize, name: String, number_of_parents: usize) -> Self {
        Self { id: id as u32, name, incident_edges: Vec::new(), subtype: NodeType::ConvolutionTreeNode { number_of_parents: number_of_parents as u32 } }
    }

    /// Adds an incident edge to the node.
    ///
    /// # Arguments
    /// * `edge` - Edge identifier.
    pub fn add_incident_edge(&mut self, edge: usize) {
        self.incident_edges.push(edge as u32);
    }

    /// Returns the node's name.
    pub fn get_name(&self) -> &str { 
        &self.name
    }

    /// Returns the node's ID.
    pub fn get_id(&self) -> usize {
        self.id as usize
    }

    /// Returns a reference to the node subtype.
    pub fn get_subtype(&self) -> &NodeType {
        &self.subtype
    }

    /// Updates the node subtype.
    ///
    /// # Arguments
    /// * `subtype` - New subtype for the node.
    pub fn set_subtype(&mut self, subtype: NodeType) {
        self.subtype = subtype;
    }

    /// Returns the number of neighbors of the node.
    pub fn neighbors_count(&self) -> usize {
        self.incident_edges.len()
    }

    /// Returns a specific incident edge by neighbor index within the nodes neighbors.
    ///
    /// # Arguments
    /// * `neighbor_id` - Index of neighbor.
    ///
    /// # Returns
    /// Edge identifier.
    pub fn get_incident_edge(&self, neighbor_id: usize) -> usize {
        self.incident_edges[neighbor_id] as usize
    }

    /// Returns all incident edges.
    pub fn get_incident_edges(&self) -> impl Iterator<Item = usize> {
        self.incident_edges.iter().map(|&edge_id| edge_id as usize)
    }

    /// Replaces incident edges with a new list.
    ///
    /// # Arguments
    /// * `new_incident_edges` - Replacement edge list.
    pub fn set_incident_edges(&mut self, new_incident_edges: impl Iterator<Item = usize>) {
        self.incident_edges = new_incident_edges.map(|x| x as u32).collect();
    }

    /// Returns the node category as a string.
    ///
    /// # Returns
    /// `"peptide"`, `"factor"`, `"taxon"`, or `"convolution_tree"`.
    pub fn category(&self) -> &str {
        match self.subtype {
            NodeType::PeptideNode { .. } => "peptide",
            NodeType::FactorNode { .. } => "factor",
            NodeType::TaxonNode { .. } => "taxon",
            NodeType::ConvolutionTreeNode { .. } => "convolution_tree"
        }
    }

    /// Checks if the node is a factor node.
    pub fn is_factor_node(&self) -> bool {
        matches!(self.subtype, NodeType::FactorNode { .. })
    }

    /// Checks if the node is a taxon node.
    pub fn is_taxon_node(&self) -> bool {
        matches!(self.subtype, NodeType::TaxonNode { .. })
    }

    /// Checks if the node is a convolution node.
    pub fn is_convolution_tree_node(&self) -> bool {
        matches!(self.subtype, NodeType::ConvolutionTreeNode { .. })
    }

    /// Updates prior belief for taxon nodes.
    ///
    /// # Arguments
    /// * `prior` - Probability to assign as active.
    pub fn fill_in_prior(&mut self, prior: f64) {
        if matches!(self.subtype, NodeType::TaxonNode { .. }) {
            self.subtype = NodeType::TaxonNode { initial_belief_0: 1.0 - prior, initial_belief_1: prior };
        }
    }

    /// Initializes factor CPD with noisy-OR parameters.
    ///
    /// # Arguments
    /// * `alpha` - Peptide detection probability.
    /// * `beta` - Noise parameter.
    /// * `regularized` - Whether to apply parent-count regularization.
    pub fn fill_in_factor(&mut self, alpha: f64, beta: f64, regularized: bool) {
        if let NodeType::FactorNode { parent_number, .. } = self.subtype {
            let degree: usize = parent_number as usize;

            let mut cpd_array: Vec<[f64; 2]> = Vec::with_capacity(degree + 1);
            let mut cpd_array_regularized = Vec::with_capacity(degree + 1);
            let exponent_array: Vec<usize> = (0..=degree).collect();
            let divide_array: Vec<f64> = std::iter::once(1usize).chain(1..=degree).map(|x| x as f64).collect();
            
            // regularize cpd priors to penalize higher number of parents
            // log domain to avoid underflow
            let mut cpd_sum: f64 = 0.0;
            let mut cpd_regularized_sum: f64 = 0.0;
            for (i, exp) in exponent_array.iter().enumerate() {
                let cpd_0 = (1.0 - alpha).powi(*exp as i32) * (1.0 - beta);
                let cpd_1 = 1.0 - cpd_0;
                cpd_sum += cpd_0 + cpd_1;
                cpd_array.push([cpd_0, cpd_1]);

                let cpd_regularized_0 = (cpd_0.powi(*exp as i32) * (1.0 - beta)) / divide_array[i];
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

    /// Normalizes CPD values in-place.
    ///
    /// # Arguments
    /// * `arr` - CPD array.
    /// * `sum` - Normalization constant.
    /// * `avoid_underflow` - If true, enforce minimum values.
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

    /// Parses key-value data from a GraphML element.
    ///
    /// # Arguments
    /// * `data` - GraphML `<data>` element.
    ///
    /// # Returns
    /// A `(key, value)` pair.
    fn parse_data(data: &Element) -> Result<(String, String), Box<dyn std::error::Error>> {
        let key = data.attr("key").ok_or("key not found while parsing Node data")?.to_string();
        let value = data.text();
    
        Ok((key, value))
    }

    /// Parses a GraphML `<node>` element into a `Node`.
    ///
    /// # Arguments
    /// * `node` - GraphML element to parse.
    /// * `id` - Node ID.
    ///
    /// # Returns
    /// A new `Node` or error if subtype is unknown.
    pub fn parse_node(node: &Element, id: usize) -> Result<Self, Box<dyn std::error::Error>> {
        // Process a node
        let name = node.attr("id").ok_or("id not found while parsing Node")?.to_string();
    
        // Initialize data for this node
        let mut current_node_data = HashMap::new();
        for data in node.children().filter(|d| d.name() == "data") {
            let (data_key, data_val) = Self::parse_data(data)?;
            current_node_data.insert(data_key, data_val);
        }
    
        let subtype: NodeType = match current_node_data.get("d2").map(String::as_str) {
            Some("factor") => {
                let parent_number: usize = current_node_data.get("d3").ok_or("d3 not found while parsing factor Node")?.parse()?;
                let initial_belief: Factor = Factor { array: Vec::new(), array_labels: Vec::new() };
                NodeType::FactorNode { parent_number: parent_number as u32, initial_belief }
            }
            Some("peptide") => {
                let initial_belief_0: f64 = current_node_data.get("d0").ok_or("d0 not found while parsing peptide Node")?.parse()?;
                let initial_belief_1: f64 = current_node_data.get("d1").ok_or("d1 not found while parsing peptide Node")?.parse()?;
                NodeType::PeptideNode { initial_belief_0, initial_belief_1 }
            }
            Some("taxon") => {
                let initial_belief_0: f64 = 0.0;
                let initial_belief_1: f64 = 0.0;
                NodeType::TaxonNode { initial_belief_0, initial_belief_1 }
            }
            _ => {
                return Err("Node data has unknown type".into());
            }
        };

        Ok(Self { id: id as u32, name, incident_edges: Vec::new(), subtype })
    }

    /// Serializes the node into GraphML string format.
    ///
    /// # Returns
    /// A GraphML representation of the node.
    pub fn to_graphml(&self) -> Result<String, Box<dyn std::error::Error>> {

        let mut graphml = String::new();

        writeln!(&mut graphml, r#"<node id="{}">"#, self.name)?;
        writeln!(&mut graphml, r#"<data key="d2">{}</data>"#, self.category())?;

        match &self.subtype {
            NodeType::PeptideNode { initial_belief_0, initial_belief_1 } => {
                writeln!(&mut graphml, r#"<data key="d0">{}</data>"#, initial_belief_0)?;
                writeln!(&mut graphml, r#"<data key="d1">{}</data>"#, initial_belief_1)?;
            },
            NodeType::FactorNode { parent_number, initial_belief: _ } => {
                writeln!(&mut graphml, r#"<data key="d3">{}</data>"#, parent_number)?;
            },
            _ => {}
        }

        writeln!(&mut graphml, r#"</node>"#)?;

        Ok(graphml)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use minidom::Element;

    fn dummy_factor_node(id: usize, parent_number: usize) -> Node {
        Node::new(
            id,
            format!("factor_{}", id),
            NodeType::FactorNode {
                parent_number,
                initial_belief: Factor { array: vec![], array_labels: vec![] },
            },
        )
    }

    #[test]
    fn test_new_and_getters() {
        let node = Node::new(1, "pep".into(), NodeType::PeptideNode { initial_belief_0: 0.3, initial_belief_1: 0.7 });
        assert_eq!(node.get_id(), 1);
        assert_eq!(node.get_name(), "pep");
        assert_eq!(node.category(), "peptide");
        assert!(!node.is_factor_node());
        assert!(!node.is_taxon_node());
    }

    #[test]
    fn test_copy_with_id() {
        let node = Node::new(1, "a".into(), NodeType::TaxonNode { initial_belief_0: 0.0, initial_belief_1: 1.0 });
        let copy = node.copy_with_id(42);
        assert_eq!(copy.get_id(), 42);
        assert_eq!(copy.get_name(), "a");
    }

    #[test]
    fn test_new_convolution_node() {
        let node = Node::new_convolution_node(10, "conv".into(), 3);
        assert_eq!(node.get_id(), 10);
        assert_eq!(node.category(), "convolution_tree");
    }

    #[test]
    fn test_incident_edges() {
        let mut node = Node::new(1, "test".into(), NodeType::TaxonNode { initial_belief_0: 0.1, initial_belief_1: 0.9 });
        node.add_incident_edge(5);
        assert_eq!(node.neighbors_count(), 1);
        assert_eq!(node.get_incident_edge(0), 5);
        assert_eq!(node.get_incident_edges(), &vec![5]);
        node.set_incident_edges(vec![7, 8]);
        assert_eq!(node.get_incident_edges(), &vec![7, 8]);
    }

    #[test]
    fn test_set_and_get_subtype() {
        let mut node = Node::new(1, "x".into(), NodeType::PeptideNode { initial_belief_0: 0.2, initial_belief_1: 0.8 });
        node.set_subtype(NodeType::TaxonNode { initial_belief_0: 0.5, initial_belief_1: 0.5 });
        assert!(matches!(node.get_subtype(), NodeType::TaxonNode { .. }));
    }

    #[test]
    fn test_fill_in_prior() {
        let mut node = Node::new(1, "tax".into(), NodeType::TaxonNode { initial_belief_0: 0.0, initial_belief_1: 0.0 });
        node.fill_in_prior(0.8);
        if let NodeType::TaxonNode { initial_belief_0, initial_belief_1 } = node.get_subtype() {
            assert!((*initial_belief_0 - 0.2).abs() < 1e-9);
            assert!((*initial_belief_1 - 0.8).abs() < 1e-9);
        } else {
            panic!("expected TaxonNode");
        }
    }

    #[test]
    fn test_fill_in_factor() {
        let mut node = dummy_factor_node(2, 2);
        node.fill_in_factor(0.5, 0.1, false);
        if let NodeType::FactorNode { initial_belief, .. } = node.get_subtype() {
            assert!(!initial_belief.array.is_empty());
        } else {
            panic!("expected FactorNode");
        }
    }

    #[test]
    fn test_normalize_cpd() {
        let mut arr = vec![[2.0, 2.0], [1.0, 3.0]];
        Node::normalize_cpd(&mut arr, 8.0, true);
        for row in arr {
            assert!(row[0] >= 1e-30);
            assert!(row[1] >= 1e-30);
            assert!((row[0] + row[1]) <= 1.0 + 1e-9);
        }
    }

    #[test]
    fn test_parse_data() {
        let xml: &str = r#"<data xmlns="ns" key="d0">0.5</data>"#;
        let elem: Element = xml.parse().unwrap();
        let (k, v) = Node::parse_data(&elem);
        assert_eq!(k, "d0");
        assert_eq!(v, "0.5");
    }

    #[test]
    fn test_parse_node_peptide() {
        let xml: &str = r#"<node xmlns="ns" id="n1">
            <data key="d2">peptide</data>
            <data key="d0">0.1</data>
            <data key="d1">0.9</data>
        </node>"#;
        let elem: Element = xml.parse().unwrap();
        let node = Node::parse_node(&elem, 7).unwrap();
        assert_eq!(node.get_id(), 7);
        assert_eq!(node.category(), "peptide");
    }

    #[test]
    fn test_to_graphml() {
        let node = Node::new(1, "pepx".into(), NodeType::PeptideNode { initial_belief_0: 0.4, initial_belief_1: 0.6 });
        let xml = node.to_graphml();
        assert!(xml.contains("<node id=\"pepx\">"));
        assert!(xml.contains("<data key=\"d2\">peptide</data>"));
    }
}
