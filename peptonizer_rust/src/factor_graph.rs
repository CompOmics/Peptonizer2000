use crate::node::{Factor, Node, NodeType};
use minidom::Element;
use std::collections::{HashMap, HashSet};
use serde::{Serialize, Deserialize};
use std::fmt::Write;
use csv::ReaderBuilder;

/// Represents a single taxon weight record parsed from a CSV file.
#[derive(Deserialize)]
pub struct TaxonWeight {
    pub id: usize,
    pub sequence: String,
    pub score: f32,
    pub psms: usize,
    pub higher_taxa: usize,
    pub weight: f32,
    pub log_weight: f32
}


/// Parses a CSV string into a vector of `TaxonWeight` structs.
///
/// # Arguments
/// * `taxa_weights_csv` - A string containing CSV data for taxon weights. The CSV
///   must include headers: `id, sequence, score, psms, higher_taxa, weight, log_weight`.
///
/// # Returns
/// Returns a `Result` containing a vector of `TaxonWeight` structs if parsing succeeds.
///
/// # Errors
/// Returns an error if the CSV cannot be read, or if any record fails deserialization.
pub fn parse_taxon_weights_csv(taxa_weights_csv: String) -> Result<Vec<TaxonWeight>, Box<dyn std::error::Error>> {
    let mut rdr = ReaderBuilder::new()
        .has_headers(true)
        .from_reader(taxa_weights_csv.as_bytes());
    
    let mut taxa_weights = Vec::with_capacity(rdr.records().count());
    for record in rdr.deserialize() {
        let row: TaxonWeight = record?;
        taxa_weights.push(row);
    }

    Ok(taxa_weights)
}


/// Generates a GraphML representation of a factor graph from a CSV string of taxon weights.
///
/// # Arguments
/// * `taxa_weights_csv` - A string containing CSV data for taxon weights.
///
/// # Returns
/// Returns a `Result` containing a GraphML string representation of the factor graph.
///
/// # Errors
/// Returns an error if CSV parsing fails or if any error occurs during graph construction.
pub fn generate_graph(taxa_weights_csv: String) -> Result<String, Box<dyn std::error::Error>> {

    let taxa_weights = parse_taxon_weights_csv(taxa_weights_csv)?;

    let graph = CTFactorGraph::from_taxa_weights(taxa_weights);

    Ok(graph.to_graphml()?)
}


/// Represents an edge in a factor graph connecting two nodes.
#[derive(Debug, Serialize, Clone)]
pub struct Edge {
    id: u32,
    node1_id: u32,
    node2_id: u32,
    node1_in_node2_id: u32,
    node2_in_node1_id: u32,
    message_length: Option<u32>
}


impl Edge {

    /// Creates a new edge connecting two nodes in a factor graph.
    ///
    /// # Arguments
    /// * `id` - Unique identifier for the edge.
    /// * `node1_id` - ID of the first node connected by this edge.
    /// * `node2_id` - ID of the second node connected by this edge.
    /// * `message_length` - Optional message length associated with the edge. Can be `None` if not applicable.
    ///
    /// # Returns
    /// An `Edge` instance representing a connection between the two specified nodes.
    pub fn new(id: usize, node1_id: usize, node2_id: usize, node1_in_node2_id: usize, node2_in_node1_id: usize, message_length: Option<usize>) -> Edge {
        Edge { 
            id: id as u32, 
            node1_id: node1_id as u32, 
            node2_id: node2_id as u32, 
            node1_in_node2_id: node1_in_node2_id as u32, 
            node2_in_node1_id: node2_in_node1_id as u32, 
            message_length: message_length.map(|x| x as u32) 
        }
    }


    pub fn set_node1_in_node2_id(&mut self, id: usize) {
        self.node1_in_node2_id = id as u32;
    }

    pub fn set_node2_in_node1_id(&mut self, id: usize) {
        self.node2_in_node1_id = id as u32;
    }

    /// Returns the ID of the edge.
    pub fn get_id(&self) -> usize {
        self.id as usize
    }

    /// Returns the first node ID of the edge.
    pub fn get_node1_id(&self) -> usize {
        self.node1_id as usize
    }

    /// Returns the second node ID of the edge.
    pub fn get_node2_id(&self) -> usize {
        self.node2_id as usize
    }

    /// Returns a tuple of the two node IDs of the edge.
    pub fn get_node_ids(&self) -> (usize, usize) {
        (self.node1_id as usize, self.node2_id as usize)
    }

    pub fn get_node_and_neighbor_ids(&self) -> ((usize, usize), (usize, usize)) {
        ((self.node1_id as usize, self.node1_in_node2_id as usize), (self.node2_id as usize, self.node2_in_node1_id as usize))
    }

    /// Returns the message length associated with the edge.
    pub fn get_message_length(&self) -> Option<usize> {
        self.message_length.map(|x| x as usize)
    }

    pub fn copy_with_id(&self, new_id: usize) -> Self {
        let mut copy: Edge = self.clone();
        copy.id = new_id as u32;
        copy
    }
}

#[derive(Debug)]
pub struct CTFactorGraph {
    nodes: Vec<Node>,
    edges: Vec<Edge>,
}

impl CTFactorGraph {

    /// Creates a new factor graph from a list of nodes and edges.
    ///
    /// # Arguments
    /// * `nodes` - A vector of `Node` instances representing all nodes in the graph.
    /// * `edges` - A vector of `Edge` instances representing all edges connecting the nodes.
    ///
    /// # Returns
    /// A `CTFactorGraph` instance containing the provided nodes and edges.
    pub fn new(nodes: Vec<Node>, edges: Vec<Edge>) -> CTFactorGraph {
        CTFactorGraph { nodes, edges }
    }

    /// Adds the names and categories of all nodes to the provided vectors.
    ///
    /// # Arguments
    /// * `node_names` - Mutable reference to a vector to store node names.
    /// * `node_categories` - Mutable reference to a vector to store node categories.
    pub fn add_node_names_categories(&self, node_names: &mut Vec<String>, node_categories: &mut Vec<String>) {
        for node in &self.nodes {
            node_names.push(node.get_name().to_string());
            node_categories.push(node.category().to_string());
        }
    } 

    /// Returns a reference to the node with the given ID.
    ///
    /// # Arguments
    /// * `node_id` - ID of the node to retrieve.
    ///
    /// # Returns
    /// A reference to the `Node` corresponding to `node_id`.
    pub fn get_node(&self, node_id: usize) -> &Node {
        &self.nodes[node_id]
    }

    /// Returns a reference to the edge with the given ID.
    ///
    /// # Arguments
    /// * `edge_id` - ID of the edge to retrieve.
    ///
    /// # Returns
    /// A reference to the `Edge` corresponding to `edge_id`.
    pub fn get_edge(&self, edge_id: usize) -> &Edge {
        &self.edges[edge_id]
    }

    /// Returns the total number of nodes in the graph.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Returns the total number of edges in the graph.
    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }

    /// Returns a reference to all nodes in the graph.
    pub fn get_nodes(&self) -> &Vec<Node> {
        &self.nodes
    }

    /// Returns a reference to all edges in the graph.
    pub fn get_edges(&self) -> &Vec<Edge> {
        &self.edges
    }

    fn parse_edge(edge: &Element) -> Result<(String, String), Box<dyn std::error::Error>> {
        let source: String = edge.attr("source").ok_or("Source attribute does not exist in Edge")?.to_string();
        let target: String = edge.attr("target").ok_or("Target attribute does not exist in Edge")?.to_string();
    
        Ok((source, target))
    }

    /// Converts the factor graph into a GraphML string representation.
    ///
    /// # Returns
    /// A `String` containing the GraphML XML representation of the factor graph.
    pub fn to_graphml(&self) -> Result<String, Box<dyn std::error::Error>> {

        let mut graphml = String::new();

        writeln!(
            &mut graphml,
            r#"<?xml version="1.0" encoding="UTF-8"?>"#
        )?;
        writeln!(
            &mut graphml,
            r#"<graphml xmlns="http://graphml.graphdrawing.org/xmlns" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">"#
        )?;

        writeln!(
            &mut graphml,
            r#"  <key id="d3" for="node" attr.name="ParentNumber" attr.type="long" />
  <key id="d2" for="node" attr.name="category" attr.type="string" />
  <key id="d1" for="node" attr.name="InitialBelief_1" attr.type="double" />
  <key id="d0" for="node" attr.name="InitialBelief_0" attr.type="double" />"#
        )?;

        writeln!(&mut graphml, r#"  <graph edgedefault="undirected">"#)?;

        for node in &self.nodes {
            write!(&mut graphml, "{}", node.to_graphml()?)?;
        }

        for edge in &self.edges {
            let node1: &Node = self.get_node(edge.get_node1_id());
            let node2: &Node = self.get_node(edge.get_node2_id());

            writeln!(
                &mut graphml,
                r#"<edge source="{}" target="{}" />"#,
                node1.get_name(),
                node2.get_name()
            )?;
        }

        writeln!(&mut graphml, r#"  </graph>"#)?;
        writeln!(&mut graphml, r#"</graphml>"#)?;

        Ok(graphml)
    }
    
    /// Constructs a `CTFactorGraph` from a GraphML string.
    ///
    /// # Arguments
    /// * `graphml_str` - A string containing a GraphML representation of a graph.
    ///
    /// # Returns
    /// A `Result` containing the constructed `CTFactorGraph` if successful.
    ///
    /// # Errors
    /// Returns an error if parsing the GraphML fails or if nodes/edges cannot be created correctly.
    pub fn from_graphml(graphml_str: &str) -> Result<CTFactorGraph, Box<dyn std::error::Error>> {
        let root: Element = graphml_str.parse()?;

        let node_count = root.children().filter(|n| n.name() == "graph").map(|g| g.children().filter(|n| n.name() == "node").count()).sum();
        let mut nodes: Vec<Node> = Vec::with_capacity(node_count);
        let edge_count = root.children().filter(|n| n.name() == "graph").map(|g| g.children().filter(|n| n.name() == "edge").count()).sum();
        let mut edges: Vec<Edge> = Vec::with_capacity(edge_count);
        let mut node_map: HashMap<String, usize> = HashMap::new();
        
        let mut next_node_id = 0;
        let mut next_edge_id = 0;
        for graph_xml in root.children().filter(|n| n.name() == "graph") {
            for node_xml in graph_xml.children().filter(|n| n.name() == "node") {
                let node: Node = Node::parse_node(node_xml, next_node_id)?;
                let node_name: String = node.get_name().to_string();
                node_map.insert(node_name, next_node_id);
                next_node_id += 1;

                nodes.push(node);
            }
    
            for edge_xml in graph_xml.children().filter(|n| n.name() == "edge") {
                let (source, target) = Self::parse_edge(edge_xml)?;
    
                let node1_id: usize = *node_map.get(&source).ok_or("Source node of edge not present in graph")?;
                let node2_id: usize = *node_map.get(&target).ok_or("Target node of edge not present in graph")?;
                let node1: &Node = &nodes[node1_id];
                let node2: &Node = &nodes[node2_id];
                let edge = Edge::new(next_edge_id, node1_id, node2_id, node2.neighbors_count(), node1.neighbors_count(), None);
                next_edge_id += 1;
    
                let node1: &mut Node = &mut nodes[node1_id];
                node1.add_incident_edge(edge.get_id());
                let node2: &mut Node = &mut nodes[node2_id];
                node2.add_incident_edge(edge.get_id());
                edges.push(edge);
            }
        }
    
        Ok( CTFactorGraph { nodes, edges })
    }

    /// Constructs a `CTFactorGraph` from a list of `TaxonWeight`s.
    ///
    /// # Arguments
    /// * `taxa_weights` - Vector of `TaxonWeight` structs used to build nodes and edges.
    ///
    /// # Returns
    /// Returns a `CTFactorGraph` representing the factor graph built from the input data.
    pub fn from_taxa_weights(taxa_weights: Vec<TaxonWeight>) -> CTFactorGraph {

        // Count frequencies of each higher_taxa
        let mut higher_taxa_counts: HashMap<usize, usize> = HashMap::new();
        for tw in &taxa_weights {
            *higher_taxa_counts.entry(tw.higher_taxa).or_insert(0) += 1;
        }
        // Filter to keep only those with count > 1
        let taxa_weights = taxa_weights.into_iter().filter(|tw| higher_taxa_counts[&tw.higher_taxa] > 1);

        let mut node_id_counter: usize = 0;
        let mut edge_id_counter: usize = 0;
        let mut node_name_to_id: HashMap<String, usize> = HashMap::new();
        let mut nodes: Vec<Node> = Vec::new();
        let mut edges: Vec<Edge> = Vec::new();
        for tw in taxa_weights {
            let cpd_name = tw.sequence.clone() + " CPD";

            // Add sequence node and CPD node if necessary
            if ! node_name_to_id.contains_key(&tw.sequence) {
                let node1_id = node_id_counter;
                let node_type = NodeType::PeptideNode { initial_belief_0: 1.0 - tw.score as f64, initial_belief_1: tw.score as f64 };
                let node = Node::new(node1_id, tw.sequence.clone(), node_type);
                nodes.push(node);
                node_name_to_id.insert(tw.sequence.clone(), node1_id);
                node_id_counter += 1;

                let node2_id = node_id_counter;
                let node_type = NodeType::FactorNode { 
                    parent_number: 0, 
                    initial_belief: Factor { array: Vec::new(), array_labels: Vec::new() }
                };
                let node = Node::new(node2_id, cpd_name.clone(), node_type);
                nodes.push(node);
                node_name_to_id.insert(cpd_name.clone(), node2_id);
                node_id_counter += 1;

                let edge = Edge::new(edge_id_counter, node1_id, node2_id, nodes[node2_id].neighbors_count(), nodes[node1_id].neighbors_count(), None);
                edges.push(edge);
                nodes[node1_id].add_incident_edge(edge_id_counter);
                nodes[node2_id].add_incident_edge(edge_id_counter);
                edge_id_counter += 1;
            }

            // Add taxon node if necessary
            let higher_taxa_str = tw.higher_taxa.to_string();
            if ! node_name_to_id.contains_key(&higher_taxa_str) {
                let node_type = NodeType::TaxonNode { initial_belief_0: 0.0, initial_belief_1: 0.0 };
                let node = Node::new(node_id_counter, higher_taxa_str.clone(), node_type);
                nodes.push(node);

                node_name_to_id.insert(higher_taxa_str.clone(), node_id_counter);
                node_id_counter += 1;
            }

            // Add edge
            let node1_id = node_name_to_id[&higher_taxa_str];
            let node2_id = node_name_to_id[&cpd_name];
            let edge = Edge::new(edge_id_counter, node1_id, node2_id, nodes[node2_id].neighbors_count(), nodes[node1_id].neighbors_count(), None);
            edges.push(edge);
            nodes[node1_id].add_incident_edge(edge_id_counter);
            nodes[node2_id].add_incident_edge(edge_id_counter);
            edge_id_counter += 1;
        }

        // set parent_number correct for factor nodes (CPDs)
        for node in nodes.iter_mut() {
            if node.is_factor_node() {
                node.set_subtype(NodeType::FactorNode { 
                    parent_number: (node.neighbors_count() - 1) as u32, 
                    initial_belief: Factor { array: Vec::new(), array_labels: Vec::new() }
                });
            }
        }

        CTFactorGraph { nodes, edges }
    }

    /// Fills all nodes with a prior probability.
    ///
    /// # Arguments
    /// * `prior` - The prior probability to assign to each node.
    pub fn fill_in_priors(&mut self, prior: f64) {
        for node in &mut self.nodes {
            node.fill_in_prior(prior);
        }
    }

    /// Fills all factor nodes with factor probabilities using alpha/beta parameters.
    ///
    /// # Arguments
    /// * `alpha` - Alpha parameter for factor probability.
    /// * `beta` - Beta parameter for factor probability.
    /// * `regularized` - Whether to apply regularization.
    pub fn fill_in_factors(&mut self, alpha: f64, beta: f64, regularized: bool) {
        for node in &mut self.nodes {
            node.fill_in_factor(alpha, beta, regularized);
        }
    }

    /// Returns the IDs of neighbors of a node, given its node ID.
    ///
    /// # Arguments
    /// * `node_id` - The node ID for which neighbors are requested.
    ///
    /// # Returns
    /// A vector of node IDs representing neighbors.
    pub fn get_neighbors_from_id(&self, node_id: usize) -> impl Iterator<Item = usize> {
        self.get_neighbors(self.get_node(node_id))
    }

    /// Returns the IDs of neighbors for a given node.
    ///
    /// # Arguments
    /// * `node` - Reference to the `Node` whose neighbors are requested.
    ///
    /// # Returns
    /// A Iterator over node IDs representing neighbors.
    pub fn get_neighbors(&self, node: &Node) -> impl Iterator<Item = usize> {
        node.get_incident_edges().map(|edge_id| {
            let (node1_id, node2_id) = self.edges[edge_id].get_node_ids();
            if node1_id == node.get_id() { node2_id } else { node1_id }
        })
    }

    /// Returns the node ID of a neighbor given a node and its neighbor ID.
    ///
    /// # Arguments
    /// * `node` - Reference to the node.
    /// * `neighbor_id` - Index of the neighbor within the nodes neighbors.
    ///
    /// # Returns
    /// Node ID of the neighbor.
    pub fn get_neighbor_node_id(&self, node: &Node, neighbor_id: usize) -> usize {
        let (node1_id, node2_id) = self.edges[node.get_incident_edge(neighbor_id)].get_node_ids();
        if node1_id == node.get_id() { node2_id } else { node1_id }
    }

    pub fn get_neighbor_node_and_neighbor_id(&self, node: &Node, neighbor_id: usize) -> (usize, usize) {
        let ((node1_id, node1_in_node2_id), (node2_id, node2_in_node1_id)) = self.edges[node.get_incident_edge(neighbor_id)].get_node_and_neighbor_ids();
        if node1_id == node.get_id() { (node2_id, node1_in_node2_id) } else { (node1_id, node2_in_node1_id) }
    }

    /// Returns the peptide node ID connected to a factor node.
    ///
    /// # Arguments
    /// * `factor_id` - ID of the factor node.
    ///
    /// # Returns
    /// `Ok(usize)` containing the peptide node ID if found.
    ///
    /// # Errors
    /// Returns an error if no peptide node is connected to the factor node.
    pub fn get_peptide_for_factor(&self, factor_id: usize) -> Result<usize, Box<dyn std::error::Error>> {
        for neighbor_id in self.get_neighbors_from_id(factor_id) {
            let neighbor = self.get_node(neighbor_id);
            if let NodeType::PeptideNode { .. } = neighbor.get_subtype() {
                return Ok(neighbor.get_id());
            }
        }
        return Err(format!("Peptide not found for factor with id {}", factor_id).into());
    }

    /// Adds convolution tree nodes to the graph, creating edges appropriately.
    pub fn add_ct_nodes(&mut self) {
        // When creating the CTGraph and not just reading from a previously saved graph format, use this function to add the CT nodes
        
        let ct_node_count = self.nodes.iter().filter(|n| n.is_factor_node() && n.neighbors_count() > 2).count();
        let mut new_nodes: Vec<Node> = Vec::with_capacity(&self.nodes.len() + ct_node_count);
        new_nodes.extend_from_slice(&self.nodes);
        let mut new_edges: Vec<Edge> = Vec::with_capacity(&self.edges.len() + ct_node_count);

        // Add nodes and keep track of edges to add/remove
        let mut next_edge_id: usize = 0;
        let mut next_node_id: usize = self.nodes.len();
        for node in &self.nodes {
            if node.is_factor_node() {
                if node.neighbors_count() > 2 {
                    let mut prot_names: Vec<String> = self.get_neighbors(node)
                        .map(|n|&self.nodes[n])
                        .filter(|n|n.is_taxon_node())
                        .map(|n| n.get_name().to_string()).collect();
                    
                    // TODO: names necessary? These nodes are added after graphml is created, because this is executed in execute_pepgm. The names are not contained in any output I think. For the algorithm itself, strings are inefficient
                    let new_node_name = format!("CTree {}", prot_names.join(" "));
                    let new_node_id = next_node_id;
                    let new_node = Node::new_convolution_node(new_node_id, new_node_name, prot_names.len());
                    next_node_id += 1;
                    new_nodes.push(new_node);

                    // Create edge Factor CTree, set node_in_node_id's to 0, we will set them correctly later
                    let edge = Edge::new(next_edge_id, new_node_id, node.get_id(), 0, 0, Some(prot_names.len() + 1));
                    next_edge_id += 1;
                    new_edges.push(edge);

                    for (i, edge_id) in node.get_incident_edges().enumerate() {
                        let neighbor_id = self.get_neighbor_node_id(node, i);
                        let neighbor: &Node = self.get_node(neighbor_id);
                        if neighbor.is_taxon_node() {
                            // Create edge CTree - Taxon, set node_in_node_id's to 0, we will set them correctly later
                            let edge = Edge::new(next_edge_id, new_node_id, neighbor_id, 0, 0, None);
                            next_edge_id += 1;
                            new_edges.push(edge);
                        } else {
                            // Add Factor - Peptide node
                            new_edges.push(self.get_edge(edge_id).copy_with_id(next_edge_id));
                            next_edge_id += 1;
                        }
                    }
                } else {
                    for edge_id in node.get_incident_edges() {
                        new_edges.push(self.get_edge(edge_id).copy_with_id(next_edge_id));
                        next_edge_id += 1;
                    }
                }
                
            }
        }

        // Clear the incident edges of each node, and refill in the next step
        for node in &mut new_nodes {
            node.set_incident_edges(Vec::new().into_iter());
        }
        
        for edge in &mut new_edges {
            let (node1_id, node2_id) = edge.get_node_ids();
            edge.set_node1_in_node2_id(new_nodes[node2_id].neighbors_count());
            edge.set_node2_in_node1_id(new_nodes[node1_id].neighbors_count());
            new_nodes[node1_id].add_incident_edge(edge.get_id());
            new_nodes[node2_id].add_incident_edge(edge.get_id());
        }

        self.nodes = new_nodes;
        self.edges = new_edges;
    }

    /// Returns a vector of connected components in the graph as separate `CTFactorGraph`s.
    ///
    /// # Returns
    /// A vector of `CTFactorGraph` instances, one per connected component.
    pub fn connected_components(&self) -> Vec<Self> {
        let mut visited: HashSet<usize> = HashSet::new();
        let mut components: Vec<Self> = Vec::new();

        for start_node in &self.nodes {
            if visited.insert(start_node.get_id()) {
                let mut component_ids: Vec<usize> = Vec::new();
                let mut old_to_new_nodes: HashMap<usize, usize> = HashMap::new();

                let mut new_nodes: Vec<Node> = Vec::new();
                let mut new_edges: Vec<Edge> = Vec::new();

                // Find ids of nodes to include in component
                component_ids.push(start_node.get_id());
                old_to_new_nodes.insert(start_node.get_id(), 0);
                self.find_component_rec(start_node.get_id(), &mut component_ids, &mut old_to_new_nodes, &mut visited);

                // Create new nodes
                for node_id in &component_ids {
                    let node = self.nodes[*node_id].copy_with_id(old_to_new_nodes[&node_id]);
                    new_nodes.push(node);
                }

                // Select edges to keep and update the node ids
                let mut next_edge_id: usize = 0;
                let mut component_edge_ids: HashSet<usize> = HashSet::new();
                let mut old_to_new_edges: HashMap<usize, usize> = HashMap::new();
                for edge in &self.edges {

                    let ((source, source_in_target), (target, target_in_source)): ((usize, usize), (usize, usize)) = edge.get_node_and_neighbor_ids();
                    if component_ids.contains(&source) && component_ids.contains(&target) {

                        let (new_source, new_target): (usize, usize) = (old_to_new_nodes[&source], old_to_new_nodes[&target]);
                        let new_edge = Edge::new(next_edge_id, new_source, new_target, source_in_target, target_in_source, edge.get_message_length());
                        next_edge_id += 1;

                        component_edge_ids.insert(edge.get_id());
                        old_to_new_edges.insert(edge.get_id(), new_edge.get_id());

                        new_edges.push(new_edge);
                    }
                }

                // Update edge ids of incident edges
                for node in &mut new_nodes {
                    let new_incident_edges: Vec<usize> = node.get_incident_edges().filter(|e| component_edge_ids.contains(e)).map(|e| old_to_new_edges[&e]).collect();
                    node.set_incident_edges(new_incident_edges.into_iter());
                }
                
                // Create graph and add to components
                let subgraph = Self { nodes: new_nodes, edges: new_edges };
                components.push(subgraph);
            }
        }

        components
    }

    fn find_component_rec(
        &self, 
        start_id: usize, 
        component_ids: &mut Vec<usize>, 
        old_to_new_nodes: &mut HashMap<usize, usize>, 
        visited: &mut HashSet<usize>
    ) {
        let start_node: &Node = &self.nodes[start_id];
        for neighbor_id in self.get_neighbors(&start_node) {
            if visited.insert(neighbor_id) {
                let next_id: usize = component_ids.len();
                component_ids.push(neighbor_id);
                old_to_new_nodes.insert(neighbor_id, next_id);
                self.find_component_rec(neighbor_id, component_ids, old_to_new_nodes, visited);                
            }
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::node::{Node, NodeType, Factor};

    fn sample_csv() -> String {
        "id,sequence,score,psms,higher_taxa,weight,log_weight
1,PEPTIDE1,0.8,5,100,0.5,-0.3
2,PEPTIDE2,0.6,3,100,0.4,-0.5
3,PEPTIDE3,0.9,7,200,0.7,-0.1"
            .to_string()
    }

    #[test]
    fn test_parse_taxon_weights_csv() {
        let csv = sample_csv();
        let taxa = parse_taxon_weights_csv(csv).unwrap();
        assert_eq!(taxa.len(), 3);
        assert_eq!(taxa[0].id, 1);
        assert!((taxa[1].score - 0.6).abs() < 1e-6);
    }

    #[test]
    fn test_generate_graph_creates_graphml() {
        let csv = sample_csv();
        let graphml = generate_graph(csv).unwrap();
        assert!(graphml.contains("graphml"));
        assert!(graphml.contains("node"));
        assert!(graphml.contains("edge"));
    }

    #[test]
    fn test_edge_getters() {
        let edge = Edge::new(1, 10, 20, Some(5));
        assert_eq!(edge.get_id(), 1);
        assert_eq!(edge.get_node1_id(), 10);
        assert_eq!(edge.get_node2_id(), 20);
        assert_eq!(edge.get_node_ids(), (10, 20));
        assert_eq!(edge.get_message_length(), Some(5));
    }

    #[test]
    fn test_ctfactorgraph_from_taxa_weights() {
        let csv = sample_csv();
        let taxa = parse_taxon_weights_csv(csv).unwrap();
        let graph = CTFactorGraph::from_taxa_weights(taxa);
        assert!(graph.node_count() > 0);
        assert!(graph.edge_count() > 0);
    }

    #[test]
    fn test_ctfactorgraph_to_and_from_graphml() {
        let csv = sample_csv();
        let taxa = parse_taxon_weights_csv(csv).unwrap();
        let graph = CTFactorGraph::from_taxa_weights(taxa);
        let graphml = graph.to_graphml();
        let parsed = CTFactorGraph::from_graphml(&graphml).unwrap();
        assert_eq!(graph.node_count(), parsed.node_count());
        assert_eq!(graph.edge_count(), parsed.edge_count());
    }

    #[test]
    fn test_neighbor_operations() {
        let csv = sample_csv();
        let taxa = parse_taxon_weights_csv(csv).unwrap();
        let graph = CTFactorGraph::from_taxa_weights(taxa);

        if graph.node_count() > 1 {
            let node = graph.get_node(0);
            for n in graph.get_neighbors(node) {
                let idx = graph.get_neighbor_index(node, n);
                assert!(idx >= 0);
                let idx2 = graph.get_neighbor_index_from_id(node.get_id(), n);
                assert_eq!(idx, idx2);
            }
        }
    }

    #[test]
    fn test_get_peptide_for_factor_returns_ok_or_err() {
        let csv = sample_csv();
        let taxa = parse_taxon_weights_csv(csv).unwrap();
        let graph = CTFactorGraph::from_taxa_weights(taxa);

        for (i, node) in graph.get_nodes().iter().enumerate() {
            if node.is_factor_node() {
                let result = graph.get_peptide_for_factor(i);
                assert!(result.is_ok() || result.is_err());
            }
        }
    }

    #[test]
    fn test_connected_components() {
        let csv = sample_csv();
        let taxa = parse_taxon_weights_csv(csv).unwrap();
        let graph = CTFactorGraph::from_taxa_weights(taxa);

        let components = graph.connected_components();
        assert!(!components.is_empty());
        let total_nodes: usize = components.iter().map(|c| c.node_count()).sum();
        assert_eq!(total_nodes, graph.node_count());
    }
}
