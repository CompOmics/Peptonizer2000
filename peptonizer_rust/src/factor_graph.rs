use crate::node::{Factor, Node, NodeType};
use minidom::Element;
use std::collections::{HashMap, HashSet};
use serde::{Serialize, Deserialize};
use std::fmt::Write;
use csv::ReaderBuilder;

#[derive(Deserialize)]
pub struct TaxonWeight {
    pub id: i32,
    pub sequence: String,
    pub score: f32,
    pub psms: i32,
    pub higher_taxa: i32,
    pub weight: f32,
    pub log_weight: f32
}


pub fn parse_taxon_weights_csv(taxa_weights_csv: String) -> Result<Vec<TaxonWeight>, Box<dyn std::error::Error>> {
    let mut rdr = ReaderBuilder::new()
        .has_headers(true)
        .from_reader(taxa_weights_csv.as_bytes());
    
    let mut taxa_weights = Vec::new();
    for record in rdr.deserialize() {
        let row: TaxonWeight = record?;
        taxa_weights.push(row);
    }

    Ok(taxa_weights)
}


pub fn generate_graph(taxa_weights_csv: String) -> Result<String, Box<dyn std::error::Error>> {

    let taxa_weights = parse_taxon_weights_csv(taxa_weights_csv)?;

    let graph = CTFactorGraph::from_taxa_weights(taxa_weights);

    Ok(graph.to_graphml())
}

#[derive(Debug, Serialize, Clone)]
pub struct Edge {
    id: i32,
    node1_id: i32,
    node2_id: i32,
    message_length: Option<i32>
}

impl Edge {

    pub fn get_id(&self) -> i32 {
        self.id
    }

    pub fn get_node1_id(&self) -> i32 {
        self.node1_id
    }

    pub fn get_node2_id(&self) -> i32 {
        self.node2_id
    }

    pub fn get_node_ids(&self) -> (i32, i32) {
        (self.node1_id, self.node2_id)
    }

    pub fn get_message_length(&self) -> Option<i32> {
        self.message_length
    }
}

#[derive(Debug)]
pub struct CTFactorGraph {
    nodes: Vec<Node>,
    edges: Vec<Edge>,
}

impl CTFactorGraph {

    pub fn add_node_names_categories(&self, node_names: &mut Vec<String>, node_categories: &mut Vec<String>) {
        for node in &self.nodes {
            node_names.push(node.get_name().to_string());
            node_categories.push(node.category().to_string());
        }
    } 

    pub fn get_node(&self, node_id: i32) -> &Node {
        &self.nodes[node_id as usize]
    }

    pub fn get_edge(&self, edge_id: i32) -> &Edge {
        &self.edges[edge_id as usize]
    }

    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }

    pub fn get_nodes(&self) -> &Vec<Node> {
        &self.nodes
    }

    pub fn get_edges(&self) -> &Vec<Edge> {
        &self.edges
    }

    fn parse_edge(edge: &Element) -> Result<(String, String), String> {
        let source: String = edge.attr("source").unwrap().to_string();
        let target: String = edge.attr("target").unwrap().to_string();
    
        Ok((source, target))
    }


    pub fn to_graphml(&self) -> String {

        let mut graphml = String::new();

        writeln!(
            &mut graphml,
            r#"<?xml version="1.0" encoding="UTF-8"?>"#
        ).unwrap();
        writeln!(
            &mut graphml,
            r#"<graphml xmlns="http://graphml.graphdrawing.org/xmlns" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">"#
        ).unwrap();

        writeln!(
            &mut graphml,
            r#"  <key id="d3" for="node" attr.name="ParentNumber" attr.type="long" />
  <key id="d2" for="node" attr.name="category" attr.type="string" />
  <key id="d1" for="node" attr.name="InitialBelief_1" attr.type="double" />
  <key id="d0" for="node" attr.name="InitialBelief_0" attr.type="double" />"#
        ).unwrap();

        writeln!(&mut graphml, r#"  <graph edgedefault="undirected">"#).unwrap();

        for node in &self.nodes {
            write!(&mut graphml, "{}", node.to_graphml()).unwrap();
        }

        for edge in &self.edges {
            let node1: &Node = self.get_node(edge.get_node1_id());
            let node2: &Node = self.get_node(edge.get_node2_id());

            writeln!(
                &mut graphml,
                r#"<edge source="{}" target="{}" />"#,
                node1.get_name(),
                node2.get_name()
            ).unwrap();
        }

        writeln!(&mut graphml, r#"  </graph>"#).unwrap();
        writeln!(&mut graphml, r#"</graphml>"#).unwrap();

        graphml
    }
    
    // Method to parse a GraphML string into the graph
    pub fn from_graphml(graphml_str: &str) -> Result<CTFactorGraph, String> {
        let mut nodes: Vec<Node> = Vec::new();
        let mut edges: Vec<Edge> = Vec::new();
        let mut node_map: HashMap<String, i32> = HashMap::new();
    
        let root: Element = graphml_str.parse().unwrap();
        
        let mut next_node_id = 0;
        let mut next_edge_id = 0;
        for graph_xml in root.children().filter(|n| n.name() == "graph") {
            for node_xml in graph_xml.children().filter(|n| n.name() == "node") {
                let node: Node = Node::parse_node(node_xml, next_node_id).unwrap();
                let node_name: String = node.get_name().to_string();
                node_map.insert(node_name, next_node_id);
                next_node_id += 1;

                nodes.push(node);
            }
    
            for edge_xml in graph_xml.children().filter(|n| n.name() == "edge") {
                let (source, target) = Self::parse_edge(edge_xml).unwrap();
    
                let node1_id: i32 = *node_map.get(&source).unwrap();
                let node2_id: i32 = *node_map.get(&target).unwrap();
                let edge: Edge = Edge { id: next_edge_id, node1_id, node2_id, message_length: None };
                next_edge_id += 1;
    
                let node1: &mut Node = &mut nodes[node1_id as usize];
                node1.add_incident_edge(edge.get_id());
                let node2: &mut Node = &mut nodes[node2_id as usize];
                node2.add_incident_edge(edge.get_id());
                edges.push(edge);
            }
        }
    
        Ok( CTFactorGraph { nodes, edges })
    }

    pub fn from_taxa_weights(taxa_weights: Vec<TaxonWeight>) -> CTFactorGraph {

        // Count frequencies of each higher_taxa
        let mut higher_taxa_counts: HashMap<i32, usize> = HashMap::new();
        for tw in &taxa_weights {
            *higher_taxa_counts.entry(tw.higher_taxa).or_insert(0) += 1;
        }
        // Filter to keep only those with count > 1
        let taxa_weights = taxa_weights.into_iter().filter(|tw| higher_taxa_counts[&tw.higher_taxa] > 1);

        let mut node_id_counter: i32 = 0;
        let mut edge_id_counter: i32 = 0;
        let mut node_name_to_id: HashMap<String, i32> = HashMap::new();
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

                let edge = Edge { id: edge_id_counter, node1_id, node2_id, message_length: None };
                edges.push(edge);
                nodes[node1_id as usize].add_incident_edge(edge_id_counter);
                nodes[node2_id as usize].add_incident_edge(edge_id_counter);
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
            let edge = Edge { id: edge_id_counter, node1_id, node2_id, message_length: None };
            edges.push(edge);
            nodes[node1_id as usize].add_incident_edge(edge_id_counter);
            nodes[node2_id as usize].add_incident_edge(edge_id_counter);
            edge_id_counter += 1;
        }

        // set parent_number correct for factor nodes (CPDs)
        for node in nodes.iter_mut() {
            if node.is_factor_node() {
                node.set_subtype(NodeType::FactorNode { 
                    parent_number: node.neighbors_count() as i32 - 1, 
                    initial_belief: Factor { array: Vec::new(), array_labels: Vec::new() }
                });
            }
        }

        CTFactorGraph { nodes, edges }
    }

    pub fn fill_in_priors(&mut self, prior: f64) {
        for node in &mut self.nodes {
            node.fill_in_prior(prior);
        }
    }

    pub fn fill_in_factors(&mut self, alpha: f64, beta: f64, regularized: bool) {
        for node in &mut self.nodes {
            node.fill_in_factor(alpha, beta, regularized);
        }
    }

    pub fn get_neighbors_from_id(&self, node_id: i32) -> Vec<i32> {
        self.get_neighbors(self.get_node(node_id))
    }

    pub fn get_neighbors(&self, node: &Node) -> Vec<i32> {
        let mut neighbors = Vec::with_capacity(node.neighbors_count() as usize);
        for edge_id in node.get_incident_edges() {
            let (node1_id, node2_id) = self.edges[*edge_id as usize].get_node_ids();
            let neighbor: i32 = if node1_id == node.get_id() { node2_id } else { node1_id };
            neighbors.push(neighbor);
        }
        
        neighbors
    }

    pub fn get_neighbor_node_id(&self, node: &Node, neighbor_id: i32) -> i32 {
        let (node1_id, node2_id) = self.edges[node.get_incident_edge(neighbor_id) as usize].get_node_ids();
        if node1_id == node.get_id() { node2_id } else { node1_id }
    }

    pub fn get_neighbor_index(&self, node: &Node, neighbor_id: i32) -> i32 {
        self.get_neighbors(node).iter().position(|id| *id == neighbor_id).expect(
            &format!("Node with id {} is not a neighbor of node with id {}", neighbor_id, node.get_id())
        ) as i32
    }

    pub fn get_neighbor_index_from_id(&self, node_id: i32, neighbor_id: i32) -> i32 {
        self.get_neighbors_from_id(node_id).iter().position(|id| *id == neighbor_id).expect(
            &format!("Node with id {} is not a neighbor of node with id {}", neighbor_id, node_id)
        ) as i32
    }

    pub fn get_peptide_for_factor(&self, factor_id: i32) -> Result<i32, Box<dyn std::error::Error>> {
        let neighbors = self.get_neighbors_from_id(factor_id);
        for neighbor_id in neighbors {
            let neighbor = self.get_node(neighbor_id);
            if let NodeType::PeptideNode { .. } = neighbor.get_subtype() {
                return Ok(neighbor.get_id());
            }
        }
        return Err(format!("Peptide not found for factor with id {}", factor_id).into());
    }

    pub fn add_ct_nodes(&mut self) {
        // When creating the CTGraph and not just reading from a previously saved graph format, use this function to add the CT nodes
        
        let mut edges_to_add: Vec<Edge> = Vec::new();
        let mut edges_to_remove: HashSet<(i32, i32)> = HashSet::new();
        let mut nodes_to_add: Vec<Node> = Vec::new();

        // Add nodes and keep track of edges to add/remove
        let mut next_edge_id: i32 = self.edges.len() as i32;
        let mut next_node_id: i32 = self.nodes.len() as i32;
        for node in &self.nodes {
            if node.is_factor_node() && node.neighbors_count() > 2 {

                let mut prot_names: Vec<String> = Vec::new();
                let mut prot_ids: Vec<i32> = Vec::new();
                for neighbor_id in self.get_neighbors(node) {
                    let neighbor: &Node = &self.nodes[neighbor_id as usize];
                    if neighbor.is_taxon_node() {
                        prot_ids.push(neighbor_id);
                        prot_names.push(neighbor.get_name().to_string());
                    }
                }
                
                // TODO: names necessary?
                let new_node_name = format!("CTree {}", prot_names.join(" "));
                let new_node_id = next_node_id;
                let new_node = Node::new_convolution_node(new_node_id, new_node_name, prot_ids.len() as i32);
                next_node_id += 1;
                nodes_to_add.push(new_node);

                let edge = Edge { id: next_edge_id, node1_id: new_node_id, node2_id: node.get_id(), message_length: Some(prot_ids.len() as i32 + 1) };
                next_edge_id += 1;
                edges_to_add.push(edge);
                for neighbor_id in prot_ids {
                    let edge = Edge { id: next_edge_id, node1_id: new_node_id, node2_id: neighbor_id, message_length: None };
                    next_edge_id += 1;
                    edges_to_add.push(edge);
                    edges_to_remove.insert((node.get_id(), neighbor_id));
                    edges_to_remove.insert((neighbor_id, node.get_id()));
                }
                
            }
        }

        // Remove edges
        let mut new_edges: Vec<Edge> = Vec::with_capacity(self.edges.len() + edges_to_add.len() - edges_to_remove.len());
        let mut next_edge_id = 0;
        for edge in &self.edges {
            if ! edges_to_remove.contains(&(edge.node1_id, edge.node2_id)) {
                let mut new_edge = edge.clone();
                new_edge.id = next_edge_id;
                next_edge_id += 1;
                new_edges.push(new_edge);
            }
        }
        // Add new edges
        for edge in edges_to_add {
            let mut new_edge = edge.clone();
            new_edge.id = next_edge_id;
            next_edge_id += 1;
            new_edges.push(new_edge);
        }

        // Update the incident edges in the nodes
        let mut new_nodes: Vec<Node> = Vec::with_capacity(self.nodes.len() + nodes_to_add.len());
        for node in &self.nodes {
            new_nodes.push(Node::new(node.get_id(), node.get_name().to_string(), node.get_subtype().clone()));
        }
        for node in nodes_to_add {
            new_nodes.push(node);
        }
        
        for edge in &new_edges {
            let (node1_id, node2_id) = edge.get_node_ids();
            new_nodes[node1_id as usize].add_incident_edge(edge.get_id());
            new_nodes[node2_id as usize].add_incident_edge(edge.get_id());
        }

        self.nodes = new_nodes;
        self.edges = new_edges;
    }

    /// Finds the connected components in an undirected graph and returns a Vec of Vecs containing nodes in each component.
    pub fn connected_components(&self) -> Vec<Self> {
        let mut visited: HashSet<i32> = HashSet::new();
        let mut components: Vec<Self> = Vec::new();

        for start_node in &self.nodes {
            if visited.insert(start_node.get_id()) {
                let mut component_ids: Vec<i32> = Vec::new();
                let mut old_to_new_nodes: HashMap<i32, i32> = HashMap::new();

                let mut new_nodes: Vec<Node> = Vec::new();
                let mut new_edges: Vec<Edge> = Vec::new();

                // Find ids of nodes to include in component
                component_ids.push(start_node.get_id());
                old_to_new_nodes.insert(start_node.get_id(), 0);
                self.find_component_rec(start_node.get_id(), &mut component_ids, &mut old_to_new_nodes, &mut visited);

                // Create new nodes
                for node_id in &component_ids {
                    let node = self.nodes[*node_id as usize].copy_with_id(old_to_new_nodes[&node_id]);
                    new_nodes.push(node);
                }

                // Select edges to keep and update the node ids
                let mut next_edge_id: i32 = 0;
                let mut component_edge_ids: HashSet<i32> = HashSet::new();
                let mut old_to_new_edges: HashMap<i32, i32> = HashMap::new();
                for edge in &self.edges {

                    let (source, target): (i32, i32) = edge.get_node_ids();
                    if component_ids.contains(&source) && component_ids.contains(&target) {

                        let (new_source, new_target): (i32, i32) = (old_to_new_nodes[&source], old_to_new_nodes[&target]);
                        let new_edge = Edge { id: next_edge_id, node1_id: new_source, node2_id: new_target, message_length: edge.get_message_length() };
                        next_edge_id += 1;

                        component_edge_ids.insert(edge.get_id());
                        old_to_new_edges.insert(edge.get_id(), new_edge.get_id());

                        new_edges.push(new_edge);
                    }
                }

                // Update edge ids of incident edges
                for node in &mut new_nodes {
                    let new_incident_edges: Vec<i32> = node.get_incident_edges().into_iter().filter(|e| component_edge_ids.contains(e)).map(|e| old_to_new_edges[e]).collect();
                    node.set_incident_edges(new_incident_edges);
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
        start_id: i32, 
        component_ids: &mut Vec<i32>, 
        old_to_new_nodes: &mut HashMap<i32, i32>, 
        visited: &mut HashSet<i32>
    ) {
        let start_node: &Node = &self.nodes[start_id as usize];
        for neighbor_id in self.get_neighbors(&start_node) {
            if visited.insert(neighbor_id) {
                let next_id: i32 = component_ids.len() as i32;
                component_ids.push(neighbor_id);
                old_to_new_nodes.insert(neighbor_id, next_id);
                self.find_component_rec(neighbor_id, component_ids, old_to_new_nodes, visited);                
            }
        }
    }
}
