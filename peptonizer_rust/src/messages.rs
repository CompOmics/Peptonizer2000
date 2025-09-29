use crate::factor_graph::CTFactorGraph;
use crate::node::{Node, NodeType};
use std::collections::{HashMap, HashSet};
use std::mem;
use crate::array_utils::*;
use crate::convolution_tree::ConvolutionTree;


/// Represents belief values stored in different node types of the factor graph.
#[derive(Debug, Clone)]
pub enum NodeBelief {
    /// Belief for peptide nodes: two probabilities (inactive, active).
    PeptideBelief(f64, f64),
    /// Belief for factor nodes: list of probability pairs.
    FactorBelief(Vec<[f64;2]>),
    /// Belief for taxon nodes: two probabilities (inactive, active).
    TaxonBelief(f64, f64),
    /// Belief placeholder for convolution tree nodes.
    ConvolutionTreeBelief
}


/// Get the initial belief for a node based on its subtype.
///
/// # Arguments
/// * `node` - Reference to the node.
///
/// # Returns
/// The corresponding `NodeBelief`.
pub fn get_initial_belief(node: &Node) -> NodeBelief {
    match node.get_subtype() {
        NodeType::PeptideNode { initial_belief_0, initial_belief_1 } => NodeBelief::PeptideBelief(*initial_belief_0, *initial_belief_1),
        NodeType::FactorNode { initial_belief, .. } => NodeBelief::FactorBelief(initial_belief.array.clone()),
        NodeType::TaxonNode { initial_belief_0, initial_belief_1 } => NodeBelief::TaxonBelief(*initial_belief_0, *initial_belief_1),
        NodeType::ConvolutionTreeNode { .. } => NodeBelief::ConvolutionTreeBelief
    }
}


impl NodeBelief {

    /// Returns all stored belief values as a flat vector.
    ///
    /// # Returns
    /// Vector of belief values.
    pub fn values(&self) -> Vec<f64> {
        match self {
            NodeBelief::PeptideBelief(a, b) => vec![*a, *b],
            NodeBelief::FactorBelief(vec) => vec.iter().flat_map(|arr| arr.to_vec()).collect(),
            NodeBelief::TaxonBelief(a, b) => vec![*a, *b],
            NodeBelief::ConvolutionTreeBelief => vec![1.0;4],
        }
    }

    /// Returns factor values if this is a factor belief.
    ///
    /// # Returns
    /// `Some(Vec<[f64; 2]>)` if factor, otherwise `None`.
    pub fn factor_values(&self) -> Option<Vec<[f64; 2]>> {
        match self {
            NodeBelief::FactorBelief(vec) => Some(vec.clone()),
            _ => None
        }
    }
}


/// Handles message passing and belief propagation in a convolution tree factor graph.
pub struct Messages {
    graph: CTFactorGraph,
    max_val: Option<(i32, i32)>,
    priorities: HashMap<(i32, i32), f64>,
    // Keeps track of residuals for duos of directed edges, indexed by [end_node][id of start node in end neighbours][id of neighbour in end node]
    total_residuals: Vec<Vec<Vec<f64>>>, 
    // Maps a node ID onto its current belief value
    current_beliefs: Vec<NodeBelief>, 
    // incoming messages for each node [end node][neighbour id]
    msg_in: Vec<Vec<Vec<f64>>>,
    msg_in_new: Vec<Vec<Vec<f64>>>,
    msg_in_log: Vec<Vec<Vec<f64>>>
}


impl Messages {

    /// Constructs a new `Messages` instance for the given graph.
    ///
    /// # Arguments
    /// * `ct_graph_in` - Factor graph to initialize messages from.
    ///
    /// # Returns
    /// Initialized `Messages` object.
    pub fn new(ct_graph_in: CTFactorGraph) -> Messages {
        let max_val: Option<(i32, i32)> = None;
        
        let priorities = HashMap::new();

        let mut total_residuals: Vec<Vec<Vec<f64>>> = Vec::with_capacity(ct_graph_in.node_count());
        for node in ct_graph_in.get_nodes() {
            
            let total_residual_node: Vec<Vec<f64>> = vec![vec![0.0; node.neighbors_count()]; node.neighbors_count()];

            total_residuals.push(total_residual_node);

        }

        let mut current_beliefs: Vec<NodeBelief> = Vec::with_capacity(ct_graph_in.node_count());
        for node in ct_graph_in.get_nodes() {
            current_beliefs.push(get_initial_belief(node));
        }

        let mut msg_in = Vec::with_capacity(ct_graph_in.node_count());
        let mut msg_in_new = Vec::with_capacity(ct_graph_in.node_count());
        for node in ct_graph_in.get_nodes() {
            
            let mut msg_in_node = Vec::with_capacity(node.neighbors_count());
            let mut msg_in_new_node = Vec::with_capacity(node.neighbors_count());
            for edge_id in node.get_incident_edges() {
                match ct_graph_in.get_edge(*edge_id).get_message_length() {
                    Some(message_length) => {
                        msg_in_node.push(vec![1.0; message_length as usize]);
                        msg_in_new_node.push(vec![1.0; message_length as usize]);
                    },
                    None => {
                        msg_in_node.push(vec![0.5, 0.5]);
                        msg_in_new_node.push(vec![0.0, 0.0]);
                    }
                }
            }

            msg_in.push(msg_in_node);
            msg_in_new.push(msg_in_new_node);
        }

        let msg_in_log = msg_in.clone();

        Messages { graph: ct_graph_in, max_val, priorities, total_residuals, current_beliefs, msg_in, msg_in_new, msg_in_log }

    }

    /// Runs zero-lookahead belief propagation until convergence or limit.
    ///
    /// # Arguments
    /// * `max_loops` - Maximum number of iterations.
    /// * `tolerance` - Residual threshold for convergence.
    ///
    /// # Returns
    /// Final node beliefs as a list of value vectors.
    pub fn zero_lookahead_bp(&mut self, max_loops: i32, tolerance: f64) -> Result<Vec<Vec<f64>>, Box<dyn std::error::Error>> {

        let mut max_residual: f64 = f64::MAX;

        // first, do 5 loops where I update all messages
        for _ in 0..5 {
            self.compute_update(false)?;

            let temp = mem::replace(&mut self.msg_in_log, mem::take(&mut self.msg_in));
            self.msg_in = mem::take(&mut self.msg_in_new);
            self.msg_in_new = temp;
        }

        // compute all residuals after 5 runs once (= initialize the residual/priorities vectors)
        for node_id in 0..self.graph.node_count() {
            let node = self.graph.get_node(node_id as i32);
            for neighbor_id in 0..node.neighbors_count() {
                let residual = self.compute_infinity_norm_residual(node_id, neighbor_id);
                self.priorities.insert((node_id as i32, neighbor_id as i32), residual);
            }
        }

        let mut k = 5;

        // keep track of the nodes of which the incoming messages have changed
        let mut prev_changed: Vec<usize> = (0..self.msg_in.len()).collect();

        while k < max_loops && max_residual > tolerance {

            // actual zero-look-ahead-BP part
            let (&(end_id, start_in_end_id), &residual) = self.priorities.iter().max_by(|a, b| a.1.partial_cmp(b.1).expect("Partial compare returned None")).ok_or("Priorities is empty")?;
            max_residual = residual;
            
            let end_node = self.graph.get_node(end_id);
            let start_id = self.graph.get_neighbor_node_id(end_node, start_in_end_id);
            
            let end_in_start_id = self.graph.get_neighbor_index_from_id(start_id, end_id);
            
            self.single_message_update(start_id, end_id, end_in_start_id, None)?;

            let priority_residual = self.compute_infinity_norm_residual(end_id as usize, start_in_end_id as usize);
            for node_id in prev_changed {
                let msg: Vec<Vec<f64>> = self.msg_in[node_id].clone();
                self.msg_in_log[node_id] = msg;
            }
            prev_changed = Vec::new();

            // if the start node is a convolution tree, all the incoming messages of the neighbours can be changed.
            let start_node = self.graph.get_node(start_id);
            match start_node.get_subtype() {
                NodeType::ConvolutionTreeNode { .. } => {
                    // TODO: shouldn't on both levels be added to prev_changed and change in msg_in?
                    for neighbor_id in self.graph.get_neighbors(start_node) {
                        let neighbor = self.graph.get_node(neighbor_id);
                        prev_changed.push(neighbor_id as usize);
                        for i in 0..neighbor.neighbors_count() {
                            self.msg_in[neighbor_id as usize][i] = self.msg_in_new[neighbor_id as usize][i].clone();
                        }
                    }
                },
                _ => {
                    self.msg_in[end_id as usize][start_in_end_id as usize] = self.msg_in_new[end_id as usize][start_in_end_id as usize].clone();
                    prev_changed.push(end_id as usize);
                }
            }

            self.compute_total_residuals(start_id, end_id, start_in_end_id, priority_residual);
            self.compute_priority(start_id, end_id, start_in_end_id)?;

            k += 1;
        }

        // marginalize once the model has converged
        for node in self.graph.get_nodes() {
            match node.get_subtype() {
                NodeType::TaxonNode { .. } | NodeType::PeptideNode { .. } => {
                    let incoming_messages: Vec<Vec<f64>> = self.msg_in[node.get_id() as usize].clone();

                    // need for logs to prevent underflow in very large multiplications
                    let incoming_messages_log: Vec<Vec<f64>> = incoming_messages.iter()
                        .map(|msg| msg.iter().map(|&x| x.ln()).collect())
                        .collect();
                    
                    let initial_belief_log: Vec<f64> = get_initial_belief(node).values().iter()
                        .map(|&b| b.ln())
                        .collect();

                    let sum_logs: Vec<f64> = incoming_messages_log.iter()
                        .fold(vec![0.0;2], |mut acc,  row| {acc[0] += row[0]; acc[1] += row[1]; acc});

                    // Compute final log-normalized message
                    let mut logged_variable_marginal: Vec<f64> = initial_belief_log.iter()
                        .zip(sum_logs.iter())
                        .map(|(&a, &b)| a + b)
                        .collect();
                    log_normalize(&mut logged_variable_marginal);

                    self.current_beliefs[node.get_id() as usize] = match node.get_subtype() {
                        NodeType::TaxonNode { .. } => NodeBelief::TaxonBelief(logged_variable_marginal[0], logged_variable_marginal[1]),
                        _ => NodeBelief::PeptideBelief(logged_variable_marginal[0], logged_variable_marginal[1]),
                    }
                },
                _ => {}
            }
        }

        Ok(self.current_beliefs.iter().map(|b| b.values()).collect())
    }

    /// Updates all outgoing messages from all nodes.
    ///
    /// # Arguments
    /// * `local_loops` - If true, update only local region around last max residual.
    fn compute_update(&mut self, local_loops: bool) -> Result<(), Box<dyn std::error::Error>> {
        let mut checked_cts: HashSet<i32> = HashSet::new();
        
        if self.max_val.is_some() && local_loops {

            let (_, start_id) = self.max_val.unwrap(); // Safe unwrap because of is_some() check
            let start_node = self.graph.get_node(start_id);
            for (end_in_start_id, end_id) in self.graph.get_neighbors(start_node).iter().enumerate() {
                self.single_message_update(start_id, *end_id, end_in_start_id as i32, Some(&mut checked_cts))?;
            }

        } else {

            for id in 0..self.graph.node_count() {
                let id = id as i32;
                let start_node = self.graph.get_node(id);
                for (end_in_start_id, end_id) in self.graph.get_neighbors(start_node).iter().enumerate() {
                    self.single_message_update(id, *end_id, end_in_start_id as i32, Some(&mut checked_cts))?;
                }
            }
        }

        Ok(())
    }

    /// Updates a single message along an edge.
    ///
    /// # Arguments
    /// * `start_id` - ID of source node.
    /// * `end_id` - ID of destination node.
    /// * `end_in_start_id` - Neighbor index of destination in source.
    /// * `checked_cts` - Optional set of already updated convolution tree nodes.
    fn single_message_update(&mut self, start_id: i32, end_id: i32, end_in_start_id: i32, mut checked_cts: Option<&mut HashSet<i32>>) -> Result<(), Box<dyn std::error::Error>> {

        let start_node: &Node = self.graph.get_node(start_id);

        match start_node.get_subtype() {
            NodeType::TaxonNode { .. } | NodeType::PeptideNode { .. } => {
                let new_message = self.compute_out_message_variable(start_id, end_id, end_in_start_id);
                let end_node: &Node = self.graph.get_node(end_id);
                let start_in_end_id: i32 = self.graph.get_neighbor_index(end_node, start_id);
                self.msg_in_new[end_id as usize][start_in_end_id as usize] = new_message;
            },
            NodeType::FactorNode { .. } => {
                let new_message = self.compute_out_message_factor(start_id, end_id, end_in_start_id)?;
                let end_node: &Node = self.graph.get_node(end_id);
                let start_in_end_id: i32 = self.graph.get_neighbor_index(end_node, start_id);
                self.msg_in_new[end_id as usize][start_in_end_id as usize] = new_message;
            },
            NodeType::ConvolutionTreeNode { number_of_parents } => 
                if checked_cts.as_ref().map_or(true, |set| ! set.contains(&start_id)) {
                    self.compute_out_messages_ct_tree(start_id, *number_of_parents)?;
                    if let Some(set) = checked_cts.as_mut() {
                        set.insert(start_id);
                    }
                }
        };

        Ok(())
    }

    /// Computes outgoing message for variable (peptide/taxon) nodes.
    ///
    /// # Arguments
    /// * `start_id` - ID of source node.
    /// * `end_id` - ID of destination node.
    /// * `end_in_start_id` - Neighbor index of destination in source.
    /// 
    /// # Returns
    /// Normalized probability vector.
    fn compute_out_message_variable(&mut self, start_id: i32, end_id: i32, end_in_start_id: i32) -> Vec<f64> {
        
        let start_node = self.graph.get_node(start_id);
        let end_node = self.graph.get_node(end_id);
        let mut incoming_messages_end: Vec<Vec<f64>> = self.msg_in[start_id as usize].clone();
        incoming_messages_end.remove(end_in_start_id as usize);
        let node_belief: Vec<f64> = self.current_beliefs[start_id as usize].values();

        if incoming_messages_end.is_empty() {
            // TODO: check in Python code the idea of the any statement (ask Tanja)
            let start_in_end_id: i32 = self.graph.get_neighbor_index(end_node, start_id);
            if node_belief == get_initial_belief(start_node).values() { 
                return node_belief; 
            } else { 
                return self.msg_in[end_id as usize][start_in_end_id as usize].clone();
            }
        }

        // need for logs to prevent underflow in very large multiplications
        let incoming_messages_end: Vec<Vec<f64>> = incoming_messages_end
        .iter()
        .map(|msg| msg.iter().map(|&x| x.ln()).collect())
        .collect();

        // Sum of incoming messages
        let sum_logs: Vec<f64> = incoming_messages_end.iter().fold(vec![0.0;2], |mut acc,  row| {acc[0] += row[0]; acc[1] += row[1]; acc});

        // Log transform node_belief
        let node_belief_log: Vec<f64> = node_belief.iter().map(|&x| x.ln()).collect();

        // Compute final log-normalized message
        let mut out_message_log: Vec<f64> = node_belief_log.iter().zip(sum_logs.iter()).map(|(&a, &b)| a + b).collect();
        log_normalize(&mut out_message_log);

        // Prevent underflow: Replace zeros with 1e-30
        avoid_underflow(&mut out_message_log);

        out_message_log
    }

    /// Computes outgoing message for factor nodes.
    ///
    /// # Arguments
    /// * `start_id` - ID of source node.
    /// * `end_id` - ID of destination node.
    /// * `end_in_start_id` - Neighbor index of destination in source.
    /// 
    /// # Returns
    /// Normalized probability vector.
    fn compute_out_message_factor(&mut self, start_id: i32, end_id: i32, end_in_start_id: i32) -> Result<Vec<f64>, Box<dyn std::error::Error>> {
        let end_node = self.graph.get_node(end_id);
        let mut incoming_messages_end: Vec<Vec<f64>> = self.msg_in[start_id as usize].clone();
        incoming_messages_end.remove(end_in_start_id as usize);
        let node_belief: Vec<[f64;2]> = self.current_beliefs[start_id as usize].factor_values().ok_or("factor_values called on a NodeBelief which is not a FactorBelief")?;

        match end_node.get_subtype() {
            NodeType::ConvolutionTreeNode { .. } => {
                // handles empty & messages with only one value
                incoming_messages_end.push(vec![1.0, 1.0]);
                
                // Product of incoming messages
                let prod: [f64;2] = incoming_messages_end.iter().fold([1.0;2], |mut acc,  row| {acc[0] *= row[0]; acc[1] *= row[1]; acc});
                
                // Compute final normalized message
                let mut out_message: Vec<[f64;2]> = node_belief.iter().map(|&a| [a[0] * prod[0], a[1] * prod[1]]).collect();
                normalize_2d(&mut out_message);
                let out_message_sum: Vec<f64> = out_message.iter().map(|[a, b]| a + b).collect();

                return Ok(out_message_sum);
            },
            _ => {
                if incoming_messages_end[0].len() > 2 {
                    let incoming_messages_log: Vec<f64> = incoming_messages_end[0].iter().map(|&x| x.ln()).collect();
                    let node_belief_log: Vec<[f64;2]> = node_belief.iter().map(|&x| [x[0].ln(), x[1].ln()]).collect();

                    let mut out_message_log: Vec<[f64;2]> = node_belief_log.iter().zip(incoming_messages_log.iter()).map(|(&a, &b)| [a[0] + b, a[1] + b]).collect();
                    log_normalize_2d(&mut out_message_log);

                    let mut out_message_sum: Vec<f64> = out_message_log.iter().take(2).map(|[a, b]| a + b).collect();
                    // Prevent underflow: Replace zeros with 1e-30
                    avoid_underflow(&mut out_message_sum);

                    return Ok(out_message_sum);

                } else {

                    incoming_messages_end.push(vec![1.0, 1.0]);

                    let prod: [f64;2] = incoming_messages_end.iter().fold([1.0;2], |mut acc,  row| {acc[0] *= row[0]; acc[1] *= row[1]; acc});

                    // Compute final normalized message
                    let mut out_message: Vec<[f64;2]> = node_belief.iter().map(|&a| [a[0] * prod[0], a[1] * prod[1]]).collect();
                    normalize_2d(&mut out_message);
                    let out_message_sum: Vec<f64> = out_message.iter().map(|[a, b]| a + b).collect();

                    return Ok(out_message_sum);
                }
            }
        }
    }

    /// Computes outgoing messages for convolution tree nodes.
    ///
    /// # Arguments
    /// * `start_id` - Node ID of convolution tree.
    /// * `number_of_parents` - Number of parent nodes.
    fn compute_out_messages_ct_tree(&mut self, start_id: i32, number_of_parents: i32) -> Result<(), Box<dyn std::error::Error>> {
        let start_node = self.graph.get_node(start_id);

        let mut prot_prob_list: Vec<Vec<f64>> = Vec::new();
        let mut old_prot_prob_list: Vec<Vec<f64>> = Vec::new();

        let mut shared_likelihoods: Vec<f64> = vec![1.0; number_of_parents as usize + 1];
        let mut old_shared_likelihoods: Vec<f64> = Vec::new();

        let mut peptides: Vec<i32> = Vec::new();
        let mut prot_list: Vec<i32> = Vec::new();

        let mut last_neighbor_id: Option<usize> = None;
        for (neighbor_id_in_start, &neighbor_id) in self.graph.get_neighbors(start_node).iter().enumerate() {
            let neighbor = self.graph.get_node(neighbor_id);
            match neighbor.get_subtype() {
                NodeType::FactorNode { .. } => {
                    peptides.push(neighbor_id);
                    
                    // Prevent underflow: Replace zeros with 1e-30
                    avoid_underflow(&mut shared_likelihoods);
                    for (a, b) in shared_likelihoods.iter_mut().zip(self.msg_in[start_id as usize][neighbor_id_in_start].iter()) {
                        *a *= b;
                    }

                    last_neighbor_id = Some(neighbor_id_in_start);
                },
                _ => {
                    prot_prob_list.push(self.msg_in[start_id as usize][neighbor_id_in_start].clone());
                    old_prot_prob_list.push(self.msg_in_log[start_id as usize][neighbor_id_in_start].clone());
                    prot_list.push(neighbor_id);
                }
            }
        }

        if let Some(neighbor_id) = last_neighbor_id {
            // Prevent underflow: Replace zeros with 1e-30
            avoid_underflow(&mut shared_likelihoods);
            old_shared_likelihoods = shared_likelihoods.iter()
                .zip(self.msg_in_log[start_id as usize][neighbor_id].iter())
                .map(|(a, b)| a * b).collect();
        }

        if old_shared_likelihoods != shared_likelihoods && 
            prot_prob_list.iter().zip(old_prot_prob_list.iter()).any(|(a, b)| a[0] != b[0]) {
            let convolution_tree = ConvolutionTree::new(shared_likelihoods, prot_prob_list)?;

            for (protein_id, protein) in prot_list.iter().enumerate() {
                let node_neighbor_index: i32 = self.graph.get_neighbor_index_from_id(*protein, start_id);
                let mut new_message = convolution_tree.message_to_variable(protein_id);
                avoid_underflow(&mut new_message);
                self.msg_in_new[*protein as usize][node_neighbor_index as usize] = new_message;
            }

            for pep in peptides {
                let node_neighbor_index: i32 = self.graph.get_neighbor_index_from_id(pep, start_id);
                let mut new_message = convolution_tree.message_to_shared_likelihood()?;
                avoid_underflow(&mut new_message);
                self.msg_in_new[pep as usize][node_neighbor_index as usize] = new_message;
            }

        } else {
            
            for protein in prot_list {
                let node_neighbor_index: i32 = self.graph.get_neighbor_index_from_id(protein, start_id);
                self.msg_in_new[protein as usize][node_neighbor_index as usize] = self.msg_in[protein as usize][node_neighbor_index as usize].clone();
            }

            for pep in peptides {
                let node_neighbor_index: i32 = self.graph.get_neighbor_index_from_id(pep, start_id);
                self.msg_in_new[pep as usize][node_neighbor_index as usize] = self.msg_in[pep as usize][node_neighbor_index as usize].clone();
            }
        }

        Ok(())
    }

    /// Computes infinity norm residual for an edge.
    ///
    /// # Arguments
    /// * `end_id` - ID of destination node.
    /// * `end_in_start_id` - Neighbor index of destination in source.
    /// 
    /// # Returns
    /// Residual as `f64`.
    fn compute_infinity_norm_residual(&mut self, end_id: usize, start_in_end_id: usize) -> f64 {
        let msg1: &mut Vec<f64> = &mut self.msg_in[end_id][start_in_end_id];
        let msg2: &mut Vec<f64> = &mut self.msg_in_log[end_id][start_in_end_id];

        if msg1.len() != msg2.len() {
            *msg2 = vec![1.0; msg1.len()];
        }

        avoid_underflow(msg1);

        let residual = msg1.iter()
                           .zip(msg2.iter())
                           .map(|(m1, m2)| (m1 / m2).ln().abs())
                           .fold(f64::NEG_INFINITY, f64::max);

        residual
    }

    /// Updates total residuals after a message update.
    /// 
    /// # Arguments
    /// * `start_id` - ID of source node.
    /// * `end_id` - ID of destination node.
    /// * `end_in_start_id` - Neighbor index of destination in source.
    /// * `current_residual` - Redidual to add to message
    fn compute_total_residuals(&mut self, start_id: i32, end_id: i32, start_in_end_id: i32, current_residual: f64) {
        let start_node = self.graph.get_node(start_id);
        let end_node = self.graph.get_node(end_id);

        let end_in_start_id: i32 = self.graph.get_neighbor_index(start_node, end_id);

        for (i, &neighbor_id) in self.graph.get_neighbors(start_node).iter().enumerate() {
            if neighbor_id != end_id {
                self.total_residuals[start_id as usize][i][end_in_start_id as usize] = 0.0;
            }
        }

        for (i, &neighbor_id) in self.graph.get_neighbors(end_node).iter().enumerate() {
            if neighbor_id != start_id {
                self.total_residuals[end_id as usize][start_in_end_id as usize][i] += current_residual;
            }
        }
    }

    /// Updates message priorities for scheduling.
    /// 
    /// # Arguments
    /// * `start_id` - ID of source node.
    /// * `end_id` - ID of destination node.
    /// * `end_in_start_id` - Neighbor index of destination in source.
    fn compute_priority(&mut self, start_id: i32, end_id: i32, start_in_end_id: i32) -> Result<(), Box<dyn std::error::Error>> {
        let end_node = self.graph.get_node(end_id);

        let priority = self.priorities.get_mut(&(end_id, start_in_end_id)).ok_or("Index out-of-bounds")?;
        *priority = 0.0;

        for (i, &neighbor_id) in self.graph.get_neighbors(end_node).iter().enumerate() {
            if neighbor_id != start_id {
                let neighbor_node = self.graph.get_node(neighbor_id);
                let end_in_neighbor_id: i32 = self.graph.get_neighbor_index(neighbor_node, end_id);
                let priority = self.priorities.get_mut(&(neighbor_id, end_in_neighbor_id)).ok_or("Index out-of-bounds")?;
                *priority = self.graph
                    .get_neighbors(end_node)
                    .iter()
                    .enumerate()
                    .map(|(j, &sum_run)| {
                        if sum_run != neighbor_id { 
                            self.total_residuals[end_id as usize][j][i] 
                        } else { 
                            0.0
                        }
                    }).sum();
            }
        }

        Ok(())
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::node::{Node, NodeType, Factor};
    use crate::factor_graph::CTFactorGraph;
    use std::collections::HashSet;
    use crate::factor_graph::Edge;

    /// Creates a minimal graph with Peptide -> Factor -> Taxon
    fn create_minimal_graph() -> CTFactorGraph {
        let mut nodes: Vec<Node> = Vec::new();
        let mut edges: Vec<Edge> = Vec::new();

        let peptide_node_1 = Node::new(
            0,
            "peptide1".to_string(),
            NodeType::PeptideNode { initial_belief_0: 0.7, initial_belief_1: 0.3 }
        );
        nodes.push(peptide_node_1);

        let factor_node_1 = Node::new(
            1,
            "factor1".to_string(),
            NodeType::FactorNode { 
                parent_number: 2, 
                initial_belief: Factor { 
                    array: vec![[0.6, 0.4], [0.2, 0.8]], 
                    array_labels: vec!["p0".to_string(), "p1".to_string()] 
                } 
            }
        );
        nodes.push(factor_node_1);

        let taxon_node_1 = Node::new(
            2,
            "taxon_1".to_string(),
            NodeType::TaxonNode { initial_belief_0: 0.5, initial_belief_1: 0.5 }
        );
        nodes.push(taxon_node_1);

        let edge0 = Edge::new(0, 0, 1, None);
        nodes[0].add_incident_edge(0);
        nodes[1].add_incident_edge(0);
        edges.push(edge0);


        let edge1 = Edge::new(1, 1, 2, None);
        nodes[1].add_incident_edge(1);
        nodes[2].add_incident_edge(1);
        edges.push(edge1);



        CTFactorGraph::new(nodes, edges)
}

    #[test]
    fn test_get_initial_belief() {
        let peptide_node = Node::new(
            0,
            "peptide1".to_string(),
            NodeType::PeptideNode { initial_belief_0: 0.6, initial_belief_1: 0.4 }
        );

        if let NodeType::PeptideNode { initial_belief_0, initial_belief_1 } = peptide_node.get_subtype() {
            assert!((initial_belief_0 - 0.6).abs() < 1e-10);
            assert!((initial_belief_1 - 0.4).abs() < 1e-10);
        } else {
            panic!("Expected PeptideNode");
        }

        let factor_node = Node::new(
            1,
            "factor1".to_string(),
            NodeType::FactorNode {
                parent_number: 2,
                initial_belief: Factor { array: vec![[0.5, 0.5], [0.3, 0.7]], array_labels: vec!["p0".to_string(), "p1".to_string()] }
            }
        );

        if let NodeType::FactorNode { parent_number, initial_belief } = factor_node.get_subtype() {
            assert_eq!(*parent_number, 2);
            assert_eq!(initial_belief.array.len(), 2);
            assert!((initial_belief.array[0][0] - 0.5).abs() < 1e-10);
            assert!((initial_belief.array[1][1] - 0.7).abs() < 1e-10);
            assert_eq!(initial_belief.array_labels, vec!["p0".to_string(), "p1".to_string()]);
        } else {
            panic!("Expected FactorNode");
        }

        let taxon_node = Node::new(
            2,
            "taxon1".to_string(),
            NodeType::TaxonNode { initial_belief_0: 0.3, initial_belief_1: 0.7 }
        );

        if let NodeType::TaxonNode { initial_belief_0, initial_belief_1 } = taxon_node.get_subtype() {
            assert!((initial_belief_0 - 0.3).abs() < 1e-10);
            assert!((initial_belief_1 - 0.7).abs() < 1e-10);
        } else {
            panic!("Expected TaxonNode");
        }

        let ct_node = Node::new_convolution_node(3, "ct1".to_string(), 4);

        if let NodeType::ConvolutionTreeNode { number_of_parents } = ct_node.get_subtype() {
            assert_eq!(*number_of_parents, 4);
        } else {
            panic!("Expected ConvolutionTreeNode");
        }
    }

    #[test]
    fn test_nodebelief_values_and_factor_values() {
        let pb = NodeBelief::PeptideBelief(0.1,0.9);
        assert_eq!(pb.values(), vec![0.1,0.9]);
        let fb = NodeBelief::FactorBelief(vec![[0.2,0.8]]);
        assert_eq!(fb.values(), vec![0.2,0.8]);
        assert_eq!(fb.factor_values(), Some(vec![[0.2,0.8]]));
        let tb = NodeBelief::TaxonBelief(0.3,0.7);
        assert_eq!(tb.values(), vec![0.3,0.7]);
        let cb = NodeBelief::ConvolutionTreeBelief;
        assert_eq!(cb.values(), vec![1.0;4]);
    }

    #[test]
    fn test_messages_zero_lookahead_bp() {
        let graph = create_minimal_graph();
        let mut messages = Messages::new(graph);
        let beliefs = messages.zero_lookahead_bp(5,1e-6);
        assert_eq!(beliefs[0].len(),2);
        assert_eq!(beliefs[2].len(),2);
        let sum: f64 = beliefs[0].iter().sum();
        assert!((sum-1.0).abs()<1e-6);
    }

    #[test]
    fn test_compute_out_message_variable() {
        let graph = create_minimal_graph();
        let mut messages = Messages::new(graph);
        let msg = messages.compute_out_message_variable(0,1,0);
        assert_eq!(msg.len(),2);
        let s: f64 = msg.iter().sum();
        assert!((s-1.0).abs()<1e-6);
    }

    #[test]
    fn test_compute_out_message_factor() {
        let graph = create_minimal_graph();
        let mut messages = Messages::new(graph);
        let msg = messages.compute_out_message_factor(1,2,0);
        assert_eq!(msg.len(),2);
        let s: f64 = msg.iter().sum();
        assert!((s-1.0).abs()<1e-6);
    }

    #[test]
    fn test_compute_infinity_norm_residual_and_total_residuals() {
        let graph = create_minimal_graph();
        let mut messages = Messages::new(graph);
        let residual = messages.compute_infinity_norm_residual(0,0);
        assert!(residual >= 0.0);
        messages.compute_total_residuals(0,1,0,0.1);
        assert!(messages.total_residuals[1][0][1] > 0.0);
    }

    #[test]
    fn test_compute_priority() {
        let graph = create_minimal_graph();
        let mut messages = Messages::new(graph);
        for node_id in 0..6 {
            for neighbor_id in 0..2 {
                messages.priorities.insert((node_id, neighbor_id), 0.1);
            }
        }
        messages.compute_priority(0,2,0);
        assert!(messages.priorities.get(&(1,0)).is_some());
    }

    #[test]
    fn test_single_message_update_and_compute_update() {
        let graph = create_minimal_graph();
        let mut messages = Messages::new(graph);
        let mut checked = HashSet::new();
        messages.single_message_update(0,1,0,Some(&mut checked)).unwrap();
        assert!(checked.is_empty() || checked.contains(&0)==false);
        messages.compute_update(false).unwrap(); 
    }
}
