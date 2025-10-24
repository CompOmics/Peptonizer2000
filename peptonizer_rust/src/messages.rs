use crate::factor_graph::CTFactorGraph;
use crate::node::{Node, NodeType};
use std::collections::HashSet;
use std::mem;
use crate::array_utils::*;
use crate::convolution_tree::ConvolutionTree;
use priority_queue::PriorityQueue;
use ordered_float::OrderedFloat;
use crate::array_utils::sum_logs_batched;

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

    /// Returns the belief values as a fixed-size array `[f64; 2]` for variable nodes (peptide or taxon).
    ///
    /// # Returns
    /// * `Some([f64; 2])` if the node is a peptide or taxon node.
    /// * `None` if the node is a factor or convolution tree node.
    pub fn variable_values(&self) -> Option<[f64; 2]> {
        match self {
            NodeBelief::PeptideBelief(a, b) | NodeBelief::TaxonBelief(a, b) => Some([*a, *b]),
            _ => None
        }
    }

    /// Returns factor values if this is a factor belief.
    ///
    /// # Returns
    /// `Some(Vec<[f64; 2]>)` if factor, otherwise `None`.
    pub fn factor_values(&self) -> Option<&Vec<[f64; 2]>> {
        match self {
            NodeBelief::FactorBelief(vec) => Some(vec),
            _ => None
        }
    }
}


#[derive(Clone, Debug)]
pub enum MessagesInNode {
    MessagesInVariable {
        messages: Vec<[f64; 2]>,
    },
    MessagesInFactor {
        ct_message: Option<Vec<f64>>,
        variable_messages: Vec<[f64; 2]>,
    },
    MessagesInCTree {
        factor_message: Vec<f64>,
        variable_messages: Vec<[f64; 2]>,
    },
}

impl MessagesInNode {
    
    /// Returns all variable messages stored in the node.
    ///
    /// # Returns
    /// Reference to a vector of `[f64; 2]` messages.
    pub fn get_messages(&self) -> &Vec<[f64;2]> {
        match self {
            MessagesInNode::MessagesInVariable { messages } => messages,
            MessagesInNode::MessagesInFactor { variable_messages, .. } => variable_messages,
            MessagesInNode::MessagesInCTree { variable_messages, .. } => variable_messages,
        }
    }

    /// Sets the message for a specific neighbor.
    ///
    /// # Arguments
    /// * `neighbor_index` - Index of the neighbor.
    /// * `message` - New message to store.
    pub fn set_message(&mut self, neighbor_index: usize, message: [f64; 2]) {
        match self {
            MessagesInNode::MessagesInVariable { messages } => messages[neighbor_index] = message,
            MessagesInNode::MessagesInFactor { variable_messages, ct_message } => {
                if ct_message.is_none() {
                    variable_messages[neighbor_index] = message;
                } else {
                    variable_messages[neighbor_index-1] = message;
                }
            }
            MessagesInNode::MessagesInCTree { variable_messages, .. } => variable_messages[neighbor_index-1] = message,
        }
    }

    /// Retrieves the message for a specific neighbor.
    ///
    /// # Arguments
    /// * `neighbor_index` - Index of the neighbor.
    ///
    /// # Returns
    /// Message as `[f64; 2]`.
    pub fn get_message(&self, neighbor_index: usize) -> [f64;2] {
        match self {
            MessagesInNode::MessagesInVariable { messages } => messages[neighbor_index],
            MessagesInNode::MessagesInFactor { variable_messages, ct_message } => {
                if ct_message.is_none() {
                    return variable_messages[neighbor_index];
                } else {
                    return variable_messages[neighbor_index - 1];
                }
            },
            MessagesInNode::MessagesInCTree { variable_messages, .. } => variable_messages[neighbor_index - 1],
        }
    }

    /// Returns the number of messages stored in the node.
    ///
    /// # Returns
    /// Message count as `usize`.
    pub fn get_message_count(&self) -> usize {
        self.get_messages().len()
    }

    /// Checks whether a message at a given index is a convolution tree message.
    ///
    /// # Arguments
    /// * `message_index` - Index of the message.
    ///
    /// # Returns
    /// `true` if it is a convolution tree message, `false` otherwise.
    pub fn is_ct_message(&self, message_index: usize) -> bool {
        if message_index != 0 {
            return false;
        }
        match self {
            MessagesInNode::MessagesInCTree { .. } => true,
            MessagesInNode::MessagesInFactor { ct_message: Some(_), .. } => true,
            _ => false,
        }
    }

    /// Retrieves the convolution tree message if it exists.
    ///
    /// # Returns
    /// `Ok(&Vec<f64>)` if message exists, error otherwise.
    pub fn get_ct_message(&self) -> Result<&Vec<f64>, Box<dyn std::error::Error>> {
        match self {
            MessagesInNode::MessagesInCTree { factor_message, .. } => Ok(factor_message),
            MessagesInNode::MessagesInFactor { ct_message, .. } => {
                if let Some(mes) = ct_message {
                    return Ok(mes);
                }
                Err("Factor Node does not have a connected CTNode".into())
            },
            _ => Err("get_ct_message is only valid on MessageInCTree and MessagesInFactor".into())
        }
    }

    /// Sets the convolution tree message.
    ///
    /// # Arguments
    /// * `message` - Message to set.
    ///
    /// # Returns
    /// `Ok(())` on success, error if node type is incompatible.
    pub fn set_ct_message(&mut self, message: Vec<f64>) -> Result<(), Box<dyn std::error::Error>> {
        match self {
            MessagesInNode::MessagesInCTree { factor_message, .. } => {
                *factor_message = message;
                return Ok(());
            },
            MessagesInNode::MessagesInFactor { ct_message, .. } => {
                *ct_message = Some(message);
                return Ok(());
            },
            _ => Err("set_ct_message is only valid on MessageInCTree and MessagesInFactor".into())
        }
    }
}


/// Handles message passing and belief propagation in a convolution tree factor graph.
pub struct Messages {
    graph: CTFactorGraph,
    max_val: Option<(usize, usize)>,
    priorities: PriorityQueue<(u32, u32), OrderedFloat<f64>>,
    // Keeps track of residuals for duos of directed edges, indexed by [end_node][id of start node in end neighbours][id of neighbour in end node]
    total_residuals: Vec<Vec<Vec<f64>>>, 
    // Maps a node ID onto its current belief value
    current_beliefs: Vec<NodeBelief>, 
    // incoming messages for each node [end node][neighbour id]
    msg_in: Vec<MessagesInNode>,
    msg_in_new: Vec<MessagesInNode>,
    msg_in_log: Vec<MessagesInNode>
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
        let max_val: Option<(usize, usize)> = None;
        
        let priorities = PriorityQueue::new();

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

            match node.get_subtype() {
                NodeType::PeptideNode { .. } | NodeType::TaxonNode { .. } => {
                    msg_in.push(MessagesInNode::MessagesInVariable { messages: vec![[0.5, 0.5]; node.neighbors_count()] });
                    msg_in_new.push(MessagesInNode::MessagesInVariable { messages: vec![[0.0, 0.0]; node.neighbors_count()]});
                },
                NodeType::FactorNode { .. } => {
                    let ct_message_length = ct_graph_in.get_edge(node.get_incident_edge(0)).get_message_length();
                    let mut variable_messages = vec![[0.5, 0.5]; node.neighbors_count()];
                    let mut variable_messages_new = vec![[0.0, 0.0]; node.neighbors_count()];
                    let mut ct_message = None;
                    let mut ct_message_new = None;

                    if let Some(message_length) = ct_message_length {
                        variable_messages = vec![[0.5, 0.5]; node.neighbors_count() - 1];
                        variable_messages_new = vec![[0.0, 0.0]; node.neighbors_count() - 1];
                        ct_message = Some(vec![1.0; message_length]);
                        ct_message_new = Some(vec![1.0; message_length]);
                    }

                    msg_in.push(MessagesInNode::MessagesInFactor { variable_messages, ct_message });
                    msg_in_new.push(MessagesInNode::MessagesInFactor { variable_messages: variable_messages_new, ct_message: ct_message_new });
                },
                NodeType::ConvolutionTreeNode { .. } => {
                    let variable_messages = vec![[0.5, 0.5]; node.neighbors_count() - 1];
                    let variable_messages_new = vec![[0.0, 0.0]; node.neighbors_count() - 1];

                    let factor_message_length = ct_graph_in.get_edge(node.get_incident_edge(0))
                        .get_message_length().expect("First incident edge should be the edge to Factor Node");

                    let factor_message = vec![1.0; factor_message_length];
                    let factor_message_new = vec![1.0; factor_message_length];

                    msg_in.push(MessagesInNode::MessagesInCTree { variable_messages: variable_messages, factor_message: factor_message });
                    msg_in_new.push(MessagesInNode::MessagesInCTree { variable_messages: variable_messages_new, factor_message: factor_message_new });
                }
            }
        }

        let msg_in_log = msg_in.clone();

        Messages { graph: ct_graph_in, max_val, priorities, total_residuals, current_beliefs, msg_in, msg_in_new, msg_in_log }

    }

    /// Adds a message priority to the scheduling queue.
    ///
    /// # Arguments
    /// * `node_id` - Node ID of the message.
    /// * `neighbor_id` - Neighbor ID of the message.
    /// * `priority` - Priority value.
    pub fn push_priority(&mut self, node_id: usize, neighbor_id: usize, priority: f64) {
        self.priorities.push((node_id as u32, neighbor_id as u32), OrderedFloat(priority));
    }

    /// Retrieves the highest priority message from the scheduling queue.
    ///
    /// # Returns
    /// Tuple of `(node_id, neighbor_id)` and priority.
    pub fn get_highest_priority(&self) -> Result<((usize, usize), f64), Box<dyn std::error::Error>> {
        let (&(end_id, start_in_end_id), &residual) = self.priorities.peek().ok_or("Priorities is empty")?;
        Ok(((end_id as usize, start_in_end_id as usize), residual.0))
    }

    /// Runs zero-lookahead belief propagation until convergence or limit.
    ///
    /// # Arguments
    /// * `max_loops` - Maximum number of iterations.
    /// * `tolerance` - Residual threshold for convergence.
    ///
    /// # Returns
    /// Final node beliefs as a list of value vectors.
    pub fn zero_lookahead_bp(&mut self, max_loops: u32, tolerance: f64) -> Result<Vec<Vec<f64>>, Box<dyn std::error::Error>> {

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
            let node = self.graph.get_node(node_id);
            for neighbor_id in 0..node.neighbors_count() {
                let residual = self.compute_infinity_norm_residual(node_id, neighbor_id);
                if residual > tolerance {
                    self.push_priority(node_id, neighbor_id, residual);
                }
            }
        }

        let mut k = 5;

        // keep track of the nodes of which the incoming messages have changed
        let mut prev_changed: Vec<usize> = (0..self.msg_in.len()).collect();
        while k < max_loops && max_residual > tolerance && ! self.priorities.is_empty() {

            // actual zero-look-ahead-BP part
            let ((end_id, start_in_end_id), max_residual) = self.get_highest_priority()?;
            
            let end_node = self.graph.get_node(end_id);
            let (start_id, end_in_start_id) = self.graph.get_neighbor_node_and_neighbor_id(end_node, start_in_end_id);
                        
            self.single_message_update(start_id, end_id, end_in_start_id, None)?;

            let priority_residual = self.compute_infinity_norm_residual(end_id, start_in_end_id);
            for node_id in prev_changed {
                let msg: MessagesInNode = self.msg_in[node_id].clone();
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
                        prev_changed.push(neighbor_id);
                        self.msg_in[neighbor_id] = self.msg_in_new[neighbor_id].clone();
                    }
                },
                _ => {
                    if self.msg_in[end_id].is_ct_message(start_in_end_id) {
                        self.msg_in[end_id].set_ct_message(self.msg_in_new[end_id].get_ct_message()?.clone())?;
                    }
                    self.msg_in[end_id].set_message(start_in_end_id, self.msg_in_new[end_id].get_message(start_in_end_id));
                    prev_changed.push(end_id);
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
                    let incoming_messages: &Vec<[f64; 2]> = self.msg_in[node.get_id()].get_messages();
                    
                    let initial_belief: [f64; 2] = get_initial_belief(node).variable_values().ok_or("Node should have PeptideBelief or TaxonBelief")?;

                    let sum_logs: [f64;2] = incoming_messages.iter()
                        .fold([0.0;2], |mut acc,  row| {acc[0] += row[0].ln(); acc[1] += row[1].ln(); acc});

                    // Compute final log-normalized message
                    let mut logged_variable_marginal: [f64; 2] = [initial_belief[0].ln() + sum_logs[0], initial_belief[1].ln() + sum_logs[1]];
                    log_normalize(&mut logged_variable_marginal);

                    self.current_beliefs[node.get_id()] = match node.get_subtype() {
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
        let mut checked_cts: HashSet<usize> = HashSet::new();
        
        if self.max_val.is_some() && local_loops {

            let (_, start_id) = self.max_val.unwrap(); // Safe unwrap because of is_some() check
            let start_node = self.graph.get_node(start_id);
            let neighbors: Vec<usize> = self.graph.get_neighbors(start_node).collect();
            for (end_in_start_id, &end_id) in neighbors.iter().enumerate() {
                self.single_message_update(start_id, end_id, end_in_start_id, Some(&mut checked_cts))?;
            }

        } else {
            for id in 0..self.graph.node_count() {
                let start_node = self.graph.get_node(id);
                let neighbors: Vec<usize> = self.graph.get_neighbors(start_node).collect();
                for (end_in_start_id, &end_id) in neighbors.iter().enumerate() {
                    self.single_message_update(id, end_id, end_in_start_id, Some(&mut checked_cts))?;
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
    fn single_message_update(&mut self, start_id: usize, end_id: usize, end_in_start_id: usize, mut checked_cts: Option<&mut HashSet<usize>>) -> Result<(), Box<dyn std::error::Error>> {

        let start_node: &Node = self.graph.get_node(start_id);
        let (_, start_in_end_id): (usize, usize) = self.graph.get_neighbor_node_and_neighbor_id(start_node, end_in_start_id);

        match start_node.get_subtype() {
            NodeType::TaxonNode { .. } | NodeType::PeptideNode { .. } => {
                let new_message = self.compute_out_message_variable(start_id, end_id, end_in_start_id);
                let end_node: &Node = self.graph.get_node(end_id);
                self.msg_in_new[end_id].set_message(start_in_end_id, new_message);
            },
            NodeType::FactorNode { .. } => {
                if self.graph.get_node(end_id).is_convolution_tree_node() {
                    let new_message = self.compute_out_message_factor_ctree(start_id, end_id, end_in_start_id)?;
                    self.msg_in_new[end_id].set_ct_message(new_message)?;
                } else {
                    let new_message = self.compute_out_message_factor(start_id, end_id, end_in_start_id)?;
                    self.msg_in_new[end_id].set_message(start_in_end_id, new_message);
                }
            },
            NodeType::ConvolutionTreeNode { number_of_parents } => 
                if checked_cts.as_ref().map_or(true, |set| ! set.contains(&start_id)) {
                    self.compute_out_messages_ct_tree(start_id, *number_of_parents as usize)?;
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
    /// Normalized probability array.
    fn compute_out_message_variable(&mut self, start_id: usize, end_id: usize, end_in_start_id: usize) -> [f64; 2] {
        
        let start_node = self.graph.get_node(start_id);
        let end_node = self.graph.get_node(end_id);

        let node_belief: [f64; 2] = self.current_beliefs[start_id].variable_values().expect("Start node belief should be a TaxonBelief or PeptideBelief");

        if self.msg_in[start_id].get_message_count() <= 1 {
            // TODO: There was an if case that was always true, because current_beliefs only updated after 
            return node_belief;
        }

        // Sum of incoming messages, need for logs to prevent underflow in very large multiplications
        let sum_logs: [f64; 2] = sum_logs_batched(self.msg_in[start_id].get_messages());

        // Compute final log-normalized message, Log transform node_belief
        let mut out_message_log: [f64; 2] = [node_belief[0].ln() + sum_logs[0], node_belief[1].ln() + sum_logs[1]];
        log_normalize(&mut out_message_log);

        // Prevent underflow: Replace zeros with 1e-30
        avoid_underflow_arr(&mut out_message_log);

        out_message_log
    }

    /// Computes outgoing message for factor nodes, for edges not going to a ct tree node.
    ///
    /// # Arguments
    /// * `start_id` - ID of source node.
    /// * `end_id` - ID of destination node.
    /// * `end_in_start_id` - Neighbor index of destination in source.
    /// 
    /// # Returns
    /// Normalized probability vector.
    fn compute_out_message_factor(&mut self, start_id: usize, end_id: usize, end_in_start_id: usize) -> Result<[f64; 2], Box<dyn std::error::Error>> {
        let incoming_messages: &MessagesInNode = &self.msg_in[start_id];
        let node_belief: &Vec<[f64;2]> = self.current_beliefs[start_id].factor_values().ok_or("factor_values called on a NodeBelief which is not a FactorBelief")?;
        
        if incoming_messages.is_ct_message(0) {
            let mut out_message: Vec<f64> = node_belief.iter().zip(incoming_messages.get_ct_message()?.iter())
                .map(|(&a, &b)| a[0] * b + a[1] * b).collect();
            normalize(&mut out_message);

            let mut out_message: [f64; 2] = [out_message[0], out_message[1]];
            // Prevent underflow: Replace zeros with 1e-30
            avoid_underflow_arr(&mut out_message);

            return Ok(out_message);
        } else {
            // Product of incoming messages
            // incoming_messages is always a 2x2, we must ignore row end_in_start_id, so product is just the other row
            let prod: [f64; 2] = incoming_messages.get_message(1-end_in_start_id);
            
            // Compute final normalized message
            let mut out_message: Vec<[f64;2]> = node_belief.iter().map(|&a| [a[0] * prod[0], a[1] * prod[1]]).collect();

            normalize_2d(&mut out_message);
            let out_message_sum: Vec<f64> = out_message.iter().map(|[a, b]| a + b).collect();

            return Ok([out_message_sum[0], out_message_sum[1]]);
        }
    }

    /// Computes outgoing message for factor nodes, for edges going to a ct tree node.
    ///
    /// # Arguments
    /// * `start_id` - ID of source node.
    /// * `end_id` - ID of destination node.
    /// * `end_in_start_id` - Neighbor index of destination in source.
    /// 
    /// # Returns
    /// Normalized probability vector.
    fn compute_out_message_factor_ctree(&mut self, start_id: usize, end_id: usize, end_in_start_id: usize) -> Result<Vec<f64>, Box<dyn std::error::Error>> {
        let incoming_messages: &MessagesInNode = &self.msg_in[start_id];
        let node_belief: &Vec<[f64;2]> = self.current_beliefs[start_id].factor_values().ok_or("factor_values called on a NodeBelief which is not a FactorBelief")?;

        // Product of incoming messages
        // Dimensions of incoming messages is always 1x2, so product is not necessary.
        let prod: [f64;2] = incoming_messages.get_message(1);
        
        // Compute final normalized message
        let mut out_message: Vec<[f64;2]> = node_belief.iter().map(|&a| [a[0] * prod[0], a[1] * prod[1]]).collect();
        normalize_2d(&mut out_message);
        let out_message_sum: Vec<f64> = out_message.iter().map(|[a, b]| a + b).collect();
        
        return Ok(out_message_sum);
    }

    /// Computes outgoing messages for convolution tree nodes.
    ///
    /// # Arguments
    /// * `start_id` - Node ID of convolution tree.
    /// * `number_of_parents` - Number of parent nodes.
    fn compute_out_messages_ct_tree(&mut self, start_id: usize, number_of_parents: usize) -> Result<(), Box<dyn std::error::Error>> {
        let start_node = self.graph.get_node(start_id);

        let neighbor_list: Vec<usize> = self.graph.get_neighbors(start_node).collect();
        let prot_list: &[usize] = &neighbor_list[1..];
        let (factor_id, start_id_in_factor) = self.graph.get_neighbor_node_and_neighbor_id(start_node, 0);
        
        let mut shared_likelihoods: Vec<f64> = self.msg_in[start_id].get_ct_message()?.clone();
        let mut old_shared_likelihoods: Vec<f64> = shared_likelihoods.iter()
            .zip(self.msg_in_log[start_id].get_ct_message()?)
            .map(|(a, b)| a * b).collect();

        let prot_prob_list: &Vec<[f64; 2]> = self.msg_in[start_id].get_messages();

        let old_prot_prob_list: &Vec<[f64; 2]> = self.msg_in_log[start_id].get_messages();

        if old_shared_likelihoods != *shared_likelihoods && 
            prot_prob_list.iter().zip(old_prot_prob_list.iter()).any(|(a, b)| a[0] != b[0]) {
            let convolution_tree = ConvolutionTree::new(shared_likelihoods, prot_prob_list.clone())?;

            for (protein_id, protein) in prot_list.iter().enumerate() {
                let (_, node_neighbor_index): (usize, usize) = self.graph.get_neighbor_node_and_neighbor_id(start_node, protein_id+1);
                let new_message = convolution_tree.message_to_variable(protein_id);
                let mut new_message = [new_message[0], new_message[1]];
                avoid_underflow_arr(&mut new_message);
                self.msg_in_new[*protein].set_message(node_neighbor_index, new_message);
            }

            let mut new_message = convolution_tree.message_to_shared_likelihood()?;

            avoid_underflow(&mut new_message);
            self.msg_in_new[factor_id].set_ct_message(new_message)?;

        } else {
            
            for (protein_id, &protein) in prot_list.iter().enumerate() {
                let (_, node_neighbor_index): (usize, usize) = self.graph.get_neighbor_node_and_neighbor_id(start_node, protein_id+1);
                self.msg_in_new[protein].set_message(node_neighbor_index, self.msg_in[protein].get_message(node_neighbor_index).clone());
            }

            let ct_msg = self.msg_in[factor_id].get_ct_message()?.clone();
            self.msg_in_new[factor_id].set_ct_message(ct_msg)?;
        }

        Ok(())
    }

    /// Computes infinity norm residual for an edge.
    ///
    /// # Arguments
    /// * `end_id` - ID of destination node.
    /// * `start_in_end_id` - Neighbor index of destination in source.
    /// 
    /// # Returns
    /// Residual as `f64`.
    fn compute_infinity_norm_residual(&mut self, end_id: usize, start_in_end_id: usize) -> f64 {

        if self.msg_in[end_id].is_ct_message(start_in_end_id) {
            let msg1 = self.msg_in[end_id].get_ct_message().unwrap();
            let msg2 = self.msg_in_log[end_id].get_ct_message().unwrap();

            let residual = msg1.iter()
                            .zip(msg2.iter())
                            .map(|(m1, m2)| (m1 / m2).ln().abs())
                            .fold(f64::NEG_INFINITY, f64::max);

            return residual;
        }

        let msg1: &mut [f64; 2] = &mut self.msg_in[end_id].get_message(start_in_end_id);
        let msg2: &[f64; 2] = &self.msg_in_log[end_id].get_message(start_in_end_id);

        (msg1[0] / msg2[0]).max(msg1[1] / msg2[1]).ln().abs()
    }

    /// Updates total residuals after a message update.
    /// 
    /// # Arguments
    /// * `start_id` - ID of source node.
    /// * `end_id` - ID of destination node.
    /// * `end_in_start_id` - Neighbor index of destination in source.
    /// * `current_residual` - Redidual to add to message
    fn compute_total_residuals(&mut self, start_id: usize, end_id: usize, start_in_end_id: usize, current_residual: f64) {
        let start_node = self.graph.get_node(start_id);
        let end_node = self.graph.get_node(end_id);
        let (_, end_in_start_id) = self.graph.get_neighbor_node_and_neighbor_id(end_node, start_in_end_id);

        for (i, neighbor_id) in self.graph.get_neighbors(start_node).enumerate() {
            if neighbor_id != end_id {
                self.total_residuals[start_id][i][end_in_start_id] = 0.0;
            }
        }

        for (i, neighbor_id) in self.graph.get_neighbors(end_node).enumerate() {
            if neighbor_id != start_id {
                self.total_residuals[end_id][start_in_end_id][i] += current_residual;
            }
        }
    }

    /// Updates message priorities for scheduling.
    /// 
    /// # Arguments
    /// * `start_id` - ID of source node.
    /// * `end_id` - ID of destination node.
    /// * `end_in_start_id` - Neighbor index of destination in source.
    fn compute_priority(&mut self, start_id: usize, end_id: usize, start_in_end_id: usize) -> Result<(), Box<dyn std::error::Error>> {
        let end_node = self.graph.get_node(end_id);

        self.priorities.remove(&(end_id as u32, start_in_end_id as u32));

        for i in 0..end_node.neighbors_count() {
            let (neighbor_id, end_in_neighbor_id) = self.graph.get_neighbor_node_and_neighbor_id(&end_node, i);
            if neighbor_id != start_id {
                let priority: f64 = self.graph
                    .get_neighbors(end_node)
                    .enumerate()
                    .map(|(j, sum_run)| {
                        if sum_run != neighbor_id { 
                            self.total_residuals[end_id][j][i]
                        } else { 
                            0.0
                        }
                    }).sum();

                if self.priorities.change_priority(&(neighbor_id as u32, end_in_neighbor_id as u32), OrderedFloat(priority)).is_none() {
                    self.priorities.push((neighbor_id as u32, end_in_neighbor_id as u32), OrderedFloat(priority));
                }
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

        let edge0 = Edge::new(0, 0, 1, 0, 0, None);
        nodes[0].add_incident_edge(0);
        nodes[1].add_incident_edge(0);
        edges.push(edge0);


        let edge1 = Edge::new(1, 1, 2, 0, 1, None);
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
        assert_eq!(fb.factor_values(), Some(&vec![[0.2,0.8]]));
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
        assert!(beliefs.is_ok());
        let beliefs = beliefs.unwrap();

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

        let s: f64 = msg[0] + msg[1];
        assert!((s-1.0).abs()<1e-6);
    }

    #[test]
    fn test_compute_out_message_factor() {
        let graph = create_minimal_graph();
        let mut messages = Messages::new(graph);
        let msg = messages.compute_out_message_factor(1,2,0);
        assert!(msg.is_ok());
        let msg = msg.unwrap();

        let s: f64 = msg[0] + msg[1];
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
                messages.priorities.push((node_id, neighbor_id), OrderedFloat(0.1));
            }
        }
        messages.compute_priority(0,2,0);
        assert!(messages.priorities.peek().is_some());
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
