use rustfft::{FftPlanner, num_complex::Complex, num_traits::Zero};
use crate::array_utils::normalize;

#[derive(Debug, Clone)]
struct CTNode {
    joint_above: Vec<f64>,
    likelihood_below: Option<Vec<f64>>,
}

impl CTNode {
    /// Creates a new CTNode with a given probability vector.
    /// 
    /// Normalizes the input vector so that it sums to 1.
    /// 
    /// # Arguments
    /// * `joint_above` - A vector of probabilities for this node.
    fn new(mut joint_above: Vec<f64>) -> Self {
        normalize(&mut joint_above);
        Self {
            joint_above: joint_above,
            likelihood_below: None,
        }
    }

    /// Creates a count node by convolving the joint distributions of two child nodes.
    /// 
    /// # Arguments
    /// * `lhs` - Left child node.
    /// * `rhs` - Right child node.
    /// 
    /// # Returns
    /// A new CTNode representing the combined distribution.
    fn create_count_node(lhs: CTNode, rhs: CTNode) -> CTNode {
        let joint_above = fft_convolve(&lhs.joint_above, &rhs.joint_above);
        let node = CTNode::new(joint_above);
        node
    }

    /// Computes the upward message by convolving the likelihood below with a sibling's distribution.
    /// 
    /// # Arguments
    /// * `answer_size` - Size of the output message.
    /// * `other_joint_vector` - The joint distribution of the sibling node.
    /// 
    /// # Returns
    /// A normalized probability vector representing the upward message.
    fn message_up(&self, answer_size: usize, other_joint_vector: &Vec<f64>) -> Vec<f64> {
        let likelihood = self.likelihood_below.as_ref().expect("Likelihood below is None!");
        let starting_point = other_joint_vector.len() - 1;
        let result = fft_convolve(
            &other_joint_vector.iter().rev().cloned().collect::<Vec<f64>>(),
            likelihood,
        );
        let mut result = result[starting_point..starting_point + answer_size].to_vec();
        normalize(&mut result);
        result
    }

    /// Returns the likelihood messages from below for this node.
    /// 
    /// # Returns
    /// A vector of probabilities from the likelihood below.
    fn messages_up(&self) -> Vec<f64> {
        self.likelihood_below.clone().expect("Likelihood below is None!")
    }
}

#[derive(Debug)]
pub struct ConvolutionTree {
    n_to_shared_likelihoods: Vec<f64>,
    log_length: usize,
    all_layers: Vec<Vec<CTNode>>,
    protein_layer: Vec<CTNode>,
    n_proteins: usize
}

impl ConvolutionTree {
    /// Constructs a ConvolutionTree given shared likelihoods and protein probability vectors.
    /// 
    /// # Arguments
    /// * `n_to_shared_likelihoods` - Shared likelihoods vector.
    /// * `proteins` - A vector of protein probability distributions.
    /// 
    /// # Returns
    /// A fully constructed ConvolutionTree with messages propagated backward.
    pub fn new(n_to_shared_likelihoods: Vec<f64>, proteins: Vec<Vec<f64>>) -> Self {
        let log_length = (proteins.len() as f64).log2().ceil() as usize;
        let mut tree = ConvolutionTree {
            n_to_shared_likelihoods,
            log_length,
            all_layers: Vec::new(),
            protein_layer: Vec::new(),
            n_proteins: proteins.len()
        };

        tree.build_first_layer(proteins);
        tree.build_remaining_layers();
        tree.propagate_backward();

        tree
    }

    /// Builds the first layer of the tree from protein probability vectors.
    /// Pads with dummy nodes to make the layer length a power of 2.
    /// 
    /// # Arguments
    /// * `proteins` - A vector of protein probability distributions.
    fn build_first_layer(&mut self, proteins: Vec<Vec<f64>>) {
        let mut layer = proteins.into_iter().map(CTNode::new).collect::<Vec<CTNode>>();

        // Pad with dummy nodes to make the length a power of 2
        while layer.len() < 2usize.pow(self.log_length as u32) {
            layer.push(CTNode::new(vec![1.0, 0.0]));
        }

        self.all_layers.push(layer);
    }

    /// Builds all remaining layers of the tree by combining nodes into count nodes.
    /// Each layer is constructed by convolving pairs of nodes from the previous layer.
    fn build_remaining_layers(&mut self) {
        for _ in 0..self.log_length {
            let most_recent_layer = self.all_layers.last().unwrap();
            let mut new_layer = Vec::new();

            for i in (0..most_recent_layer.len()).step_by(2) {
                let left = most_recent_layer[i].clone();
                let right = most_recent_layer[i + 1].clone();
                new_layer.push(CTNode::create_count_node(left, right));
            }

            self.all_layers.push(new_layer);
        }

        let mut likelihood_below = self.n_to_shared_likelihoods.clone();
        normalize(&mut likelihood_below);
        self.all_layers.last_mut().unwrap()[0].likelihood_below = Some(likelihood_below);
    }

    /// Propagates likelihood messages from the root down to the protein nodes.
    /// Updates each node's `likelihood_below` using upward messages from parent nodes.
    fn propagate_backward(&mut self) {
        for l in (1..=self.log_length).rev() {
            for i in 0..self.all_layers[l].len() {
                
                let left_parent = &self.all_layers[l-1][2*i];
                let right_parent = &self.all_layers[l-1][2*i + 1];
                let node = &self.all_layers[l][i];

                let likelihood_below_left = Some(node.message_up(left_parent.joint_above.len(), &right_parent.joint_above));
                let likelihood_below_right = Some(node.message_up(right_parent.joint_above.len(), &left_parent.joint_above));

                self.all_layers[l-1][2*i].likelihood_below = likelihood_below_left;
                self.all_layers[l-1][2*i+1].likelihood_below = likelihood_below_right;
            }
        }

        self.protein_layer = self.all_layers[0].clone();
    }

    /// Retrieves the message to a variable node (protein).
    /// 
    /// # Arguments
    /// * `prot_idx` - Index of the protein node.
    /// 
    /// # Returns
    /// A probability vector representing the message to that protein.
    pub fn message_to_variable(&self, prot_idx: usize) -> Vec<f64> {
        self.protein_layer[prot_idx].messages_up()
    }

    /// Retrieves the message to the shared likelihood node.
    /// 
    /// # Returns
    /// A probability vector representing the message to the shared likelihood.
    pub fn message_to_shared_likelihood(&self) -> Vec<f64> {

        // Extract the required range
        self.all_layers.last().unwrap()[0].joint_above[..=self.n_proteins].to_vec()
    }
}


/// Performs FFT-based convolution of two probability vectors.
/// 
/// # Arguments
/// * `a` - First probability vector.
/// * `b` - Second probability vector.
/// 
/// # Returns
/// A new vector representing the convolution of `a` and `b`.
fn fft_convolve(a: &Vec<f64>, b: &Vec<f64>) -> Vec<f64> {
    let len = a.len() + b.len() - 1;
    let fft_size = len.next_power_of_two();

    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(fft_size);
    let ifft = planner.plan_fft_inverse(fft_size);

    let mut a_complex: Vec<Complex<f64>> = a.iter().map(|&x| Complex::new(x, 0.0)).collect();
    let mut b_complex: Vec<Complex<f64>> = b.iter().map(|&x| Complex::new(x, 0.0)).collect();
    a_complex.resize(fft_size, Complex::zero());
    b_complex.resize(fft_size, Complex::zero());

    fft.process(&mut a_complex);
    fft.process(&mut b_complex);

    let mut result_complex: Vec<Complex<f64>> = a_complex.iter().zip(&b_complex).map(|(x, y)| x * y).collect();
    ifft.process(&mut result_complex);

    result_complex.iter().take(len).map(|c| c.re / fft_size as f64).collect()
}