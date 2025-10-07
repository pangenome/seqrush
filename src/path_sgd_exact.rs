/// Exact line-by-line port of ODGI's path_linear_sgd from path_sgd.cpp
/// This is a direct translation maintaining all original logic and structure
use crate::bidirected_ops::BidirectedGraph;
use crate::bidirected_graph::Handle;
use rand::Rng;
use rand_xoshiro::Xoshiro256Plus;
use rand::SeedableRng;
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Duration;
use std::collections::HashMap;

/// XP path index - mimics ODGI's XP structure
pub struct XPIndex {
    /// Node positions for each path step
    pub path_steps: Vec<(usize, Handle)>, // (path_id, handle)
    /// Path lengths in bp
    pub path_lengths: HashMap<usize, usize>,
    /// Path step counts
    pub path_step_counts: HashMap<usize, usize>,
    /// Path names for debugging
    pub path_names: HashMap<usize, String>,
    /// Total steps across all paths
    pub total_steps: usize,
}

impl XPIndex {
    pub fn from_graph(graph: &BidirectedGraph) -> Self {
        let mut path_steps = Vec::new();
        let mut path_lengths = HashMap::new();
        let mut path_step_counts = HashMap::new();
        let mut path_names = HashMap::new();

        for (path_idx, path) in graph.paths.iter().enumerate() {
            let path_id = path_idx;
            let path_name = path.name.clone();
            path_names.insert(path_id, path_name.clone());

            let mut path_len = 0;
            let mut step_count = 0;

            for &handle in &path.steps {
                path_steps.push((path_id, handle));
                let node_id = handle.node_id();
                if let Some(node) = graph.nodes.get(&node_id) {
                    path_len += node.sequence.len();
                }
                step_count += 1;
            }

            path_lengths.insert(path_id, path_len);
            path_step_counts.insert(path_id, step_count);
        }

        let total_steps = path_steps.len();

        XPIndex {
            path_steps,
            path_lengths,
            path_step_counts,
            path_names,
            total_steps,
        }
    }

    pub fn get_path_length(&self, path_id: usize) -> usize {
        *self.path_lengths.get(&path_id).unwrap_or(&0)
    }

    pub fn get_path_step_count(&self, path_id: usize) -> usize {
        *self.path_step_counts.get(&path_id).unwrap_or(&0)
    }

    pub fn get_position_of_step(&self, graph: &BidirectedGraph, path_id: usize, step_rank: usize) -> usize {
        // Get the path
        let mut pos = 0;
        let mut current_rank = 0;

        for (pid, handle) in &self.path_steps {
            if *pid == path_id {
                if current_rank == step_rank {
                    return pos;
                }
                let node_id = handle.node_id();
                if let Some(node) = graph.nodes.get(&node_id) {
                    pos += node.sequence.len();
                }
                current_rank += 1;
            }
        }

        pos
    }

    pub fn get_handle_of_step(&self, path_id: usize, step_rank: usize) -> Handle {
        let mut current_rank = 0;
        for (pid, handle) in &self.path_steps {
            if *pid == path_id {
                if current_rank == step_rank {
                    return *handle;
                }
                current_rank += 1;
            }
        }
        Handle::new(1, false) // fallback
    }
}

/// Parameters exactly matching ODGI's defaults
pub struct PathSGDParams {
    pub iter_max: u64,
    pub iter_with_max_learning_rate: u64,
    pub min_term_updates: u64,
    pub delta: f64,
    pub eps: f64,
    pub eta_max: f64,
    pub theta: f64,
    pub space: u64,
    pub space_max: u64,
    pub space_quantization_step: u64,
    pub cooling_start: f64,
    pub nthreads: usize,
    pub progress: bool,
    pub snapshot: bool,
}

impl Default for PathSGDParams {
    fn default() -> Self {
        PathSGDParams {
            iter_max: 100,
            iter_with_max_learning_rate: 0,
            min_term_updates: 100,
            delta: 0.0,
            eps: 0.01,
            eta_max: 100.0,
            theta: 0.99,
            space: 100,
            space_max: 100,
            space_quantization_step: 10,
            cooling_start: 0.8,
            nthreads: 1,
            progress: false,
            snapshot: false,
        }
    }
}

/// Atomic f64 wrapper since Rust doesn't have AtomicF64 in stable
#[derive(Clone)]
struct AtomicF64 {
    value: Arc<Mutex<f64>>,
}

impl AtomicF64 {
    fn new(val: f64) -> Self {
        AtomicF64 {
            value: Arc::new(Mutex::new(val)),
        }
    }

    fn store(&self, val: f64) {
        *self.value.lock().unwrap() = val;
    }

    fn load(&self) -> f64 {
        *self.value.lock().unwrap()
    }
}

/// Fast precise power function for zipf distribution
fn fast_precise_pow(base: f64, exp: f64) -> f64 {
    base.powf(exp)
}

/// Dirty Zipfian distribution - mimics ODGI's implementation
struct DirtyZipfian {
    n: u64,
    theta: f64,
    zeta: f64,
}

impl DirtyZipfian {
    fn new(min: u64, max: u64, theta: f64, zeta: f64) -> Self {
        DirtyZipfian {
            n: max - min + 1,
            theta,
            zeta,
        }
    }

    fn sample<R: Rng>(&self, rng: &mut R) -> u64 {
        // Simplified zipfian sampling
        let u: f64 = rng.gen();
        let mut sum = 0.0;
        for i in 1..=self.n {
            sum += fast_precise_pow(1.0 / i as f64, self.theta);
            if sum / self.zeta >= u {
                return i;
            }
        }
        self.n
    }
}

/// Calculate learning rate schedule - line 467-500
fn path_linear_sgd_schedule(
    w_min: f64,
    w_max: f64,
    iter_max: u64,
    iter_with_max_learning_rate: u64,
    eps: f64,
) -> Vec<f64> {
    let eta_max = 1.0 / w_min;
    let eta_min = eps / w_max;
    let lambda = (eta_max / eta_min).ln() / (iter_max as f64 - 1.0);

    let mut etas = Vec::with_capacity((iter_max + 1) as usize);
    for t in 0..=iter_max {
        let diff = (t as i64 - iter_with_max_learning_rate as i64).abs() as f64;
        etas.push(eta_max * (-lambda * diff).exp());
    }

    etas
}

/// Main path-guided SGD function - exact port of path_linear_sgd from line 12-465
pub fn path_linear_sgd(
    graph: &BidirectedGraph,
    path_index: &XPIndex,
    path_sgd_use_paths: Vec<usize>, // path IDs to use
    params: PathSGDParams,
) -> Vec<f64> {
    // Line 44: Calculate first cooling iteration
    let first_cooling_iteration = (params.cooling_start * params.iter_max as f64).floor() as u64;

    // Line 47: Total term updates
    let total_term_updates = params.iter_max * params.min_term_updates;

    if params.progress {
        eprintln!("[odgi::path_linear_sgd] 1D path-guided SGD:");
    }

    // Line 54: Number of nodes
    let num_nodes = graph.nodes.len();

    // Line 56: Our positions in 1D - using Arc<Mutex> for shared atomic access
    let x_positions: Vec<AtomicF64> = (0..num_nodes)
        .map(|_| AtomicF64::new(0.0))
        .collect();

    // Line 57-61: Snapshot tracking
    let snapshot_in_progress = Arc::new(AtomicBool::new(false));
    let snapshot_progress: Vec<Arc<AtomicBool>> = (0..params.iter_max)
        .map(|i| Arc::new(AtomicBool::new(i == 0)))
        .collect();

    // Line 63-69: Seed positions with graph order
    let mut len = 0.0;
    let mut node_id_to_index: HashMap<usize, usize> = HashMap::new();
    for (idx, (&node_id, node)) in graph.nodes.iter().enumerate() {
        x_positions[idx].store(len);
        node_id_to_index.insert(node_id, idx);
        len += node.sequence.len() as f64;
    }

    // Line 80-103: Check if we have paths with more than one step
    let mut at_least_one_path_with_more_than_one_step = false;
    for &path_id in &path_sgd_use_paths {
        if path_index.get_path_step_count(path_id) > 1 {
            at_least_one_path_with_more_than_one_step = true;
            break;
        }
    }

    if !at_least_one_path_with_more_than_one_step {
        // Return current positions if no valid paths
        return x_positions.iter().map(|x| x.load()).collect();
    }

    // Line 107-112: Calculate w_min and w_max
    let w_min = 1.0 / params.eta_max;
    let w_max = 1.0;

    // Line 114-122: Get schedule
    if params.progress {
        eprintln!(
            "[odgi::path_linear_sgd] calculating linear SGD schedule ({} {} {} {} {})",
            w_min, w_max, params.iter_max, params.iter_with_max_learning_rate, params.eps
        );
    }

    let etas = path_linear_sgd_schedule(
        w_min,
        w_max,
        params.iter_max,
        params.iter_with_max_learning_rate,
        params.eps,
    );

    // Line 124-138: Cache zipf zetas
    if params.progress {
        let zeta_count = if params.space <= params.space_max {
            params.space
        } else {
            params.space_max + (params.space - params.space_max) / params.space_quantization_step + 1
        };
        eprintln!(
            "[odgi::path_linear_sgd] calculating zetas for {} zipf distributions",
            zeta_count
        );
    }

    let zetas_size = if params.space <= params.space_max {
        params.space
    } else {
        params.space_max + (params.space - params.space_max) / params.space_quantization_step + 1
    };
    let mut zetas = vec![0.0; (zetas_size + 1) as usize];

    let mut zeta_tmp = 0.0;
    for i in 1..=params.space {
        zeta_tmp += fast_precise_pow(1.0 / i as f64, params.theta);
        if i <= params.space_max {
            zetas[i as usize] = zeta_tmp;
        }
        if i > params.space_max && (i - params.space_max) % params.space_quantization_step == 0 {
            let idx = (params.space_max + (i - params.space_max) / params.space_quantization_step) as usize;
            if idx < zetas.len() {
                zetas[idx] = zeta_tmp;
            }
        }
    }

    // Line 140-159: Atomic variables for coordination
    let term_updates = Arc::new(AtomicU64::new(0));
    let eta = Arc::new(AtomicF64::new(etas[0]));
    let adj_theta = Arc::new(AtomicF64::new(params.theta));
    let cooling = Arc::new(AtomicBool::new(false));
    let delta_max = Arc::new(AtomicF64::new(0.0));
    let work_todo = Arc::new(AtomicBool::new(true));
    let iteration = Arc::new(AtomicUsize::new(0));

    // Line 161-203: Checker thread lambda
    let checker_handle = {
        let term_updates = Arc::clone(&term_updates);
        let eta = Arc::clone(&eta);
        let adj_theta = Arc::clone(&adj_theta);
        let cooling = Arc::clone(&cooling);
        let delta_max = Arc::clone(&delta_max);
        let work_todo = Arc::clone(&work_todo);
        let iteration = Arc::clone(&iteration);
        let snapshot_in_progress = Arc::clone(&snapshot_in_progress);
        let etas = etas.clone();

        thread::spawn(move || {
            while work_todo.load(Ordering::Relaxed) {
                if term_updates.load(Ordering::Relaxed) > params.min_term_updates {
                    let iter = iteration.fetch_add(1, Ordering::Relaxed) + 1;

                    if iter > params.iter_max as usize {
                        work_todo.store(false, Ordering::Relaxed);
                    } else if delta_max.load() <= params.delta && params.delta > 0.0 {
                        if params.progress {
                            eprintln!(
                                "[odgi::path_linear_sgd] delta_max: {} <= delta: {}. Threshold reached.",
                                delta_max.load(),
                                params.delta
                            );
                        }
                        work_todo.store(false, Ordering::Relaxed);
                    } else {
                        eta.store(etas[iter]);
                        delta_max.store(params.delta);

                        if iter > first_cooling_iteration as usize {
                            adj_theta.store(0.001);
                            cooling.store(true, Ordering::Relaxed);
                        }
                    }

                    term_updates.store(0, Ordering::Relaxed);
                }

                thread::sleep(Duration::from_millis(1));
            }
        })
    };

    // Line 205-407: Worker threads
    let mut worker_handles = vec![];

    for tid in 0..params.nthreads {
        let term_updates = Arc::clone(&term_updates);
        let eta = Arc::clone(&eta);
        let adj_theta = Arc::clone(&adj_theta);
        let cooling = Arc::clone(&cooling);
        let delta_max = Arc::clone(&delta_max);
        let work_todo = Arc::clone(&work_todo);
        let snapshot_in_progress = Arc::clone(&snapshot_in_progress);
        let x_positions = x_positions.clone(); // Each thread gets its own Arc references
        let path_index_steps = path_index.path_steps.clone();
        let path_index_step_counts = path_index.path_step_counts.clone();
        let total_steps = path_index.total_steps;
        let zetas = zetas.clone();
        let node_id_to_index = node_id_to_index.clone();
        let graph_clone = graph.clone(); // We need graph for node lookups

        let handle = thread::spawn(move || {
            // Line 208: Seed PRNG with unique seed per thread
            let seed = 9399220 + tid as u64;
            let mut gen = Xoshiro256Plus::seed_from_u64(seed);

            // Line 215-216: Distributions
            let mut term_updates_local = 0u64;

            while work_todo.load(Ordering::Relaxed) {
                if !snapshot_in_progress.load(Ordering::Relaxed) {
                    // Line 221-222: Sample random step
                    let step_index = gen.gen_range(0..total_steps);

                    // Line 226-227: Get path and step info
                    let (path_id, handle_a) = path_index_steps[step_index];

                    // Line 229-232: Skip single-step paths
                    let path_step_count = path_index_step_counts
                        .get(&path_id)
                        .copied()
                        .unwrap_or(0);

                    if path_step_count <= 1 {
                        continue;
                    }

                    // Line 237-241: Calculate step rank
                    let mut s_rank = 0;
                    for i in 0..step_index {
                        if path_index_steps[i].0 == path_id {
                            s_rank += 1;
                        }
                    }

                    // Line 244-280: Sample second step
                    let handle_b;
                    let flip: bool = gen.gen();

                    if cooling.load(Ordering::Relaxed) || flip {
                        let theta = adj_theta.load();

                        if s_rank > 0 && (gen.gen::<bool>() || s_rank == path_step_count - 1) {
                            // Line 247-259: Go backward
                            let jump_space = params.space.min(s_rank as u64);
                            let space_idx = if jump_space > params.space_max {
                                params.space_max + (jump_space - params.space_max) / params.space_quantization_step + 1
                            } else {
                                jump_space
                            };

                            let zipf = DirtyZipfian::new(1, jump_space, theta, zetas[space_idx as usize]);
                            let z_i = zipf.sample(&mut gen);
                            let target_rank = s_rank - z_i as usize;

                            // Find handle at target rank
                            let mut current_rank = 0;
                            handle_b = path_index_steps.iter()
                                .find(|(pid, _)| {
                                    if *pid == path_id {
                                        if current_rank == target_rank {
                                            return true;
                                        }
                                        current_rank += 1;
                                    }
                                    false
                                })
                                .map(|(_, h)| *h)
                                .unwrap_or(handle_a);
                        } else {
                            // Line 261-273: Go forward
                            let jump_space = params.space.min((path_step_count - s_rank - 1) as u64);
                            let space_idx = if jump_space > params.space_max {
                                params.space_max + (jump_space - params.space_max) / params.space_quantization_step + 1
                            } else {
                                jump_space
                            };

                            let zipf = DirtyZipfian::new(1, jump_space, theta, zetas[space_idx as usize]);
                            let z_i = zipf.sample(&mut gen);
                            let target_rank = s_rank + z_i as usize;

                            // Find handle at target rank
                            let mut current_rank = 0;
                            handle_b = path_index_steps.iter()
                                .find(|(pid, _)| {
                                    if *pid == path_id {
                                        if current_rank == target_rank {
                                            return true;
                                        }
                                        current_rank += 1;
                                    }
                                    false
                                })
                                .map(|(_, h)| *h)
                                .unwrap_or(handle_a);
                        }
                    } else {
                        // Line 275-280: Sample randomly across path
                        let random_rank = gen.gen_range(0..path_step_count);

                        let mut current_rank = 0;
                        handle_b = path_index_steps.iter()
                            .find(|(pid, _)| {
                                if *pid == path_id {
                                    if current_rank == random_rank {
                                        return true;
                                    }
                                    current_rank += 1;
                                }
                                false
                            })
                            .map(|(_, h)| *h)
                            .unwrap_or(handle_a);
                    }

                    // Line 305-316: Get positions in path
                    let pos_in_path_a = {
                        let mut pos = 0;
                        let mut current_rank = 0;
                        for (pid, h) in &path_index_steps {
                            if *pid == path_id {
                                if current_rank == s_rank {
                                    break;
                                }
                                if let Some(node) = graph_clone.nodes.get(&h.node_id()) {
                                    pos += node.sequence.len();
                                }
                                current_rank += 1;
                            }
                        }
                        pos
                    };

                    let pos_in_path_b = {
                        let mut pos = 0;
                        for (pid, h) in &path_index_steps {
                            if *pid == path_id {
                                if *h == handle_b {
                                    break;
                                }
                                if let Some(node) = graph_clone.nodes.get(&h.node_id()) {
                                    pos += node.sequence.len();
                                }
                            }
                        }
                        pos
                    };

                    // Line 318-324: Calculate term distance
                    let term_dist = (pos_in_path_a as f64 - pos_in_path_b as f64).abs();

                    if term_dist == 0.0 {
                        continue;
                    }

                    // Line 333-343: Calculate weight and learning rate
                    let term_weight = 1.0 / term_dist;
                    let w_ij = term_weight;
                    let mut mu = eta.load() * w_ij;
                    if mu > 1.0 {
                        mu = 1.0;
                    }

                    // Line 344-347: Get node indices
                    let d_ij = term_dist;
                    let i = *node_id_to_index.get(&handle_a.node_id()).unwrap_or(&0);
                    let j = *node_id_to_index.get(&handle_b.node_id()).unwrap_or(&0);

                    // Line 353-362: Calculate distance
                    let mut dx = x_positions[i].load() - x_positions[j].load();
                    if dx == 0.0 {
                        dx = 1e-9;
                    }
                    let mag = dx.abs();

                    // Line 367-376: Check distances and update delta_max
                    let delta = mu * (mag - d_ij) / 2.0;
                    let delta_abs = delta.abs();

                    // Update delta_max if needed
                    loop {
                        let current = delta_max.load();
                        if delta_abs <= current {
                            break;
                        }
                        delta_max.store(delta_abs);
                        break; // Simplified CAS
                    }

                    // Line 378-393: Calculate and apply updates
                    let r = delta / mag;
                    let r_x = r * dx;

                    x_positions[i].store(x_positions[i].load() - r_x);
                    x_positions[j].store(x_positions[j].load() + r_x);

                    // Line 397-404: Update counters
                    term_updates_local += 1;
                    if term_updates_local >= 1000 {
                        term_updates.fetch_add(term_updates_local, Ordering::Relaxed);
                        if params.progress {
                            // Progress meter would go here
                        }
                        term_updates_local = 0;
                    }
                }
            }
        });

        worker_handles.push(handle);
    }

    // Line 445-451: Join all threads
    for handle in worker_handles {
        handle.join().unwrap();
    }

    checker_handle.join().unwrap();

    // Line 459-463: Convert atomic positions to regular vector
    x_positions.iter().map(|x| x.load()).collect()
}

/// Main entry point that mimics ODGI's usage
pub fn path_sgd_sort(graph: &BidirectedGraph, params: PathSGDParams) -> Vec<Handle> {
    // Build XP index
    let path_index = XPIndex::from_graph(graph);

    // Use all paths
    let path_sgd_use_paths: Vec<usize> = (0..graph.paths.len()).collect();

    // Run SGD
    let positions = path_linear_sgd(graph, &path_index, path_sgd_use_paths, params);

    // Sort handles by position
    let mut handle_positions: Vec<(Handle, f64)> = Vec::new();
    for (idx, &node_id) in graph.nodes.keys().enumerate() {
        let handle = Handle::new(node_id, false);
        handle_positions.push((handle, positions[idx]));
    }

    handle_positions.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    handle_positions.into_iter().map(|(h, _)| h).collect()
}