/// Exact implementation of ODGI's path_linear_sgd from path_sgd.cpp
/// This is the algorithm used by `odgi sort -p Ygs`
use crate::bidirected_ops::BidirectedGraph;
use crate::bidirected_graph::Handle;
use rand::distributions::{Distribution, Uniform};
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256Plus;
use std::collections::HashMap;
use std::sync::{Arc, atomic::{AtomicBool, AtomicU64, Ordering}};
use std::thread;
use std::time::Duration;

/// Path index structure - simplified version of ODGI's XP
/// In the real implementation, this uses SDSL bit vectors for efficiency
pub struct PathIndex {
    /// For each path step, stores the handle it refers to
    step_to_handle: Vec<Handle>,
    /// For each path step, stores its position in the path (in bp)
    step_to_position: Vec<usize>,
    /// For each path step, stores which path it belongs to
    step_to_path: Vec<usize>,
    /// For each path step, stores its rank in the path (0-indexed)
    step_to_rank: Vec<usize>,
    /// Path metadata
    paths: Vec<PathInfo>,
    /// Map from path name to index
    path_name_to_idx: HashMap<String, usize>,
}

#[derive(Clone)]
struct PathInfo {
    name: String,
    step_count: usize,
    length: usize, // in bp
    first_step: usize, // index in step arrays
}

impl PathIndex {
    pub fn from_graph(graph: &BidirectedGraph) -> Self {
        let mut step_to_handle = Vec::new();
        let mut step_to_position = Vec::new();
        let mut step_to_path = Vec::new();
        let mut step_to_rank = Vec::new();
        let mut paths = Vec::new();
        let mut path_name_to_idx = HashMap::new();

        for (path_idx, path) in graph.paths.iter().enumerate() {
            path_name_to_idx.insert(path.name.clone(), path_idx);
            let first_step = step_to_handle.len();
            let mut position = 0;

            for (rank, &handle) in path.steps.iter().enumerate() {
                step_to_handle.push(handle);
                step_to_position.push(position);
                step_to_path.push(path_idx);
                step_to_rank.push(rank);

                // Add node length to position
                if let Some(node) = graph.nodes.get(&handle.node_id()) {
                    position += node.sequence.len();
                }
            }

            paths.push(PathInfo {
                name: path.name.clone(),
                step_count: path.steps.len(),
                length: position,
                first_step,
            });
        }

        PathIndex {
            step_to_handle,
            step_to_position,
            step_to_path,
            step_to_rank,
            paths,
            path_name_to_idx,
        }
    }

    pub fn get_total_steps(&self) -> usize {
        self.step_to_handle.len()
    }

    pub fn get_handle_of_step(&self, step_idx: usize) -> Handle {
        self.step_to_handle[step_idx]
    }

    pub fn get_position_of_step(&self, step_idx: usize) -> usize {
        self.step_to_position[step_idx]
    }

    pub fn get_path_of_step(&self, step_idx: usize) -> usize {
        self.step_to_path[step_idx]
    }

    pub fn get_rank_of_step(&self, step_idx: usize) -> usize {
        self.step_to_rank[step_idx]
    }

    pub fn get_path_step_count(&self, path_idx: usize) -> usize {
        self.paths[path_idx].step_count
    }

    pub fn get_step_at_path_position(&self, path_idx: usize, rank: usize) -> usize {
        self.paths[path_idx].first_step + rank
    }

    pub fn num_paths(&self) -> usize {
        self.paths.len()
    }

    pub fn get_path_length(&self, path_idx: usize) -> usize {
        self.paths[path_idx].length
    }
}

/// Dirty Zipfian distribution with cached zeta values
/// This is a port of ODGI's dirty_zipfian_int_distribution
struct DirtyZipfian {
    min: u64,
    max: u64,
    theta: f64,
    zeta: f64,
}

impl DirtyZipfian {
    fn new(min: u64, max: u64, theta: f64, zeta: f64) -> Self {
        DirtyZipfian { min, max, theta, zeta }
    }

    fn sample(&self, rng: &mut impl Rng) -> u64 {
        let u: f64 = rng.gen();
        let mut sum = 0.0;

        for i in self.min..=self.max {
            sum += (1.0 / fast_pow(i as f64, self.theta)) / self.zeta;
            if sum >= u {
                return i;
            }
        }
        self.max
    }
}

/// Fast power function for Zipfian distribution
fn fast_pow(base: f64, exp: f64) -> f64 {
    base.powf(exp)
}

/// Convert f64 to u64 bits for atomic operations
fn f64_to_u64(f: f64) -> u64 {
    f.to_bits()
}

/// Convert u64 bits to f64
fn u64_to_f64(u: u64) -> f64 {
    f64::from_bits(u)
}

/// Path-guided stochastic gradient descent parameters
#[derive(Clone)]
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
}

impl Default for PathSGDParams {
    fn default() -> Self {
        // These are ODGI's default parameters from odgi sort -p Ygs
        PathSGDParams {
            iter_max: 100,  // ODGI default for Ygs
            iter_with_max_learning_rate: 0,  // ODGI default
            min_term_updates: 100,
            delta: 0.0,
            eps: 0.01,
            eta_max: 100.0,  // ODGI default (w_min = 1/100)
            theta: 0.99,  // ODGI default for Ygs
            space: 100,
            space_max: 100,
            space_quantization_step: 10,
            cooling_start: 0.8,
            nthreads: 1,
            progress: false,
        }
    }
}

/// Main path-guided SGD implementation (exact port of ODGI's path_linear_sgd)
pub fn path_linear_sgd(
    graph: Arc<BidirectedGraph>,
    params: PathSGDParams,
) -> HashMap<usize, f64> {
    let num_nodes = graph.nodes.len();
    if num_nodes == 0 {
        return HashMap::new();
    }

    // Build path index
    let path_index = PathIndex::from_graph(&graph);

    // Check if we have any paths with more than one step
    let mut has_valid_paths = false;
    for path in &path_index.paths {
        if path.step_count > 1 {
            has_valid_paths = true;
            break;
        }
    }

    if !has_valid_paths {
        eprintln!("[path_sgd] No paths with multiple steps found");
        return HashMap::new();
    }

    // Initialize positions based on current graph order
    let x: Vec<AtomicU64> = (0..num_nodes)
        .map(|_| AtomicU64::new(0))
        .collect();

    // Seed positions with graph layout
    // IMPORTANT: Sort nodes by ID to ensure deterministic ordering!
    let mut sorted_nodes: Vec<_> = graph.nodes.iter().collect();
    sorted_nodes.sort_by_key(|(node_id, _)| *node_id);

    let mut len = 0u64;
    let mut handle_to_idx: HashMap<Handle, usize> = HashMap::new();
    let mut idx = 0;
    for (&node_id, node) in sorted_nodes {
        let handle = Handle::forward(node_id);
        handle_to_idx.insert(handle, idx);
        x[idx].store(f64_to_u64(len as f64), Ordering::Relaxed);
        len += node.sequence.len() as u64;
        idx += 1;
    }

    // Calculate first cooling iteration
    let first_cooling_iteration = (params.cooling_start * params.iter_max as f64).floor() as u64;

    // Calculate learning rate schedule
    let w_min = 1.0 / params.eta_max;
    let w_max = 1.0;
    let etas = path_linear_sgd_schedule(
        w_min,
        w_max,
        params.iter_max,
        params.iter_with_max_learning_rate,
        params.eps,
    );

    // Pre-calculate zetas for Zipfian distribution
    let zeta_size = if params.space <= params.space_max {
        params.space as usize
    } else {
        params.space_max as usize + (params.space - params.space_max) as usize / params.space_quantization_step as usize + 1
    } + 1;

    let mut zetas = vec![0.0; zeta_size];
    let mut zeta_tmp = 0.0;
    for i in 1..=params.space {
        zeta_tmp += fast_pow(1.0 / i as f64, params.theta);
        if i <= params.space_max {
            zetas[i as usize] = zeta_tmp;
        }
        if i >= params.space_max && (i - params.space_max) % params.space_quantization_step == 0 {
            let idx = params.space_max as usize + 1 + ((i - params.space_max) / params.space_quantization_step) as usize;
            if idx < zetas.len() {
                zetas[idx] = zeta_tmp;
            }
        }
    }

    // Shared state for threads
    let x = Arc::new(x);
    let path_index = Arc::new(path_index);
    let handle_to_idx = Arc::new(handle_to_idx);
    let zetas = Arc::new(zetas);
    let etas = Arc::new(etas);

    let term_updates = Arc::new(AtomicU64::new(0));
    let iteration = Arc::new(AtomicU64::new(0));
    let eta = Arc::new(AtomicU64::new(f64_to_u64(etas[0])));
    let adj_theta = Arc::new(AtomicU64::new(f64_to_u64(params.theta)));  // Adaptive theta
    let cooling = Arc::new(AtomicBool::new(false));
    let work_todo = Arc::new(AtomicBool::new(true));
    let delta_max = Arc::new(AtomicU64::new(0));

    let total_term_updates = params.iter_max * params.min_term_updates;

    // Progress tracking
    if params.progress {
        eprintln!("[path_sgd] Starting with {} iterations, {} term updates per iteration",
                 params.iter_max, params.min_term_updates);
    }

    // Checker thread - monitors progress and updates learning rate
    let checker_handle = {
        let term_updates = Arc::clone(&term_updates);
        let iteration = Arc::clone(&iteration);
        let work_todo = Arc::clone(&work_todo);
        let eta = Arc::clone(&eta);
        let adj_theta = Arc::clone(&adj_theta);
        let cooling = Arc::clone(&cooling);
        let etas = Arc::clone(&etas);
        let delta_max = Arc::clone(&delta_max);
        let min_term_updates = params.min_term_updates;
        let iter_max = params.iter_max;
        let delta = params.delta;

        thread::spawn(move || {
            while work_todo.load(Ordering::Relaxed) {
                let curr_updates = term_updates.load(Ordering::Relaxed);

                // ODGI checks: have we done enough updates for THIS iteration?
                if curr_updates > min_term_updates {
                    let curr_iter = iteration.load(Ordering::Relaxed);
                    iteration.fetch_add(1, Ordering::Relaxed);
                    let new_iter = curr_iter + 1;

                    // Check stopping conditions BEFORE updating eta
                    if new_iter > iter_max {
                        work_todo.store(false, Ordering::Relaxed);
                    } else {
                        let curr_delta = u64_to_f64(delta_max.load(Ordering::Relaxed));
                        if delta > 0.0 && curr_delta <= delta {
                            work_todo.store(false, Ordering::Relaxed);
                        } else {
                            // Update learning rate
                            if (new_iter as usize) < etas.len() {
                                eta.store(f64_to_u64(etas[new_iter as usize]), Ordering::Relaxed);
                            }

                            // Reset delta_max to delta threshold
                            delta_max.store(f64_to_u64(delta), Ordering::Relaxed);

                            // Check if we're in cooling phase
                            if new_iter > first_cooling_iteration {
                                adj_theta.store(f64_to_u64(0.001), Ordering::Relaxed);
                                cooling.store(true, Ordering::Relaxed);
                            }
                        }
                    }

                    // CRITICAL: Reset term_updates after each iteration (ODGI does this!)
                    term_updates.store(0, Ordering::Relaxed);
                }

                thread::sleep(Duration::from_millis(1));
            }
        })
    };

    // Worker threads
    let mut handles = vec![];

    for tid in 0..params.nthreads {
        let x = Arc::clone(&x);
        let path_index = Arc::clone(&path_index);
        let handle_to_idx = Arc::clone(&handle_to_idx);
        let zetas = Arc::clone(&zetas);
        let term_updates = Arc::clone(&term_updates);
        let work_todo = Arc::clone(&work_todo);
        let eta = Arc::clone(&eta);
        let adj_theta = Arc::clone(&adj_theta);
        let cooling = Arc::clone(&cooling);
        let delta_max = Arc::clone(&delta_max);
        let space = params.space;
        let space_max = params.space_max;
        let space_quantization_step = params.space_quantization_step;
        let min_term_updates = params.min_term_updates;

        let handle = thread::spawn(move || {
            let seed = 9399220 + tid as u64;
            let mut rng = Xoshiro256Plus::seed_from_u64(seed);

            let total_steps = path_index.get_total_steps();
            let step_dist = Uniform::new(0, total_steps);
            let flip_dist = Uniform::new(0, 2);

            // Track local updates, batch them to global counter
            let mut term_updates_local = 0u64;

            // ODGI workers just check work_todo, no per-thread limit
            while work_todo.load(Ordering::Relaxed) {
                // Sample a random step
                let step_idx = step_dist.sample(&mut rng);
                let path_idx = path_index.get_path_of_step(step_idx);
                let path_step_count = path_index.get_path_step_count(path_idx);

                if path_step_count == 1 {
                    continue;
                }

                let rank_a = path_index.get_rank_of_step(step_idx);
                let mut rank_b = rank_a;

                // Decide how to sample the second step
                if cooling.load(Ordering::Relaxed) || flip_dist.sample(&mut rng) == 1 {
                    // Use Zipfian distribution with adaptive theta
                    let current_theta = u64_to_f64(adj_theta.load(Ordering::Relaxed));

                    if rank_a > 0 && (flip_dist.sample(&mut rng) == 1 || rank_a == path_step_count - 1) {
                        // Go backward
                        let jump_space = space.min(rank_a as u64);
                        let space_idx = if jump_space > space_max {
                            space_max as usize + ((jump_space - space_max) / space_quantization_step) as usize + 1
                        } else {
                            jump_space as usize
                        };

                        let space_idx = space_idx.min(zetas.len() - 1);
                        let zipf = DirtyZipfian::new(1, jump_space, current_theta, zetas[space_idx]);
                        let z_i = zipf.sample(&mut rng);
                        rank_b = rank_a.saturating_sub(z_i as usize);
                    } else if rank_a < path_step_count - 1 {
                        // Go forward
                        let jump_space = space.min((path_step_count - rank_a - 1) as u64);
                        let space_idx = if jump_space > space_max {
                            space_max as usize + ((jump_space - space_max) / space_quantization_step) as usize + 1
                        } else {
                            jump_space as usize
                        };

                        let space_idx = space_idx.min(zetas.len() - 1);
                        let zipf = DirtyZipfian::new(1, jump_space, current_theta, zetas[space_idx]);
                        let z_i = zipf.sample(&mut rng);
                        rank_b = (rank_a + z_i as usize).min(path_step_count - 1);
                    }
                } else {
                    // Sample randomly across the path
                    let rank_dist = Uniform::new(0, path_step_count);
                    rank_b = rank_dist.sample(&mut rng);
                }

                if rank_a == rank_b {
                    continue;
                }

                // Get handles for the terms
                let step_a_idx = path_index.get_step_at_path_position(path_idx, rank_a);
                let step_b_idx = path_index.get_step_at_path_position(path_idx, rank_b);

                let term_i = path_index.get_handle_of_step(step_a_idx);
                let term_j = path_index.get_handle_of_step(step_b_idx);

                // Get positions in path
                let pos_a = path_index.get_position_of_step(step_a_idx) as f64;
                let pos_b = path_index.get_position_of_step(step_b_idx) as f64;

                // Calculate term distance
                let term_dist = (pos_a - pos_b).abs();
                if term_dist == 0.0 {
                    continue;
                }

                let term_weight = 1.0 / term_dist;
                let mu = u64_to_f64(eta.load(Ordering::Relaxed)) * term_weight;
                let mu = mu.min(1.0);

                // Get node indices from our mapping
                // IMPORTANT: Always use forward orientation when looking up indices,
                // since handle_to_idx only contains forward handles
                let i = handle_to_idx.get(&Handle::forward(term_i.node_id())).copied().unwrap_or(0);
                let j = handle_to_idx.get(&Handle::forward(term_j.node_id())).copied().unwrap_or(0);

                // Calculate position difference
                let x_i = u64_to_f64(x[i].load(Ordering::Relaxed));
                let x_j = u64_to_f64(x[j].load(Ordering::Relaxed));
                let mut dx = x_i - x_j;

                // ODGI uses epsilon to avoid NaN, not continue
                if dx == 0.0 {
                    dx = 1e-9;
                }

                // Calculate update
                let mag = dx.abs();
                let delta_update = mu * (mag - term_dist) / 2.0;

                // Update delta_max
                let delta_abs = delta_update.abs();
                let mut current = delta_max.load(Ordering::Relaxed);
                while delta_abs > u64_to_f64(current) {
                    match delta_max.compare_exchange_weak(
                        current,
                        f64_to_u64(delta_abs),
                        Ordering::Relaxed,
                        Ordering::Relaxed,
                    ) {
                        Ok(_) => break,
                        Err(x) => current = x,
                    }
                }

                // Apply update
                let r = delta_update / mag;
                let r_x = r * dx;

                // Update positions atomically
                let new_i = x_i - r_x;
                let new_j = x_j + r_x;
                x[i].store(f64_to_u64(new_i), Ordering::Relaxed);
                x[j].store(f64_to_u64(new_j), Ordering::Relaxed);

                // ODGI batches updates to reduce atomic contention
                term_updates_local += 1;
                if term_updates_local >= 1000 {
                    term_updates.fetch_add(term_updates_local, Ordering::Relaxed);
                    term_updates_local = 0;
                }
            }

            // Flush remaining local updates
            if term_updates_local > 0 {
                term_updates.fetch_add(term_updates_local, Ordering::Relaxed);
            }
        });

        handles.push(handle);
    }

    // Wait for all threads
    for handle in handles {
        handle.join().unwrap();
    }

    work_todo.store(false, Ordering::Relaxed);
    checker_handle.join().unwrap();

    // Convert atomic positions to final positions
    let mut positions = HashMap::new();
    for (idx, pos) in x.iter().enumerate() {
        positions.insert(idx, u64_to_f64(pos.load(Ordering::Relaxed)));
    }

    if params.progress {
        eprintln!("[path_sgd] Complete: {} term updates", term_updates.load(Ordering::Relaxed));
    }

    positions
}

/// Learning rate schedule (exact port from ODGI)
fn path_linear_sgd_schedule(
    w_min: f64,
    w_max: f64,
    iter_max: u64,
    iter_with_max_learning_rate: u64,
    eps: f64,
) -> Vec<f64> {
    let mut etas = Vec::new();

    // ODGI formula (from path_sgd.cpp lines 478-493)
    let eta_max = 1.0 / w_min;
    let eta_min = eps / w_max;
    let lambda = (eta_max / eta_min).ln() / (iter_max as f64 - 1.0);

    // Note: ODGI loops from 0 to iter_max (inclusive), so iter_max+1 values
    for t in 0..=iter_max {
        let eta = eta_max * (-lambda * ((t as i64 - iter_with_max_learning_rate as i64).abs() as f64)).exp();
        etas.push(eta);
    }

    etas
}

/// Apply path-guided SGD to graph and return sorted handles
pub fn path_sgd_sort(graph: &BidirectedGraph, params: PathSGDParams) -> Vec<Handle> {
    let graph_arc = Arc::new(graph.clone());
    let positions = path_linear_sgd(graph_arc, params);

    // Create mapping from index to handle
    // IMPORTANT: Sort nodes by ID to ensure deterministic ordering!
    let mut sorted_node_ids: Vec<_> = graph.nodes.keys().copied().collect();
    sorted_node_ids.sort();

    let mut idx_to_handle: HashMap<usize, Handle> = HashMap::new();
    for (idx, node_id) in sorted_node_ids.iter().enumerate() {
        idx_to_handle.insert(idx, Handle::forward(*node_id));
    }

    // Sort nodes by position
    let mut node_positions: Vec<(usize, f64)> = positions.into_iter().collect();
    node_positions.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    // Map back to handles
    node_positions.into_iter()
        .map(|(idx, _)| idx_to_handle.get(&idx).copied().unwrap_or(Handle::forward(0)))
        .collect()
}