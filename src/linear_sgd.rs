use crate::bidirected_ops::BidirectedGraph;
use crate::bidirected_graph::Handle;
use rand::distributions::{Distribution, Uniform};
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256Plus;
use std::sync::{Arc, atomic::{AtomicBool, AtomicU64, Ordering}};
use std::sync::atomic::AtomicU64 as AtomicF64; // We'll use u64 for atomic float operations
use std::thread;
use std::time::Duration;
use std::collections::HashMap;

/// Zipfian distribution for sampling jumps along paths
/// This is a simplified version - ODGI uses a more complex cached version
struct ZipfianDistribution {
    n: u64,
    theta: f64,
    zeta: f64,
}

impl ZipfianDistribution {
    fn new(n: u64, theta: f64) -> Self {
        let mut zeta = 0.0;
        for i in 1..=n {
            zeta += 1.0 / (i as f64).powf(theta);
        }
        ZipfianDistribution { n, theta, zeta }
    }

    fn sample(&self, rng: &mut impl Rng) -> u64 {
        let u: f64 = rng.gen();
        let mut sum = 0.0;
        for i in 1..=self.n {
            sum += (1.0 / (i as f64).powf(self.theta)) / self.zeta;
            if sum >= u {
                return i;
            }
        }
        self.n
    }
}

/// Path-guided stochastic gradient descent for graph linearization
pub struct LinearSGD {
    graph: Arc<BidirectedGraph>,
    bandwidth: u64,
    sampling_rate: f64,
    use_paths: bool,
    t_max: u64,
    eps: f64,
    delta: f64,
    nthreads: usize,
}

impl LinearSGD {
    pub fn new(
        graph: Arc<BidirectedGraph>,
        bandwidth: u64,
        sampling_rate: f64,
        use_paths: bool,
        t_max: u64,
        eps: f64,
        delta: f64,
        nthreads: usize,
    ) -> Self {
        LinearSGD {
            graph,
            bandwidth,
            sampling_rate,
            use_paths,
            t_max,
            eps,
            delta,
            nthreads,
        }
    }

    /// Run the SGD algorithm and return node positions (ODGI-style implementation)
    pub fn run(&self) -> HashMap<usize, f64> {
        let num_nodes = self.graph.nodes.len();
        if num_nodes == 0 {
            return HashMap::new();
        }

        // Initialize positions based on current graph order (matching ODGI)
        let x = Arc::new(self.initialize_atomic_positions());

        // Check if we have any paths with more than one step
        let mut at_least_one_path_with_steps = false;
        for path in &self.graph.paths {
            if path.steps.len() > 1 {
                at_least_one_path_with_steps = true;
                break;
            }
        }

        if !at_least_one_path_with_steps {
            eprintln!("[SGD] No paths with multiple steps found, skipping optimization");
            return x.iter().map(|(k, v)| (*k, u64_to_f64(v.load(Ordering::Relaxed)))).collect();
        }

        // Build a mapping of all path steps and their positions for random sampling
        let mut all_steps = Vec::new();
        let mut path_step_positions = Vec::new(); // Positions for each path

        for (path_idx, path) in self.graph.paths.iter().enumerate() {
            let mut step_positions = Vec::new();
            let mut cumulative_pos = 0.0;

            for (step_idx, handle) in path.steps.iter().enumerate() {
                step_positions.push(cumulative_pos);
                all_steps.push((path_idx, step_idx));

                if let Some(seq) = self.graph.get_sequence(*handle) {
                    cumulative_pos += seq.len() as f64;
                }
            }
            path_step_positions.push(step_positions);
        }

        if all_steps.is_empty() {
            eprintln!("[SGD] No path steps found");
            return x.iter().map(|(k, v)| (*k, u64_to_f64(v.load(Ordering::Relaxed)))).collect();
        }

        let all_steps = Arc::new(all_steps);
        let path_step_positions = Arc::new(path_step_positions);
        let min_term_updates = 1000; // ODGI default
        let total_term_updates = self.t_max * min_term_updates;

        // Parameters matching ODGI defaults
        let theta = 0.99; // ODGI default for Zipfian distribution
        let space = 100; // ODGI uses larger space by default
        let cooling_start = 0.5; // Start cooling at 50% of iterations
        let first_cooling_iteration = (cooling_start * self.t_max as f64) as u64;

        // Shared state
        let term_updates = Arc::new(AtomicU64::new(0));
        let iteration = Arc::new(AtomicU64::new(0));
        let work_todo = Arc::new(AtomicBool::new(true));
        let cooling = Arc::new(AtomicBool::new(false));
        let delta_max = Arc::new(AtomicF64::new(f64_to_u64(0.0)));

        // Learning rate schedule (matching ODGI)
        let w_min = 1.0 / (space as f64 * space as f64);
        let w_max = 1.0;
        let eta_max = 1.0 / w_min;
        let eta_min = self.eps / w_max;
        let lambda = (eta_max / eta_min).ln() / (self.t_max as f64 - 1.0);
        let eta = Arc::new(AtomicF64::new(f64_to_u64(eta_max)));

        eprintln!("[SGD] Starting path-guided SGD with {} path steps, {} iterations", all_steps.len(), self.t_max);

        // Checker thread (matching ODGI)
        let checker_handle = {
            let term_updates = Arc::clone(&term_updates);
            let iteration = Arc::clone(&iteration);
            let work_todo = Arc::clone(&work_todo);
            let cooling = Arc::clone(&cooling);
            let delta_max = Arc::clone(&delta_max);
            let eta = Arc::clone(&eta);
            let t_max = self.t_max;
            let delta = self.delta;

            thread::spawn(move || {
                while work_todo.load(Ordering::Relaxed) {
                    let updates = term_updates.load(Ordering::Relaxed);

                    if updates >= min_term_updates {
                        let iter = iteration.fetch_add(1, Ordering::Relaxed);

                        if iter >= first_cooling_iteration {
                            cooling.store(true, Ordering::Relaxed);
                        }

                        let current_delta_max = u64_to_f64(delta_max.load(Ordering::Relaxed));

                        if iter >= t_max || (iter > 10 && current_delta_max <= delta) {
                            work_todo.store(false, Ordering::Relaxed);
                        } else {
                            // Update learning rate
                            let new_eta = eta_max * (-lambda * iter as f64).exp();
                            eta.store(f64_to_u64(new_eta), Ordering::Relaxed);
                            delta_max.store(f64_to_u64(0.0), Ordering::Relaxed);
                        }

                        term_updates.store(0, Ordering::Relaxed);
                    }

                    thread::sleep(Duration::from_millis(1));
                }
            })
        };

        // Worker threads (matching ODGI)
        let mut worker_handles = vec![];
        for tid in 0..self.nthreads {
            let x = Arc::clone(&x);
            let all_steps = Arc::clone(&all_steps);
            let path_step_positions = Arc::clone(&path_step_positions);
            let term_updates = Arc::clone(&term_updates);
            let work_todo = Arc::clone(&work_todo);
            let cooling = Arc::clone(&cooling);
            let delta_max = Arc::clone(&delta_max);
            let eta = Arc::clone(&eta);
            let graph = Arc::clone(&self.graph);

            let handle = thread::spawn(move || {
                // Each thread gets its own PRNG with unique seed (matching ODGI)
                let seed = 9399220 + tid as u64;
                let mut rng = Xoshiro256Plus::seed_from_u64(seed);
                let step_dist = Uniform::new(0, all_steps.len());
                let flip_dist = Uniform::new(0, 2);

                // Simple Zipfian for now (ODGI uses a more complex cached version)
                let zipf = ZipfianDistribution::new(space.min(100), if theta > 0.0 { theta } else { 0.01 });

                // Debug counters
                let mut updates_made = 0;
                let mut total_error = 0.0;

                while work_todo.load(Ordering::Relaxed) {
                    // Sample a random step from all path steps
                    let step_idx = step_dist.sample(&mut rng);
                    let (path_idx, step_a_idx) = all_steps[step_idx];
                    let path = &graph.paths[path_idx];

                    if path.steps.len() <= 1 {
                        continue;
                    }

                    // Determine step_b based on cooling and Zipfian distribution
                    let step_b_idx = if cooling.load(Ordering::Relaxed) || flip_dist.sample(&mut rng) == 1 {
                        // Use Zipfian distribution to favor nearby nodes
                        if step_a_idx > 0 && (flip_dist.sample(&mut rng) == 1 || step_a_idx == path.steps.len() - 1) {
                            // Go backward
                            let jump_space = (step_a_idx as u64).min(space);
                            let z = zipf.sample(&mut rng).min(jump_space);
                            step_a_idx.saturating_sub(z as usize)
                        } else {
                            // Go forward
                            let remaining = path.steps.len() - step_a_idx - 1;
                            let jump_space = (remaining as u64).min(space);
                            let z = zipf.sample(&mut rng).min(jump_space);
                            (step_a_idx + z as usize).min(path.steps.len() - 1)
                        }
                    } else {
                        // Random step in the path
                        let uniform = Uniform::new(0, path.steps.len());
                        uniform.sample(&mut rng)
                    };

                    if step_a_idx == step_b_idx {
                        continue;
                    }

                    // Get the handles
                    let handle_i = path.steps[step_a_idx];
                    let handle_j = path.steps[step_b_idx];
                    let node_i = handle_i.node_id();
                    let node_j = handle_j.node_id();

                    if node_i == node_j {
                        continue;
                    }

                    // Get precomputed positions (matching ODGI's get_position_of_step)
                    let pos_a = path_step_positions[path_idx][step_a_idx];
                    let pos_b = path_step_positions[path_idx][step_b_idx];

                    // Distance is the absolute difference in positions (matching ODGI)
                    let d_ij = (pos_a - pos_b).abs();

                    if d_ij == 0.0 {
                        continue;
                    }

                    // Weight based on distance (matching ODGI: 1/d not 1/d²)
                    let w_ij = 1.0 / d_ij;

                    // Get current positions
                    let xi = x.get(&node_i).map(|p| u64_to_f64(p.load(Ordering::Relaxed))).unwrap_or(0.0);
                    let xj = x.get(&node_j).map(|p| u64_to_f64(p.load(Ordering::Relaxed))).unwrap_or(0.0);

                    // Calculate update (matching ODGI exactly)
                    let dx = xi - xj;
                    let mag = if dx < 0.0 { -dx } else { dx };

                    // Learning rate
                    let eta_val = u64_to_f64(eta.load(Ordering::Relaxed));
                    let mu = (eta_val * w_ij).min(1.0);

                    // SGD update
                    let delta = mu * (mag - d_ij) / 2.0;
                    let delta_abs = if delta < 0.0 { -delta } else { delta };

                    // Track max delta
                    loop {
                        let current_max = u64_to_f64(delta_max.load(Ordering::Relaxed));
                        if delta_abs <= current_max {
                            break;
                        }
                        if delta_max.compare_exchange_weak(
                            f64_to_u64(current_max),
                            f64_to_u64(delta_abs),
                            Ordering::Relaxed,
                            Ordering::Relaxed,
                        ).is_ok() {
                            break;
                        }
                    }

                    // Debug: Track error
                    total_error += (mag - d_ij).abs();
                    updates_made += 1;

                    // Apply update
                    if mag > 1e-9 {
                        let r = delta / mag;
                        let r_x = r * dx;

                        if let Some(pos_i) = x.get(&node_i) {
                            loop {
                                let old_u64 = pos_i.load(Ordering::Relaxed);
                                let old = u64_to_f64(old_u64);
                                let new = old - r_x;
                                if pos_i.compare_exchange_weak(
                                    old_u64,
                                    f64_to_u64(new),
                                    Ordering::Relaxed,
                                    Ordering::Relaxed,
                                ).is_ok() {
                                    break;
                                }
                            }
                        }
                        if let Some(pos_j) = x.get(&node_j) {
                            loop {
                                let old_u64 = pos_j.load(Ordering::Relaxed);
                                let old = u64_to_f64(old_u64);
                                let new = old + r_x;
                                if pos_j.compare_exchange_weak(
                                    old_u64,
                                    f64_to_u64(new),
                                    Ordering::Relaxed,
                                    Ordering::Relaxed,
                                ).is_ok() {
                                    break;
                                }
                            }
                        }
                    }

                    term_updates.fetch_add(1, Ordering::Relaxed);
                }

                // Debug output removed for cleaner runs
            });

            worker_handles.push(handle);
        }

        // Wait for completion
        for handle in worker_handles {
            handle.join().unwrap();
        }
        checker_handle.join().unwrap();

        let final_iter = iteration.load(Ordering::Relaxed);
        let final_delta = u64_to_f64(delta_max.load(Ordering::Relaxed));
        let total_updates = term_updates.load(Ordering::Relaxed);
        eprintln!("[SGD] Complete: {} iterations, {} term updates, final delta: {:.6}",
                 final_iter, total_updates, final_delta);

        // Extract final positions and check if they make sense
        let final_positions: HashMap<usize, f64> = x.iter()
            .map(|(k, v)| (*k, u64_to_f64(v.load(Ordering::Relaxed))))
            .collect();

        // Debug: Check the range of positions
        let mut min_pos = f64::MAX;
        let mut max_pos = f64::MIN;
        let mut negative_count = 0;
        for &pos in final_positions.values() {
            min_pos = min_pos.min(pos);
            max_pos = max_pos.max(pos);
            if pos < 0.0 {
                negative_count += 1;
            }
        }
        eprintln!("[SGD] Final position range: {:.2} to {:.2} (span: {:.2}), {} negative positions",
                 min_pos, max_pos, max_pos - min_pos, negative_count);

        // Debug: Check how many nodes moved significantly
        let mut moved_count = 0;
        let mut total_movement = 0.0;
        for (node_id, &final_pos) in &final_positions {
            // Find initial position
            let mut initial_pos = 0.0;
            let mut found = false;
            let mut cumulative = 0.0;
            let mut sorted_ids: Vec<_> = self.graph.nodes.iter()
                .enumerate()
                .filter_map(|(id, n)| if n.is_some() { Some(id) } else { None })
                .collect();
            sorted_ids.sort();
            for id in sorted_ids {
                if id == *node_id {
                    initial_pos = cumulative;
                    found = true;
                    break;
                }
                if let Some(seq) = self.graph.get_sequence(Handle::new(id, false)) {
                    cumulative += seq.len() as f64;
                }
            }
            if found {
                let movement = (final_pos - initial_pos).abs();
                if movement > 1.0 {
                    moved_count += 1;
                    total_movement += movement;
                }
            }
        }
        eprintln!("[SGD] {} nodes moved >1.0 units, avg movement: {:.2}",
                 moved_count,
                 if moved_count > 0 { total_movement / moved_count as f64 } else { 0.0 });

        final_positions
    }

    /// Initialize atomic node positions based on graph order (matching ODGI)
    fn initialize_atomic_positions(&self) -> HashMap<usize, AtomicF64> {
        let mut positions = HashMap::new();
        let mut current_pos = 0.0;

        // Get nodes in sorted order (matching ODGI's approach)
        let mut node_ids: Vec<usize> = self.graph.nodes.iter()
            .enumerate()
            .filter_map(|(id, n)| if n.is_some() { Some(id) } else { None })
            .collect();
        node_ids.sort();

        for node_id in node_ids {
            positions.insert(node_id, AtomicF64::new(f64_to_u64(current_pos)));
            let handle = Handle::new(node_id, false);
            if let Some(seq) = self.graph.get_sequence(handle) {
                current_pos += seq.len() as f64;
            }
        }

        positions
    }


    /// Get the sorted node order after SGD
    pub fn get_order(&self) -> Vec<Handle> {
        let positions = self.run();

        // Create pairs of (position, handle)
        let mut layout_handles: Vec<(f64, Handle)> = positions
            .iter()
            .map(|(&node_id, &pos)| {
                let handle = Handle::new(node_id, false);
                (pos, handle)
            })
            .collect();

        // Sort by position
        layout_handles.sort_by(|a, b| {
            a.0.partial_cmp(&b.0).unwrap()
                .then_with(|| a.1.node_id().cmp(&b.1.node_id()))
        });

        layout_handles.into_iter().map(|(_, h)| h).collect()
    }
}

// Helper functions for atomic float operations
fn f64_to_u64(f: f64) -> u64 {
    f.to_bits()
}

fn u64_to_f64(u: u64) -> f64 {
    f64::from_bits(u)
}