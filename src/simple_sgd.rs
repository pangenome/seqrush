use crate::bidirected_ops::BidirectedGraph;
use crate::bidirected_graph::Handle;
use rand::{thread_rng, Rng};
use std::collections::HashMap;
use std::sync::{Arc, atomic::{AtomicBool, AtomicU64, Ordering}};
use std::sync::atomic::AtomicU64 as AtomicF64;
use std::thread;
use std::time::Duration;

/// Simple SGD that just tries to place connected nodes close together
pub struct SimpleSGD {
    graph: Arc<BidirectedGraph>,
    t_max: u64,
    eps: f64,
    delta: f64,
    nthreads: usize,
}

impl SimpleSGD {
    pub fn new(
        graph: Arc<BidirectedGraph>,
        t_max: u64,
        eps: f64,
        delta: f64,
        nthreads: usize,
    ) -> Self {
        SimpleSGD {
            graph,
            t_max,
            eps,
            delta,
            nthreads,
        }
    }

    pub fn run(&self) -> HashMap<usize, f64> {
        // Initialize positions linearly
        let x = Arc::new(self.initialize_positions());

        // Collect all edges
        let edges: Vec<_> = self.graph.edges.iter().cloned().collect();
        if edges.is_empty() {
            return x.iter().map(|(k, v)| (*k, u64_to_f64(v.load(Ordering::Relaxed)))).collect();
        }

        let edges = Arc::new(edges);
        eprintln!("[SimpleSGD] Optimizing {} edges for {} iterations", edges.len(), self.t_max);

        // Learning rate
        let eta = Arc::new(AtomicF64::new(f64_to_u64(1.0)));
        let delta_max = Arc::new(AtomicF64::new(f64_to_u64(f64::MAX)));
        let work_todo = Arc::new(AtomicBool::new(true));
        let iteration = Arc::new(AtomicU64::new(0));

        // Checker thread
        let checker = {
            let delta_max = Arc::clone(&delta_max);
            let work_todo = Arc::clone(&work_todo);
            let iteration = Arc::clone(&iteration);
            let eta = Arc::clone(&eta);
            let t_max = self.t_max;
            let delta = self.delta;

            thread::spawn(move || {
                // Reset delta_max to 0 initially
                delta_max.store(f64_to_u64(0.0), Ordering::Relaxed);

                while work_todo.load(Ordering::Relaxed) {
                    thread::sleep(Duration::from_millis(500));

                    let iter = iteration.fetch_add(1, Ordering::Relaxed);
                    let current_delta = u64_to_f64(delta_max.load(Ordering::Relaxed));

                    if iter % 50 == 0 {
                        eprintln!("[SimpleSGD] Iteration {} delta: {:.6}", iter, current_delta);
                    }

                    if iter >= t_max || (iter > 5 && current_delta <= delta) {
                        eprintln!("[SimpleSGD] Stopping at iteration {} with delta {:.6}", iter, current_delta);
                        work_todo.store(false, Ordering::Relaxed);
                        break;
                    }

                    // Update learning rate (exponential decay)
                    let new_eta = 1.0 * (0.95_f64).powi(iter as i32);
                    eta.store(f64_to_u64(new_eta), Ordering::Relaxed);

                    // Reset delta for next iteration
                    delta_max.store(f64_to_u64(0.0), Ordering::Relaxed);
                }
            })
        };

        // Worker threads
        let mut workers = vec![];
        for tid in 0..self.nthreads {
            let x = Arc::clone(&x);
            let edges = Arc::clone(&edges);
            let work_todo = Arc::clone(&work_todo);
            let delta_max = Arc::clone(&delta_max);
            let eta = Arc::clone(&eta);

            let handle = thread::spawn(move || {
                let mut rng = thread_rng();
                let mut updates = 0u64;

                while work_todo.load(Ordering::Relaxed) {
                    // Pick random edge
                    let edge = &edges[rng.gen_range(0..edges.len())];
                    let i = edge.from.node_id();
                    let j = edge.to.node_id();

                    // Get positions
                    let xi = x.get(&i).map(|p| u64_to_f64(p.load(Ordering::Relaxed))).unwrap_or(0.0);
                    let xj = x.get(&j).map(|p| u64_to_f64(p.load(Ordering::Relaxed))).unwrap_or(0.0);

                    // We want connected nodes to be distance 1 apart
                    let desired_distance = 1.0;
                    let actual_distance = (xi - xj).abs();

                    // SGD update
                    let learning_rate = u64_to_f64(eta.load(Ordering::Relaxed));
                    let error = actual_distance - desired_distance;
                    let update = learning_rate * error * 0.5;

                    // Track max delta
                    let update_abs = update.abs();
                    loop {
                        let current_max = u64_to_f64(delta_max.load(Ordering::Relaxed));
                        if update_abs <= current_max {
                            break;
                        }
                        if delta_max.compare_exchange_weak(
                            f64_to_u64(current_max),
                            f64_to_u64(update_abs),
                            Ordering::Relaxed,
                            Ordering::Relaxed,
                        ).is_ok() {
                            break;
                        }
                    }

                    // Apply update (pull nodes together or push apart)
                    if actual_distance > 0.01 {
                        let direction = if xi < xj { 1.0 } else { -1.0 };

                        if let Some(pos_i) = x.get(&i) {
                            loop {
                                let old_u64 = pos_i.load(Ordering::Relaxed);
                                let old = u64_to_f64(old_u64);
                                let new = old + direction * update;
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

                        if let Some(pos_j) = x.get(&j) {
                            loop {
                                let old_u64 = pos_j.load(Ordering::Relaxed);
                                let old = u64_to_f64(old_u64);
                                let new = old - direction * update;
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

                    updates += 1;
                }
            });

            workers.push(handle);
        }

        // Wait for completion
        for w in workers {
            w.join().unwrap();
        }
        checker.join().unwrap();

        // Return final positions
        x.iter().map(|(k, v)| (*k, u64_to_f64(v.load(Ordering::Relaxed)))).collect()
    }

    fn initialize_positions(&self) -> HashMap<usize, AtomicF64> {
        let mut positions = HashMap::new();
        let mut current_pos = 0.0;

        let mut node_ids: Vec<_> = self.graph.nodes.iter()
            .enumerate()
            .filter_map(|(id, n)| if n.is_some() { Some(id) } else { None })
            .collect();
        node_ids.sort();

        for node_id in node_ids {
            positions.insert(node_id, AtomicF64::new(f64_to_u64(current_pos)));
            let handle = Handle::new(node_id, false);
            if let Some(seq) = self.graph.get_sequence(handle) {
                current_pos += seq.len() as f64;
            } else {
                current_pos += 1.0;
            }
        }

        positions
    }

    pub fn get_order(&self) -> Vec<Handle> {
        let positions = self.run();

        let mut layout_handles: Vec<_> = positions
            .iter()
            .map(|(&node_id, &pos)| (pos, Handle::new(node_id, false)))
            .collect();

        layout_handles.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        layout_handles.into_iter().map(|(_, h)| h).collect()
    }
}

// Helper functions
fn f64_to_u64(f: f64) -> u64 {
    f.to_bits()
}

fn u64_to_f64(u: u64) -> f64 {
    f64::from_bits(u)
}