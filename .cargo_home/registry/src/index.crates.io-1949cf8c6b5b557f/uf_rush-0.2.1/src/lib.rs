use std::sync::atomic::{AtomicUsize, Ordering};

/// Constant defining the number of `rank` bits in a node represented as a `usize`.
const RANK_BITS: u32 = usize::BITS.ilog2();

/// Constant defining the number of `parent` bits in a node represented as a `usize`.
const PARENT_BITS: u32 = usize::BITS - RANK_BITS;

/// Maximum allowable size of a lock-free union-find data structure.
pub const MAX_SIZE: usize = usize::MAX >> RANK_BITS;

/// Thread-safe and lock-free implementation of a union-find (also known as disjoint set) data
/// structure.
///
/// This implementation is based on the algorithm presented in
///
/// > "Wait-free Parallel Algorithms for the Union-Find Problem" \
/// > by Richard J. Anderson and Heather Woll.
pub struct UFRush {
    /// List of nodes in the union-find structure, represented by atomic unsigned integers.
    nodes: Vec<AtomicUsize>,
}

/// Implementation block for the UFRush struct.
impl UFRush {
    /// Creates a new union-find data structure with a specified number of elements.
    ///
    /// # Arguments
    /// * `size` - Number of elements in the union-find structure.
    ///
    /// # Returns
    /// An instance of [`UFRush`].
    ///
    /// # Panics
    /// This method will panic if the `size` exceeds the [`MAX_SIZE`].
    pub fn new(size: usize) -> Self {
        assert!(size <= MAX_SIZE);

        Self {
            nodes: (0..size).map(AtomicUsize::new).collect(),
        }
    }

    /// Returns the total number of elements in the union-find structure.
    ///
    /// # Returns
    /// The total number of elements.
    pub fn size(&self) -> usize {
        self.nodes.len()
    }

    /// Determines whether the elements `x` and `y` belong to the same subset.
    ///
    /// # Arguments
    /// * `x` - The first element.
    /// * `y` - The second element.
    ///
    /// # Returns
    /// [`true`] if `x` and `y` belong to the same subset; [`false`] otherwise.
    ///
    /// # Panics
    /// This method will panic if `x` or `y` are out of bounds.
    ///
    /// # Note
    /// The same operation checks whether two elements belong to the same subset. In a sequential
    /// scenario, this operation could be considered redundant, as it can be constructed from a pair
    /// of find operations. However, when it comes to concurrent environments, providing same as a
    /// basic operation is crucial. This is because in such scenarios, the identifiers of subsets
    /// might change dynamically due to concurrent union operations, making it challenging to reliably
    /// determine if a pair of elements belong to the same subset solely based on the outcomes of
    /// individual `find` operations.
    pub fn same(&self, x: usize, y: usize) -> bool {
        loop {
            let x_rep = self.find(x);
            let y_rep = self.find(y);
            if x_rep == y_rep {
                return true;
            }
            let x_node = self.nodes[x_rep].load(Ordering::Relaxed);
            if x_rep == parent(x_node) {
                return false;
            }
        }
    }

    /// Finds the representative of the subset that `x` belongs to.
    ///
    /// # Arguments
    /// * `x` - The element to find the representative for.
    ///
    /// # Returns
    /// The representative element of the subset that contains `x`.
    ///
    /// # Panics
    /// This method will panic if `x` is out of bounds.
    ///
    /// # Note
    /// The find operation uses the "path halving" technique, an intermediate strategy between full
    /// path compression and no compression at all.
    ///
    /// In the path halving technique, instead of making every node in the path point directly to the
    /// root as in full path compression, we only change the parent of every other node in the path to
    /// point to its grandparent. This is achieved by skipping over the parent node on each iteration
    /// during the find operation. Despite not fully compressing the path, this strategy is still
    /// effective in flattening the tree structure over time, thus accelerating future operations.
    ///
    /// The advantage of path halving is that it achieves a good balance between the speed of the find
    /// operation and the amount of modification it makes to the tree structure, avoiding a potential
    /// slowdown due to excessively frequent writes in highly concurrent scenarios. Therefore, it is
    /// particularly suitable for lock-free data structures like [`UFRush`], where minimizing
    /// write contention is crucial for performance.
    pub fn find(&self, mut x: usize) -> usize {
        assert!(x < self.size());

        let mut x_node = self.nodes[x].load(Ordering::Relaxed);
        while x != parent(x_node) {
            let x_parent = parent(x_node);
            let x_parent_node = self.nodes[x_parent].load(Ordering::Relaxed);
            let x_parent_parent = parent(x_parent_node);

            let x_new_node = encode(x_parent_parent, rank(x_node));
            let _ = self.nodes[x].compare_exchange_weak(
                x_node,
                x_new_node,
                Ordering::Release,
                Ordering::Relaxed,
            );

            x = x_parent_parent;
            x_node = self.nodes[x].load(Ordering::Relaxed);
        }
        x
    }

    /// Unites the subsets that contain `x` and `y`.
    ///
    /// If `x` and `y` are already in the same subset, no action is performed.
    ///
    /// # Arguments
    /// * `x` - The first element.
    /// * `y` - The second element.
    ///
    /// # Returns
    /// [`true`] if `x` and `y` were in different subsets and a union operation was performed;
    /// [`false`] if `x` and `y` were already in the same subset.
    ///
    /// # Panics
    /// This method will panic if `x` or `y` are out of bounds.
    ///
    /// # Note
    /// The unite operation utilizes a Union-Find algorithm that adopts the "union by rank"
    /// strategy for its union operation.
    ///
    /// In "union by rank", each node holds a rank, and when two sets are united, the set with the
    /// smaller rank becomes a subset of the set with the larger rank. If both sets have the same
    /// rank, either one can become a subset of the other, but the rank of the new root is incremented
    /// by one. This strategy ensures that the tree representing the set does not become excessively
    /// deep, which helps keep the operation's time complexity nearly constant.
    pub fn unite(&self, x: usize, y: usize) -> bool {
        loop {
            // Load representative for x and y
            let mut x_rep = self.find(x);
            let mut y_rep = self.find(y);

            // If they are already part of the same set, return false
            if x_rep == y_rep {
                return false;
            }

            // Load the encoded representation of the representatives
            let x_node = self.nodes[x_rep].load(Ordering::Relaxed);
            let y_node = self.nodes[y_rep].load(Ordering::Relaxed);

            let mut x_rank = rank(x_node);
            let mut y_rank = rank(y_node);

            // Swap the elements around to always make x the smaller one
            if x_rank > y_rank || (x_rank == y_rank && x_rep > y_rep) {
                std::mem::swap(&mut x_rep, &mut y_rep);
                std::mem::swap(&mut x_rank, &mut y_rank);
            }

            // x_rep is a root
            let cur_value = encode(x_rep, x_rank);
            // assign the new root to be y
            let new_value = encode(y_rep, x_rank);
            // change the value of the smaller subtree root to point to the other one
            if self.nodes[x_rep]
                .compare_exchange(cur_value, new_value, Ordering::Release, Ordering::Acquire)
                .is_ok()
            {
                // x_repr now points to y_repr
                // If the subtrees has the same height, increase the rank of the new root
                if x_rank == y_rank {
                    let cur_value = encode(y_rep, y_rank);
                    let new_value = encode(y_rep, y_rank + 1);
                    let _ = self.nodes[y_rep].compare_exchange_weak(
                        cur_value,
                        new_value,
                        Ordering::Release,
                        Ordering::Relaxed,
                    );
                }
                return true;
            }
            // A different thread has already merged modified the value of x_repr -> repeat
        }
    }

    /// Clears the union-find structure, making every element a separate subset.
    pub fn clear(&mut self) {
        self.nodes
            .iter_mut()
            .enumerate()
            .for_each(|(i, node)| node.store(i, Ordering::Relaxed));
    }
}

/// This unsafe implementation indicate that [`UFRush`] can safely be shared
/// across threads (`Sync`).
unsafe impl Sync for UFRush {}

/// This unsafe implementation indicate that [`UFRush`] is safe to transfer
/// the ownership between threads (`Send`).
unsafe impl Send for UFRush {}

/// Encodes the parent node and rank into a single `usize`.
fn encode(parent: usize, rank: usize) -> usize {
    parent | (rank << PARENT_BITS)
}

/// Retrieves the parent node from an encoded `usize`.
fn parent(n: usize) -> usize {
    n & MAX_SIZE
}

/// Retrieves the rank from an encoded `usize`.
fn rank(n: usize) -> usize {
    n >> PARENT_BITS
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::prelude::*;
    use std::collections::HashSet;
    use std::sync::Arc;
    use std::thread;

    #[test]
    fn test_new() {
        let uf = UFRush::new(10);
        assert_eq!(uf.size(), 10);
    }

    #[test]
    fn test_find() {
        let uf = UFRush::new(10);
        assert_eq!(uf.find(5), 5);
    }

    #[test]
    fn test_same() {
        let uf = UFRush::new(10);
        assert!(!uf.same(1, 2));
    }

    #[test]
    fn test_unite() {
        let uf = UFRush::new(10);
        assert!(!uf.same(1, 2));
        assert!(uf.unite(1, 2));
        assert!(uf.same(1, 2));
    }

    #[test]
    fn test_unite_already_united() {
        let uf = UFRush::new(10);
        assert!(uf.unite(1, 2));
        assert!(!uf.unite(1, 2));
    }

    #[test]
    fn test_clear() {
        let mut uf = UFRush::new(10);
        assert!(uf.unite(1, 2));
        assert!(uf.same(1, 2));
        uf.clear();
        assert!(!uf.same(1, 2));
    }

    #[test]
    fn test_multithreaded_build_cyclic_graph() {
        let vertices = 100;
        let uf = Arc::new(UFRush::new(vertices));

        // Spawn threads, one thread per nodes
        let handles: Vec<_> = (0..vertices)
            .map(|n| {
                let uf = Arc::clone(&uf);
                thread::spawn(move || {
                    // Unite the current node with the next one, creating a cycle
                    uf.unite(n, (n + 1) % vertices);
                })
            })
            .collect();

        // Wait for all threads to finish
        for handle in handles {
            handle.join().unwrap();
        }

        // Check results - all nodes should be in the same subset
        for n in 0..vertices - 1 {
            assert!(uf.same(n, (n + 1) % vertices));
        }
    }

    #[test]
    fn test_multithreaded_cyclic_graph() {
        assert!(is_cyclic(3, [(0, 1), (1, 2), (2, 0)]));
    }

    #[test]
    fn test_multithreaded_acyclic_graph() {
        assert!(!is_cyclic(4, [(0, 1), (1, 2), (2, 3)]));
    }

    #[test]
    fn stress_test() {
        let vertices = 100;
        let mut edges = HashSet::with_capacity(5 * vertices);
        let mut rng = rand::thread_rng();

        // Add edges to form a cycle
        edges.extend((0..vertices).map(|n| (n, (n + 1) % vertices)));

        // Add some extra random edges
        for _ in edges.len()..edges.capacity() {
            let u = rng.gen_range(0..vertices);
            let v = rng.gen_range(0..vertices);
            if u != v {
                edges.insert((u, v));
            }
        }

        // Validate the number of successful unite actions for some random arrangement of edges
        let mut edges: Vec<_> = edges.into_iter().collect();
        for _ in 0..100 {
            let uf = Arc::new(UFRush::new(vertices));

            // Spawn threads, one thread per edge
            let handles: Vec<_> = edges
                .iter()
                .map(|&(u, v)| {
                    let uf = Arc::clone(&uf);
                    thread::spawn(move || uf.unite(u, v))
                })
                .collect();

            // Wait for all threads to finish; count the number of united edges
            let total_united = handles
                .into_iter()
                .map(|handle| handle.join().unwrap())
                .filter(|&united| united)
                .count();

            // Ensure the number of successful unite actions is as expected
            assert_eq!(total_united, vertices - 1);

            // Shuffle the edges for the next iteration
            edges.shuffle(&mut rng);
        }
    }

    fn is_cyclic<I>(vertices: usize, edges: I) -> bool
    where
        I: IntoIterator<Item = (usize, usize)>,
    {
        let uf = Arc::new(UFRush::new(vertices));

        // Spawn threads, one thread per edge
        let handles: Vec<_> = edges
            .into_iter()
            .map(|(u, v)| {
                let uf = Arc::clone(&uf);
                thread::spawn(move || {
                    if uf.same(u, v) {
                        // If two nodes are in the same set, we've found a cycle
                        true
                    } else {
                        // Otherwise, unite them and continue
                        uf.unite(u, v);
                        false
                    }
                })
            })
            .collect();

        // Wait for all threads to finish and return if there was any cycle
        handles
            .into_iter()
            .map(|handle| handle.join().unwrap())
            .any(|cyclic| cyclic)
    }
}
