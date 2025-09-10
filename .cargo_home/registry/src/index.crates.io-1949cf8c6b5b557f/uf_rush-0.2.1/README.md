# UFRush

UFRush is a lock-free, thread-safe implementation of the Union-Find (Disjoint-Set) data structure in Rust. This data structure allows efficient computation of equivalence relations (represented as disjoint sets) and supports the following operations: `find`, `unite`, and `same`.

## Background

The Union-Find data structure provides an efficient way to manage a partition of a set into disjoint subsets. It supports two main operations: union (merging two subsets into one) and find (determining which subset a particular element is in).

Lock-free programming, or non-blocking programming, is a methodology applied in concurrent computing where multiple threads operate and manipulate shared data without the need for locks. Lock-free programming achieves its goals through the use of atomic operations which can be performed independently without interruption. In the context of UFRush, being lock-free means that multiple threads can perform operations on the Union-Find structure concurrently, and the system guarantees that at least one thread will make progress, improving the throughput and performance in a multithreaded environment.

The algorithm used in UFRush is based on the research paper ["Wait-Free Parallel Algorithms for the Union-Find Problem"](https://dl.acm.org/doi/pdf/10.1145/103418.103458) by Richard J. Anderson and Heather Woll.

## Methods

UFRush offers the following methods:

- `new(size: usize) -> Self` \
  Creates a new union-find data structure with a specified number of elements.

- `size() -> usize` \
  Returns the total number of elements in the union-find structure.

- `same(x: usize, y: usize) -> bool` \
  Determines whether the elements `x` and `y` belong to the same subset.

- `find(x: usize) -> usize` \
  Finds the representative of the subset that `x` belongs to.

- `unite(x: usize, y: usize) -> bool` \
  Unites the subsets that contain `x` and `y`.

- `clear()` \
  Clears the union-find structure, making every element a separate subset.

## Usage

Here's an example of how to use UFRush:

```rust
let n = 10;
let uf = UFRush::new(n);
uf.unite(0, 1);
uf.unite(1, 2);
uf.unite(3, 4);
uf.unite(5, 6);
uf.unite(6, 7);
uf.unite(7, 8);

// Test if these elements are in the same set
assert!(uf.same(0, 2)); // true, as 0 and 2 are connected through 1
assert!(uf.same(3, 4)); // true
assert!(uf.same(5, 8)); // true, as 5, 6, 7, and 8 are all connected

// Test some elements that are not in the same set
assert!(!uf.same(0, 3)); // false, as there is no connection between 0 and 3
assert!(!uf.same(2, 5)); // false

// Now, we connect some additional elements
uf.unite(2, 5);

// And re-test
assert!(uf.same(0, 5)); // true, as now 0 and 5 are connected through 1->2->5
```

## Testing
UFRush includes comprehensive unit tests, including tests for multithreaded usage.

Run tests with:
```bash
cargo test
```

## Disclaimer

Please note that while every effort has been made to ensure the accuracy, completeness, and reliability of the UFRush library, the developer makes no warranty, guarantee or promise (express or implied) concerning the content or functionality of this library.

The implementation of the library is based on complex lock-free programming techniques and is believed to be correct to the best of the developer's knowledge. However, due to the inherent complexity of lock-free programming, there may be unknown issues or edge cases that haven't been discovered yet. As such, the library is provided "as is", without warranty of any kind.

The developer disclaims all liability for any damages, direct or indirect, special or incidental, resulting from the use of the UFRush library. Users of the library are encouraged to thoroughly test and validate its performance and correctness in their specific use cases.

If any issues are discovered, please submit an issue with detailed steps for reproducing the bug. Contributions, in the form of code improvements or suggestions, are always welcome.

Finally, while the UFRush library aims to provide high performance in multithreaded environments, it is important to understand that the actual performance may vary based on specific workload characteristics and hardware configurations. Users are encouraged to conduct performance testing under their own workloads to ensure that the library meets their performance requirements.

## License
This project is licensed under the MIT License - see the LICENSE file for details.