use seqrush::graph_ops::{Graph, Node, Edge};

fn main() {
    let mut graph = Graph::new();
    
    // Create the same graph as test_partial_chain
    for i in 1..=7 {
        graph.nodes.insert(i, Node {
            id: i,
            sequence: vec![b'A' + (i as u8 - 1)],
            rank: i as f64,
        });
    }
    
    graph.edges.insert(Edge { from: 1, to: 2 });
    graph.edges.insert(Edge { from: 2, to: 3 });
    graph.edges.insert(Edge { from: 3, to: 4 });
    graph.edges.insert(Edge { from: 4, to: 5 });
    graph.edges.insert(Edge { from: 5, to: 6 });
    graph.edges.insert(Edge { from: 7, to: 2 });
    graph.edges.insert(Edge { from: 7, to: 6 });
    
    graph.paths.push(("path1".to_string(), vec![1, 2, 3, 4, 5, 6]));
    graph.paths.push(("path2".to_string(), vec![7, 2, 3, 4, 5, 6]));
    graph.paths.push(("path3".to_string(), vec![7, 6]));
    
    println!("Before compaction: {} nodes", graph.nodes.len());
    
    // Try to find simple components
    let components = graph.find_simple_components();
    println!("Found {} components", components.len());
    for (i, comp) in components.iter().enumerate() {
        println!("  Component {}: {:?}", i, comp);
    }
    
    let compacted = graph.compact_nodes();
    println!("Compacted {} nodes", compacted);
    println!("After compaction: {} nodes", graph.nodes.len());
}