use seqrush::embedded_graph::*;
use seqrush::bidirected_graph::Handle;

fn main() {
    let mut graph = EmbeddedGraph::new();
    
    // Create a self-loop: A -> A -> B
    let a = graph.add_node(b"A".to_vec());
    let b = graph.add_node(b"B".to_vec());
    
    let path = graph.add_path("test".to_string());
    graph.extend_path(path, Handle::forward(a)).unwrap();
    graph.extend_path(path, Handle::forward(a)).unwrap(); // Self-loop
    graph.extend_path(path, Handle::forward(b)).unwrap();
    
    // Debug output
    println!("Node A path steps:");
    if let Some(node_a) = graph.nodes.get(&a) {
        for ((pid, is_rev), step) in &node_a.path_steps {
            println!("  Path {:?}, reverse={}, prev={:?}, next={:?}", 
                     pid, is_rev, step.prev, step.next);
        }
    }
    
    // A should have itself as both next and prev
    let a_next = graph.get_next_steps(Handle::forward(a));
    let a_prev = graph.get_prev_steps(Handle::forward(a));
    
    println!("\nA next steps: {:?}", a_next);
    println!("A prev steps: {:?}", a_prev);
    
    // Check what we actually have
    println!("\nExpected a_next to contain both A and B");
    println!("A in next? {}", a_next.contains(&Handle::forward(a)));
    println!("B in next? {}", a_next.contains(&Handle::forward(b)));
}