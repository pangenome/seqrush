use seqrush::embedded_graph::*;
use seqrush::bidirected_graph::Handle;

fn main() {
    let mut graph = EmbeddedGraph::new();
    
    // Create a circular structure: A -> B -> C -> A
    let a = graph.add_node(b"A".to_vec());
    let b = graph.add_node(b"B".to_vec());
    let c = graph.add_node(b"C".to_vec());
    
    let path = graph.add_path("circular".to_string());
    graph.extend_path(path, Handle::forward(a)).unwrap();
    graph.extend_path(path, Handle::forward(b)).unwrap();
    graph.extend_path(path, Handle::forward(c)).unwrap();
    graph.extend_path(path, Handle::forward(a)).unwrap(); // Back to A
    
    // Debug: print connections
    println!("Connections:");
    println!("A -> {:?}", graph.get_next_steps(Handle::forward(a)));
    println!("B -> {:?}", graph.get_next_steps(Handle::forward(b)));
    println!("C -> {:?}", graph.get_next_steps(Handle::forward(c)));
    
    println!("\nPrev connections:");
    println!("A <- {:?}", graph.get_prev_steps(Handle::forward(a)));
    println!("B <- {:?}", graph.get_prev_steps(Handle::forward(b)));
    println!("C <- {:?}", graph.get_prev_steps(Handle::forward(c)));
    
    // Check perfect neighbors
    println!("\nPerfect neighbors:");
    println!("A->B? {}", graph.are_perfect_neighbors(Handle::forward(a), Handle::forward(b)));
    println!("B->C? {}", graph.are_perfect_neighbors(Handle::forward(b), Handle::forward(c)));
    println!("C->A? {}", graph.are_perfect_neighbors(Handle::forward(c), Handle::forward(a)));
}