use handlegraph::handle::{Edge, Handle, NodeId};
use handlegraph::handlegraph::*;
use handlegraph::hashgraph::HashGraph;
use handlegraph::mutablehandlegraph::*;
use handlegraph::pathhandlegraph::embedded_paths::*;

fn main() {
    println!("Testing handlegraph traits...");

    let mut graph = HashGraph::new();

    // Test what methods are available
    println!("Available methods on HashGraph:");

    // Try to create nodes and see what works
    let node1 = NodeId::from(1);
    let node2 = NodeId::from(2);

    // Test creating handles - this should work
    let h1 = graph.create_handle(b"A", node1);
    let h2 = graph.create_handle(b"T", node2);

    println!("Created handles: {:?} and {:?}", h1, h2);

    // Test creating edge
    graph.create_edge(&Edge(h1, h2));
    println!("Created edge");

    // Test node iteration
    println!("Nodes:");
    for node_id in graph.handles() {
        println!("  Node: {}", node_id.as_integer());
    }

    // Test getting sequence
    let seq = graph.sequence(h1);
    println!("Sequence of h1: {:?}", seq);

    // Test path creation
    match graph.create_path(b"test", false) {
        Ok(path_id) => {
            println!("Created path: {:?}", path_id);

            // Add steps
            graph.path_append_step(&path_id, h1);
            graph.path_append_step(&path_id, h2);

            // Test path iteration
            println!("Path steps:");
            for handle in graph.path(&path_id) {
                println!("  Step: {:?}", handle);
            }
        }
        Err(e) => println!("Failed to create path: {:?}", e),
    }

    println!("Done!");
}
