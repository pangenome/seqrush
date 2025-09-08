use handlegraph::handle::{Edge, Handle, NodeId};
use handlegraph::handlegraph::HandleGraph;
use handlegraph::hashgraph::HashGraph;
use handlegraph::mutablehandlegraph::MutableHandleGraph;
use handlegraph::pathhandlegraph::embedded_paths::{GraphPaths, MutableGraphPaths};

fn main() {
    println!("Testing handlegraph API...");

    let mut graph = HashGraph::new();

    // Create some nodes
    let n1 = graph.create_handle(b"A", NodeId::from(1));
    let n2 = graph.create_handle(b"T", NodeId::from(2));
    let n3 = graph.create_handle(b"G", NodeId::from(3));

    println!("Created nodes: {:?}, {:?}, {:?}", n1, n2, n3);

    // Create edges
    graph.create_edge(&Edge(n1, n2));
    graph.create_edge(&Edge(n2, n3));

    println!("Created edges");

    // Create a path
    let path_name = b"test_path";
    match graph.create_path(path_name, false) {
        Ok(path_id) => {
            println!("Created path: {:?}", path_id);

            // Append steps
            graph.append_step(&path_id, n1);
            graph.append_step(&path_id, n2);
            graph.append_step(&path_id, n3);

            println!("Added steps to path");
        }
        Err(e) => println!("Failed to create path: {:?}", e),
    }

    // Test graph properties
    println!("Node count: {}", graph.node_count());
    println!("Edge count: {}", graph.edge_count());
    println!("Path count: {}", graph.path_count());

    // Test iteration
    println!("\nNodes:");
    for node_id in graph.handles_iter() {
        let handle = graph.handle_from_id(node_id);
        let seq = graph.sequence_vec(handle);
        println!(
            "  Node {}: {}",
            node_id.as_integer(),
            String::from_utf8_lossy(&seq)
        );
    }

    println!("\nEdges:");
    for edge in graph.edges_iter() {
        println!(
            "  {} -> {}",
            edge.0.id().as_integer(),
            edge.1.id().as_integer()
        );
    }
}
