use handlegraph::handle::Handle;
fn main() {
    for node_id in 1..=5 {
        let fwd = Handle::forward(node_id);
        let rev = Handle::reverse(node_id);
        println!("Node {}: fwd={}, rev={}", node_id, fwd.as_integer(), rev.as_integer());
    }
}
