use rhai::plugin::*;

#[export_module]
pub mod test_module {
    #[rhai_fn(name = "test")]
    #[rhai_fn(pure)]
    pub fn test_fn(input: Point) -> bool {
        input.x > input.y
    }
}

fn main() {
    if test_module::test_fn(n) {
        println!("yes");
    } else {
        println!("no");
    }
}
