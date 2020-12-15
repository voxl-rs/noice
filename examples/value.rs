//! An example of using value noise
use noice::{utils::*, Value};

fn main() {
    PlaneMapBuilder::new(&Value::new())
        .build()
        .write_to_file("value.png");
}
