//! An example of using the Billow noise function
use noice::{utils::*, Billow};

fn main() {
    PlaneMapBuilder::new(&Billow::new())
        .build()
        .write_to_file("billow.png");
}
