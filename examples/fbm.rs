//! An example of using the fBm noise function
use noice::{utils::*, Fbm};

fn main() {
    let fbm = Fbm::new();

    PlaneMapBuilder::new(&fbm).build().write_to_file("fbm.png");
}
