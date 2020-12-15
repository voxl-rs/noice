//! An example of generating constant valued noise
use noice::{utils::*, Checkerboard};

fn main() {
    let checker = Checkerboard::default();

    PlaneMapBuilder::new(&checker)
        .build()
        .write_to_file("checkerboard.png");
}
