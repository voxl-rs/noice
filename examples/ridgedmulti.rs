//! An example of using the `RidgedMulti` noise function
use noice::{utils::*, RidgedMulti};

fn main() {
    let ridged_multi = RidgedMulti::new();

    PlaneMapBuilder::new(&ridged_multi)
        .build()
        .write_to_file("ridged_multi.png");
}
