use crate::{
    gradient, math,
    noise_fns::{NoiseFn, Seedable},
    permutationtable::PermutationTable,
};

/// Noise function that outputs 2/3/4-dimensional Simplex noise.
#[derive(Clone, Copy, Debug)]
pub struct Simplex {
    seed: u32,
    perm_table: PermutationTable,
}

impl Simplex {
    pub const DEFAULT_SEED: u32 = 0;

    pub fn new() -> Self {
        Simplex {
            seed: Self::DEFAULT_SEED,
            perm_table: PermutationTable::new(Self::DEFAULT_SEED),
        }
    }
}

impl Default for Simplex {
    fn default() -> Self {
        Self::new()
    }
}

impl Seedable for Simplex {
    /// Sets the seed value for Simplex noise
    fn set_seed(self, seed: u32) -> Self {
        // If the new seed is the same as the current seed, just return self.
        if self.seed == seed {
            return self;
        }

        // Otherwise, regenerate the permutation table based on the new seed.
        Simplex {
            seed,
            perm_table: PermutationTable::new(seed),
        }
    }

    fn seed(&self) -> u32 {
        self.seed
    }
}

// Skew Value
//
//     sqrt(n + 1) - 1
// F = ---------------
//            n
pub fn skew_factor(n: usize) -> f64 {
    let n = n as f64;

    ((n + 1.0).sqrt() - 1.0) / n
}

//  Unskew Value
//
//     1 - 1 / sqrt(n + 1)
// G = -------------------
//             n
pub fn unskew_factor(n: usize) -> f64 {
    let n = n as f64;

    (1.0 - (1.0 / (n + 1.0).sqrt())) / n
}

/// 2-dimensional Simplex noise
impl NoiseFn<[f64; 2]> for Simplex {
    fn get(&self, point: [f64; 2]) -> f64 {
        #[inline]
        fn surflet(gradient_index: usize, distance: [f64; 2]) -> f64 {
            let mut t = 0.5 - distance[0] * distance[0] - distance[1] * distance[1];

            if t < 0.0 {
                0.0
            } else {
                t *= t;
                t * t * math::dot2(gradient::get2(gradient_index), distance)
            }
        }

        /// Skew the input point per the following formula:
        /// x' = x + (x + y) * F
        /// y' = y + (x + y) * F
        fn skew_point(point: [f64; 2], factor: f64) -> [f64; 2] {
            math::add2(point, [(point[0] + point[1]) * factor; 2])
        }

        /// Unskew the input point per the following formula:
        /// x = x' - (x' + y') * G
        /// y = y' - (x' + y`) * G
        fn unskew_point(skewed_point: [f64; 2], factor: f64) -> [f64; 2] {
            math::sub2(
                skewed_point,
                [(skewed_point[0] + skewed_point[1]) * factor; 2],
            )
        }

        let skew = skew_factor(2);
        let unskew = unskew_factor(2);

        let skewed_input = skew_point(point, skew);

        // Floor the skewed coordinate to determine which skewed unit cell the point is in.
        let floored = math::to_isize2(math::map2(skewed_input, f64::floor));

        let cell = unskew_point(math::to_f64_2(floored), unskew);

        // Calculate the vector from the cell's minimum corner to the point.
        let distance = math::sub2(point, cell);

        // Sort the coordinates in decreasing order based on the largest component of the distance
        // vector.
        let offsets = if distance[0] > distance[1] {
            [1, 0]
        } else {
            [0, 1]
        };

        let corner2 = math::add2(math::sub2(distance, math::to_f64_2(offsets)), [unskew; 2]);

        let corner3 = math::add2(math::sub2(distance, [1.0; 2]), [2.0 * unskew; 2]);

        let gi0 = self.perm_table.get2(floored);
        let gi1 = self.perm_table.get2(math::add2(floored, offsets));
        let gi2 = self.perm_table.get2(math::add2(floored, [1; 2]));

        let n0 = surflet(gi0, distance);
        let n1 = surflet(gi1, corner2);
        let n2 = surflet(gi2, corner3);

        70.0 * (n0 + n1 + n2)
    }
}

/// 3-dimensional Simplex noise
impl NoiseFn<[f64; 3]> for Simplex {
    fn get(&self, point: [f64; 3]) -> f64 {
        #[inline]
        fn surflet(gradient_index: usize, distance: [f64; 3]) -> f64 {
            let mut t = 0.5
                - distance[0] * distance[0]
                - distance[1] * distance[1]
                - distance[2] * distance[2];

            if t < 0.0 {
                0.0
            } else {
                t *= t;
                t * t * math::dot3(gradient::get3(gradient_index), distance)
            }
        }

        /// Skew the input point per the following formula:
        /// x' = x + (x + y + ...) * F
        /// y' = y + (x + y + ...) * F
        /// :
        fn skew_point(point: [f64; 3], factor: f64) -> [f64; 3] {
            math::add3(point, [(point[0] + point[1] + point[2]) * factor; 3])
        }

        /// Unskew the input point per the following formula:
        /// x = x' - (x' + y' + ...) * G
        /// y = y' - (x' + y` + ...) * G
        /// :
        fn unskew_point(skewed_point: [f64; 3], factor: f64) -> [f64; 3] {
            math::sub3(
                skewed_point,
                [(skewed_point[0] + skewed_point[1] + skewed_point[2]) * factor; 3],
            )
        }

        // Skew Value
        let skew = skew_factor(3);
        // Unskew value
        let unskew = unskew_factor(3);

        let skewed_input = skew_point(point, skew);

        let floored = math::to_isize3(math::map3(skewed_input, f64::floor));

        let cell = unskew_point(math::to_f64_3(floored), unskew);

        let distance = math::sub3(point, cell);

        let offset1;
        let offset2;

        if distance[0] >= distance[1] {
            if distance[1] >= distance[2] {
                offset1 = [1, 0, 0];
                offset2 = [1, 1, 0];
            } else if distance[0] >= distance[2] {
                offset1 = [1, 0, 0];
                offset2 = [1, 0, 1];
            } else {
                offset1 = [0, 0, 1];
                offset2 = [1, 0, 1];
            }
        } else if distance[2] >= distance[1] {
            offset1 = [0, 0, 1];
            offset2 = [0, 1, 1];
        } else if distance[2] >= distance[0] {
            offset1 = [0, 1, 0];
            offset2 = [0, 1, 1];
        } else {
            offset1 = [0, 1, 0];
            offset2 = [1, 1, 0];
        }

        let offset3 = [1; 3];

        let corner2 = math::add3(math::sub3(distance, math::to_f64_3(offset1)), [unskew; 3]);

        let corner3 = math::add3(
            math::sub3(distance, math::to_f64_3(offset2)),
            [2.0 * unskew; 3],
        );

        let corner4 = math::add3(math::sub3(distance, [1.0; 3]), [3.0 * unskew; 3]);

        let gi0 = self.perm_table.get3(floored);
        let gi1 = self.perm_table.get3(math::add3(floored, offset1));
        let gi2 = self.perm_table.get3(math::add3(floored, offset2));
        let gi3 = self.perm_table.get3(math::add3(floored, offset3));

        let n0 = surflet(gi0, distance);
        let n1 = surflet(gi1, corner2);
        let n2 = surflet(gi2, corner3);
        let n3 = surflet(gi3, corner4);

        32.0 * (n0 + n1 + n2 + n3)
    }
}

/// 4-dimensional Simplex noise
impl NoiseFn<[f64; 4]> for Simplex {
    fn get(&self, point: [f64; 4]) -> f64 {
        #[inline]
        fn surflet(gradient_index: usize, distance: [f64; 4]) -> f64 {
            let mut t = 0.5
                - distance[0] * distance[0]
                - distance[1] * distance[1]
                - distance[2] * distance[2]
                - distance[3] * distance[3];

            if t < 0.0 {
                0.0
            } else {
                t *= t;
                t * t * math::dot4(gradient::get4(gradient_index), distance)
            }
        }

        /// Skew the input point per the following formula:
        /// x' = x + (x + y + ...) * F
        /// y' = y + (x + y + ...) * F
        /// :
        fn skew_point(point: [f64; 4], factor: f64) -> [f64; 4] {
            math::add4(
                point,
                [(point[0] + point[1] + point[2] + point[3]) * factor; 4],
            )
        }

        /// Unskew the input point per the following formula:
        /// x = x' - (x' + y' + ...) * G
        /// y = y' - (x' + y` + ...) * G
        /// :
        fn unskew_point(skewed_point: [f64; 4], factor: f64) -> [f64; 4] {
            math::sub4(
                skewed_point,
                [(skewed_point[0] + skewed_point[1] + skewed_point[2] + skewed_point[3]) * factor;
                    4],
            )
        }

        // Skew Value
        let skew: f64 = skew_factor(4);
        // Unskew Value
        let unskew: f64 = unskew_factor(4);

        let skewed_input = skew_point(point, skew);

        let floored = math::to_isize4(math::map4(skewed_input, f64::floor));

        let cell = unskew_point(math::to_f64_4(floored), unskew);

        let distance = math::sub4(point, cell);

        let mut rank_x: u8 = 0;
        let mut rank_y: u8 = 0;
        let mut rank_z: u8 = 0;
        let mut rank_w: u8 = 0;

        if distance[0] > distance[1] {
            rank_x += 1;
        } else {
            rank_y += 1;
        };
        if distance[0] > distance[2] {
            rank_x += 1;
        } else {
            rank_z += 1;
        };
        if distance[0] > distance[3] {
            rank_x += 1;
        } else {
            rank_w += 1;
        };
        if distance[1] > distance[2] {
            rank_y += 1;
        } else {
            rank_z += 1;
        };
        if distance[1] > distance[3] {
            rank_y += 1;
        } else {
            rank_w += 1;
        };
        if distance[2] > distance[3] {
            rank_z += 1;
        } else {
            rank_w += 1;
        };

        let mut offset1 = [0; 4];
        let mut offset2 = [0; 4];
        let mut offset3 = [0; 4];

        if rank_x >= 3 {
            offset1[0] = 1
        };
        if rank_y >= 3 {
            offset1[1] = 1
        };
        if rank_z >= 3 {
            offset1[2] = 1
        };
        if rank_w >= 3 {
            offset1[3] = 1
        };

        if rank_x >= 2 {
            offset2[0] = 1
        };
        if rank_y >= 2 {
            offset2[1] = 1
        };
        if rank_z >= 2 {
            offset2[2] = 1
        };
        if rank_w >= 2 {
            offset2[3] = 1
        };

        if rank_x >= 1 {
            offset3[0] = 1
        };
        if rank_y >= 1 {
            offset3[1] = 1
        };
        if rank_z >= 1 {
            offset3[2] = 1
        };
        if rank_w >= 1 {
            offset3[3] = 1
        };

        let offset4 = [1; 4];

        let corner2 = math::add4(math::sub4(distance, math::to_f64_4(offset1)), [unskew; 4]);

        let corner3 = math::add4(
            math::sub4(distance, math::to_f64_4(offset2)),
            [2.0 * unskew; 4],
        );

        let corner4 = math::add4(
            math::sub4(distance, math::to_f64_4(offset3)),
            [3.0 * unskew; 4],
        );

        let corner5 = math::add4(math::sub4(distance, [1.0; 4]), [4.0 * unskew; 4]);

        let gi0 = self.perm_table.get4(floored);
        let gi1 = self.perm_table.get4(math::add4(floored, offset1));
        let gi2 = self.perm_table.get4(math::add4(floored, offset2));
        let gi3 = self.perm_table.get4(math::add4(floored, offset3));
        let gi4 = self.perm_table.get4(math::add4(floored, offset4));

        let n0 = surflet(gi0, distance);
        let n1 = surflet(gi1, corner2);
        let n2 = surflet(gi2, corner3);
        let n3 = surflet(gi3, corner4);
        let n4 = surflet(gi4, corner5);

        27.0 * (n0 + n1 + n2 + n3 + n4)
    }
}
