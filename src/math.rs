/*
 * Copyright (c) 2018-2020 Thomas Kramer.
 *
 * This file is part of LibrEDA
 * (see https://codeberg.org/libreda/iron-shapes).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

//! Math helper functions.

use std::io::BufRead;
use num_traits::Num;

pub trait FastInvSqrt {
    /// Fast approximate computation of 1/sqrt(x).
    ///
    /// See: https://en.wikipedia.org/wiki/Fast_inverse_square_root
    /// And: http://www.lomont.org/papers/2003/InvSqrt.pdf
    #[inline]
    fn fast_invsqrt(x: Self) -> Self;
}

impl FastInvSqrt for f32 {
    fn fast_invsqrt(x: Self) -> Self {
        fast_invsqrt_32(x)
    }
}

impl FastInvSqrt for f64 {
    fn fast_invsqrt(x: Self) -> Self {
        fast_invsqrt_64(x)
    }
}

/// Fast approximate computation of 1/sqrt(x).
/// The error should be below 0.2%.
///
/// See: https://en.wikipedia.org/wiki/Fast_inverse_square_root
/// And: http://www.lomont.org/papers/2003/InvSqrt.pdf
pub fn fast_invsqrt_32(x: f32) -> f32 {
    let approx = 0x5f375a86u32;
    let xhalf = 0.5 * x;
    let i = x.to_bits();
    let i = approx - (i >> 1);
    let y = f32::from_bits(i);
    y * (1.5 - xhalf * y * y)
}

/// Fast approximate computation of 1/sqrt(x).
/// The error should be below 0.2%.
///
/// See: https://en.wikipedia.org/wiki/Fast_inverse_square_root
/// And: http://www.lomont.org/papers/2003/InvSqrt.pdf
pub fn fast_invsqrt_64(x: f64) -> f64 {
    let approx = 0x5fe6eb50c7aa19f9u64;
    let xhalf = 0.5 * x;
    let i = x.to_bits();
    let i = approx - (i >> 1);
    let y = f64::from_bits(i);
    y * (1.5 - xhalf * y * y)
}

#[test]
fn test_fast_invsqrt_32() {
    for i in 1..1000000 {
        let x = i as f32;
        let inv_sqrt_approx = fast_invsqrt_32(x);
        let inv_sqrt = 1. / x.sqrt();
        let abs_diff = (inv_sqrt - inv_sqrt_approx).abs();
        assert!(abs_diff/inv_sqrt < 0.002, "Error should be below 0.2%.")
    }
}

#[test]
fn test_fast_invsqrt_64() {
    for i in 1..1000000 {
        let x = i as f64;
        let inv_sqrt_approx = fast_invsqrt_64(x);
        let inv_sqrt = 1. / x.sqrt();
        let abs_diff = (inv_sqrt - inv_sqrt_approx).abs();
        assert!(abs_diff/inv_sqrt < 0.002, "Error should be below 0.2%.")
    }
}