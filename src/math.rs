// Copyright (c) 2018-2020 Thomas Kramer.
// SPDX-FileCopyrightText: 2018-2022 Thomas Kramer
//
// SPDX-License-Identifier: AGPL-3.0-or-later

//! Math helper functions.

use num_traits::PrimInt;

/// Implement the fast approximate computation of '1/sqrt(x)' for a type.
pub trait FastInvSqrt {
    /// Fast approximate computation of 1/sqrt(x).
    ///
    /// See: <https://en.wikipedia.org/wiki/Fast_inverse_square_root>
    /// And: <http://www.lomont.org/papers/2003/InvSqrt.pdf>
    fn fast_invsqrt(x: Self) -> Self;
}

impl FastInvSqrt for f32 {
    #[inline]
    fn fast_invsqrt(x: Self) -> Self {
        fast_invsqrt_32(x)
    }
}

impl FastInvSqrt for f64 {
    #[inline]
    fn fast_invsqrt(x: Self) -> Self {
        fast_invsqrt_64(x)
    }
}


/// Fast approximate computation of 1/sqrt(x).
/// The error should be below 0.2%.
///
/// See: <https://en.wikipedia.org/wiki/Fast_inverse_square_root>
/// And: <http://www.lomont.org/papers/2003/InvSqrt.pdf>
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
/// See: <https://en.wikipedia.org/wiki/Fast_inverse_square_root>
/// And: <http://www.lomont.org/papers/2003/InvSqrt.pdf>
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
        assert!(abs_diff / inv_sqrt < 0.002, "Error should be below 0.2%.")
    }
}

#[test]
fn test_fast_invsqrt_64() {
    for i in 1..1000000 {
        let x = i as f64;
        let inv_sqrt_approx = fast_invsqrt_64(x);
        let inv_sqrt = 1. / x.sqrt();
        let abs_diff = (inv_sqrt - inv_sqrt_approx).abs();
        assert!(abs_diff / inv_sqrt < 0.002, "Error should be below 0.2%.")
    }
}

/// Compute square root of integers using Newtons method.
/// Returns the biggest integer which is smaller or equal to the actual square root of `n`.
/// Similar to `(i as f64).sqrt().floor() as T` but without type conversions.
///
/// See: <https://en.wikipedia.org/wiki/Integer_square_root>
///
/// # Panics
/// Panics when given a negative number.
///
/// # Example
/// ```
/// use iron_shapes::math::int_sqrt_floor;
/// assert_eq!(int_sqrt_floor(16), 4);
/// assert_eq!(int_sqrt_floor(17), 4);
/// assert_eq!(int_sqrt_floor(24), 4);
/// assert_eq!(int_sqrt_floor(25), 5);
/// ```
pub fn int_sqrt_floor<T: PrimInt>(n: T) -> T {
    assert!(n >= T::zero(), "Cannot compute the square root of a negative number.");
    let _1 = T::one();
    let _2 = _1 + _1;
    let mut x = n;
    let mut y = (x + _1) / _2;
    while y < x {
        x = y;
        y = (x + (n / x)) / _2;
    }
    x
}


#[test]
pub fn test_int_sqrt_floor_i32() {
    for i in 0..100000i32 {
        let sqrt = int_sqrt_floor(i);
        let sqrt_reference = (i as f64).sqrt().floor() as i32;
        assert_eq!(sqrt, sqrt_reference);

        let i_lower = sqrt * sqrt;
        let i_upper = (sqrt + 1) * (sqrt + 1);
        assert!(i_lower <= i && i < i_upper);
    }
}
