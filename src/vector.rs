// Copyright (c) 2018-2020 Thomas Kramer.
// SPDX-FileCopyrightText: 2018-2022 Thomas Kramer
//
// SPDX-License-Identifier: AGPL-3.0-or-later

//! Two dimensional vectors are a core data type for Euclidean geometry in the plane.
//! `Vector`s consist of an `x` and `y` coordinate value and describe a translation in the plane.

use std::fmt;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

pub use num_traits::Zero;
use num_traits::{Float, NumCast};

use crate::point::Point;
use crate::traits::{MapPointwise, TryCastCoord};
pub use crate::traits::{Mirror, RotateOrtho};
pub use crate::types::Angle;
pub use crate::types::Orientation;
pub use crate::CoordinateType;

/// [`Vector`] defines a two dimensional vector with x and y components in the Euclidean plane.
#[derive(Copy, Clone, Default, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Vector<T> {
    /// `x` coordinate.
    pub x: T,
    /// `y` coordinate.
    pub y: T,
}

/// Shorthand notation for creating a vector.
///
/// # Example
/// ```
/// # #[macro_use]
/// # extern crate iron_shapes;
/// # fn main() {
/// use iron_shapes::prelude::*;
/// let v = vector!(1, 2);
/// assert_eq!(v, Vector::new(1, 2));
/// # }
/// ```
#[macro_export]
macro_rules! vector {
    ($x:expr, $y:expr) => {
        Vector::new($x, $y)
    };
}

impl<T: Copy> Into<(T, T)> for Vector<T> {
    #[inline]
    fn into(self) -> (T, T) {
        (self.x, self.y)
    }
}

impl<T: Copy> Into<(T, T)> for &Vector<T> {
    #[inline]
    fn into(self) -> (T, T) {
        (self.x, self.y)
    }
}

impl<T: Copy> From<(T, T)> for Vector<T> {
    #[inline]
    fn from(coords: (T, T)) -> Self {
        Vector::new(coords.0, coords.1)
    }
}

impl<'a, T: Copy> From<&'a (T, T)> for Vector<T> {
    #[inline]
    fn from(coords: &'a (T, T)) -> Self {
        Vector::new(coords.0, coords.1)
    }
}

impl<'a, T: Copy> From<&'a Vector<T>> for Vector<T> {
    #[inline]
    fn from(v: &'a Vector<T>) -> Self {
        Vector::new(v.x, v.y)
    }
}

impl<T: Copy> From<[T; 2]> for Vector<T> {
    #[inline]
    fn from(coords: [T; 2]) -> Self {
        Vector::new(coords[0], coords[1])
    }
}

impl<T> fmt::Debug for Vector<T>
where
    T: fmt::Debug,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({:?}, {:?})", &self.x, &self.y)
    }
}

impl<T> fmt::Display for Vector<T>
where
    T: fmt::Display + Copy,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({}, {})", self.x, self.y)
    }
}

impl<T: Zero> Zero for Vector<T> {
    /// Get zero-vector.
    ///
    /// # Examples
    /// ```
    /// use iron_shapes::vector::{Vector, Zero};
    ///
    /// let a = Vector::zero();
    /// let b = Vector::new(0, 0);
    ///
    /// assert_eq!(a, b);
    /// ```
    #[inline]
    fn zero() -> Self {
        Vector {
            x: T::zero(),
            y: T::zero(),
        }
    }

    /// Check if this is the zero-vector.
    ///
    /// # Examples
    /// ```
    /// use iron_shapes::vector::{Vector, Zero};
    ///
    /// assert!(Vector::<usize>::zero().is_zero());
    /// ```
    #[inline]
    fn is_zero(&self) -> bool {
        self.x.is_zero() && self.y.is_zero()
    }
}

impl<T> Vector<T> {
    /// Create a new vector with `x` and `y` coordinates.
    /// # Examples
    /// ```
    /// use iron_shapes::vector::Vector;
    /// let a = Vector::new(2, 3);
    /// assert_eq!(a.x, 2);
    /// assert_eq!(a.y, 3);
    /// ```
    pub fn new(x: T, y: T) -> Self {
        Vector { x, y }
    }
}

impl<T: Copy + Zero + PartialOrd + Sub<Output = T>> Vector<T> {
    /// Get 1-norm of vector, i.e. the sum of the absolute values of its components.
    ///
    /// # Examples
    /// ```
    /// use iron_shapes::vector::Vector;
    /// let a = Vector::new(-2, 3);
    /// assert_eq!(a.norm1(), 5);
    /// ```
    #[inline]
    pub fn norm1(&self) -> T {
        let xabs = if self.x < T::zero() {
            T::zero() - self.x
        } else {
            self.x
        };
        let yabs = if self.y < T::zero() {
            T::zero() - self.y
        } else {
            self.y
        };
        xabs + yabs
    }
}

impl<T: Copy + Zero + PartialOrd + Mul<Output = T> + Sub<Output = T>> Vector<T> {
    /// Check if `other` is oriented clockwise or counter-clockwise respective to `self`.
    ///
    /// # Examples
    ///
    /// ```
    /// use iron_shapes::vector::Vector;
    /// use iron_shapes::types::Orientation;
    ///
    /// let a = Vector::new(1, 0);
    /// let b = Vector::new(1, 1);
    /// let c = Vector::new(1, -1);
    /// let d = Vector::new(2, 0);
    ///
    /// assert_eq!(a.orientation_of(b), Orientation::CounterClockWise);
    /// assert_eq!(a.orientation_of(c), Orientation::ClockWise);
    /// assert_eq!(a.orientation_of(d), Orientation::Straight);
    /// ```
    pub fn orientation_of(&self, other: Self) -> Orientation {
        let p = self.cross_prod(other);

        if p < T::zero() {
            Orientation::ClockWise
        } else if p > T::zero() {
            Orientation::CounterClockWise
        } else {
            Orientation::Straight
        }
    }
}

impl<T: Copy + Mul<Output = T> + Add<Output = T>> Vector<T> {
    /// Get squared 2-norm of vector.
    ///
    /// # Examples
    /// ```
    /// use iron_shapes::vector::Vector;
    /// let a = Vector::new(2, 3);
    /// assert_eq!(a.norm2_squared(), 2*2+3*3);
    /// ```
    #[inline]
    pub fn norm2_squared(&self) -> T {
        self.x * self.x + self.y * self.y
    }

    /// Calculate scalar product.
    ///
    /// # Examples
    ///
    /// ```
    /// use iron_shapes::vector::Vector;
    ///
    /// let a = Vector::new(1, 2);
    /// let b = Vector::new(3, 4);
    ///
    /// assert_eq!(a.dot(b), 1*3 + 2*4);
    /// ```
    pub fn dot(&self, other: Self) -> T {
        self.x * other.x + self.y * other.y
    }
}

impl<T: Copy + Mul<Output = T> + Sub<Output = T>> Vector<T> {
    /// Calculate cross product.
    ///
    /// # Examples
    ///
    /// ```
    /// use iron_shapes::vector::Vector;
    ///
    /// let a = Vector::new(2, 0);
    /// let b = Vector::new(0, 2);
    ///
    /// assert_eq!(a.cross_prod(b), 4);
    /// assert_eq!(b.cross_prod(a), -4);
    /// ```
    #[inline]
    pub fn cross_prod(&self, other: Self) -> T {
        self.x * other.y - other.x * self.y
    }
}

impl<T: CoordinateType + NumCast> Vector<T> {
    /// Convert vector into a vector with floating point data type.
    #[inline]
    pub fn cast_to_float<F: CoordinateType + Float + NumCast>(&self) -> Vector<F> {
        // TODO: find conversion that does not panic for sure.
        Vector {
            x: F::from(self.x).unwrap(),
            y: F::from(self.y).unwrap(),
        }
    }
}

impl<T: Copy + NumCast, Dst: Copy + NumCast> TryCastCoord<T, Dst> for Vector<T> {
    type Output = Vector<Dst>;

    /// Try to cast to vector of target data type.
    ///
    /// Conversion from float to int can fail and will return `None`.
    /// Float values like infinity or non-a-number
    /// have no integer representation.
    ///
    /// # Examples
    ///
    /// ```
    /// use iron_shapes::vector::Vector;
    /// use iron_shapes::traits::TryCastCoord;
    ///
    /// let v_int = Vector::new(1,2);
    /// let maybe_v_float: Option<Vector<f64>> = v_int.try_cast();
    ///
    /// assert_eq!(maybe_v_float, Some(Vector::new(1.0, 2.0)));
    ///
    /// // Conversion from float to int can fail.
    ///
    /// let w_float = Vector::new(42.0, 0. / 0.);
    /// let maybe_w_int: Option<Vector<i32>> = w_float.try_cast();
    ///
    /// assert_eq!(maybe_w_int, None);
    /// ```
    #[inline]
    fn try_cast(&self) -> Option<Self::Output> {
        match (Dst::from(self.x), Dst::from(self.y)) {
            (Some(x), Some(y)) => Some(Vector::new(x, y)),
            _ => None,
        }
    }
}

impl<T> Vector<T>
where
    T: Copy + Float + Mul<Output = T> + Add<Output = T>,
{
    /// Get 2-norm of vector (length of vector).
    ///
    /// # Examples
    /// ```
    /// use iron_shapes::vector::Vector;
    /// let a = Vector::new(2.0, 3.0);
    /// let norm2 = a.norm2();
    /// let norm2_sq = norm2 * norm2;
    /// let expected = a.norm2_squared();
    /// assert!(norm2_sq < expected + 1e-12);
    /// assert!(norm2_sq > expected - 1e-12);
    /// ```
    #[inline]
    pub fn norm2(&self) -> T {
        (self.norm2_squared()).sqrt()
    }

    /// Return a vector with the same direction but length 1.
    ///
    /// # Panics
    /// Panics if the vector has length 0.
    #[inline]
    pub fn normalized(&self) -> Self {
        *self / self.norm2()
    }

    /// Return the normal vector onto this vector.
    /// The normal has length `1`.
    ///
    /// # Panics
    /// Panics if the vector has length 0.
    #[inline]
    pub fn normal(&self) -> Self {
        self.normalized().rotate_ortho(Angle::R90)
    }
}

impl<T: Copy + NumCast> Vector<T> {
    /// Calculate length of vector.
    ///
    /// Similar to `Vector::norm2` but does potentially return another data type for the length.
    ///
    /// # Examples
    /// ```
    /// use iron_shapes::vector::Vector;
    /// let a = Vector::new(2.0, 3.0);
    /// let length: f64 = a.length();
    /// let norm2_sq = length * length;
    /// let expected = a.norm2_squared();
    /// assert!(norm2_sq < expected + 1e-12);
    /// assert!(norm2_sq > expected - 1e-12);
    /// ```
    #[inline]
    pub fn length<F: Float>(&self) -> F {
        let x = F::from(self.x).unwrap();
        let y = F::from(self.y).unwrap();

        x.hypot(y)
    }
}

/// Vector addition.
impl<T> Add for Vector<T>
where
    T: Add<Output = T>,
{
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self {
        Vector {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

impl<T> AddAssign for Vector<T>
where
    T: AddAssign,
{
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.x += rhs.x;
        self.y += rhs.y;
    }
}

/// Vector subtraction.
impl<T> Sub for Vector<T>
where
    T: Sub<Output = T>,
{
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Vector {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
        }
    }
}

impl<T> SubAssign for Vector<T>
where
    T: SubAssign,
{
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        self.x -= rhs.x;
        self.y -= rhs.y;
    }
}

impl<T> Neg for Vector<T>
where
    T: Neg<Output = T>,
{
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        Vector {
            x: -self.x,
            y: -self.y,
        }
    }
}

/// Scalar multiplication.
impl<T, M> Mul<M> for Vector<T>
where
    T: CoordinateType + Mul<M, Output = T>,
    M: Copy,
{
    type Output = Self;

    fn mul(self, rhs: M) -> Self {
        Vector {
            x: self.x * rhs,
            y: self.y * rhs,
        }
    }
}

/// In-place scalar multiplication.
impl<T, M> MulAssign<M> for Vector<T>
where
    T: Copy + MulAssign<M>,
    M: Copy,
{
    #[inline]
    fn mul_assign(&mut self, rhs: M) {
        self.x *= rhs;
        self.y *= rhs;
    }
}

/// Scalar division.
impl<T, D> Div<D> for Vector<T>
where
    T: Copy + Div<D, Output = T>,
    D: Copy,
{
    type Output = Self;

    #[inline]
    fn div(self, rhs: D) -> Self {
        Vector {
            x: self.x / rhs,
            y: self.y / rhs,
        }
    }
}

/// Assigning scalar division.
impl<T, D> DivAssign<D> for Vector<T>
where
    T: Copy + DivAssign<D>,
    D: Copy,
{
    #[inline]
    fn div_assign(&mut self, rhs: D) {
        self.x /= rhs;
        self.y /= rhs;
    }
}

impl<T: CoordinateType> MapPointwise<T> for Vector<T> {
    #[inline]
    fn transform<F>(&self, transformation: F) -> Self
    where
        F: Fn(Point<T>) -> Point<T>,
    {
        *transformation(Point::from(self))
    }
}

/// Compute the sum of all vectors in the iterator.
/// If the iterator is empty, (0, 0) is returned.
impl<T> std::iter::Sum for Vector<T>
where
    T: Zero + Add<Output = T>,
{
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, v| acc + v)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_rational::Rational64;
    use num_traits::Zero;

    #[test]
    fn test_rational_vector() {
        let _ = Vector::new(Rational64::zero(), Rational64::zero());
    }
}
