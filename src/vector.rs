/*
 * Copyright (c) 2018-2020 Thomas Kramer.
 *
 * This file is part of LibrEDA 
 * (see https://codeberg.org/libreda).
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
use std::fmt;
use std::ops::{Add, AddAssign, Sub, SubAssign, Neg, Mul, MulAssign, Div};

use num_traits::{Float, NumCast};
pub use num_traits::Zero;

pub use crate::CoordinateType;
pub use crate::types::Orientation;
pub use crate::traits::{Angle, Mirror, RotateOrtho};

/// [`Vector`] defines a two dimensional vector with x and y components in the Euclidean plane.
#[derive(Copy, Clone, Default, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Vector<T> {
    pub x: T,
    pub y: T,
}

/// Shorthand notation for creating a vector.
macro_rules! vector {
 ($x:expr, $y:expr) => {Vector::new($x, $y)}
}

impl<T: CoordinateType> Into<(T, T)> for Vector<T> {
    fn into(self) -> (T, T) {
        (self.x, self.y)
    }
}

impl<T: CoordinateType> Into<(T, T)> for &Vector<T> {
    fn into(self) -> (T, T) {
        (self.x, self.y)
    }
}

impl<T: CoordinateType> From<(T, T)> for Vector<T> {
    fn from(coords: (T, T)) -> Self {
        Vector::new(coords.0, coords.1)
    }
}

impl<'a, T: CoordinateType> From<&'a (T, T)> for Vector<T> {
    fn from(coords: &'a (T, T)) -> Self {
        Vector::new(coords.0, coords.1)
    }
}

impl<'a, T: CoordinateType> From<&'a Vector<T>> for Vector<T> {
    fn from(v: &'a Vector<T>) -> Self {
        Vector::new(v.x, v.y)
    }
}

impl<T: CoordinateType> From<[T; 2]> for Vector<T> {
    fn from(coords: [T; 2]) -> Self {
        Vector::new(coords[0], coords[1])
    }
}

impl<T> fmt::Debug for Vector<T>
    where T: fmt::Debug + Copy {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({:?}, {:?})", self.x, self.y)
    }
}

impl<T> fmt::Display for Vector<T>
    where T: fmt::Display + Copy
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({}, {})", self.x, self.y)
    }
}

impl<T: CoordinateType> Zero for Vector<T> {
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
    fn zero() -> Self { Vector { x: T::zero(), y: T::zero() } }

    /// Check if this is the zero-vector.
    ///
    /// # Examples
    /// ```
    /// use iron_shapes::vector::{Vector, Zero};
    ///
    /// assert!(Vector::<usize>::zero().is_zero());
    /// ```
    fn is_zero(&self) -> bool {
        self.x.is_zero() && self.y.is_zero()
    }
}

impl<T: CoordinateType> Vector<T> {
    pub fn new(x: T, y: T) -> Self {
        Vector { x, y }
    }

    /// Get squared 2-norm of vector.
    ///
    /// # Examples
    /// ```
    /// use iron_shapes::vector::Vector;
    /// let a = Vector::new(2, 3);
    /// assert_eq!(a.norm2_squared(), 2*2+3*3);
    /// ```
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
    pub fn cross_prod(&self, other: Self) -> T {
        self.x * other.y - other.x * self.y
    }

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

impl<T: CoordinateType + NumCast> Vector<T> {
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
    pub fn try_cast<Target>(&self) -> Option<Vector<Target>>
        where Target: CoordinateType + NumCast {
        match (Target::from(self.x), Target::from(self.y)) {
            (Some(x), Some(y)) => Some(Vector { x, y }),
            _ => None
        }
    }

    /// Cast to vector of target data type.
    ///
    /// Conversion from float to int can fail and panic because float values
    /// like infinity or non-a-number have no integer representation.
    ///
    /// # Panics
    /// Panics if casting of the coordinate values fails. For instance when trying to cast a 'NaN' float into a integer.
    /// Use `try_cast` to detect and handle this failure.
    ///
    /// # Examples
    ///
    /// ```
    /// use iron_shapes::vector::Vector;
    ///
    /// let v_int = Vector::new(1,2);
    /// let v_float: Vector<f64> = v_int.cast();
    ///
    /// assert_eq!(v_float, Vector::new(1.0, 2.0));
    /// ```
    pub fn cast<S>(&self) -> Vector<S>
        where S: CoordinateType + NumCast {
        self.try_cast().unwrap()
    }

    /// Convert vector into a vector with floating point data type.
    pub fn cast_to_float<F: CoordinateType + Float + NumCast>(&self) -> Vector<F> {
        // TODO: find conversion that does not panic for sure.
        Vector {
            x: F::from(self.x).unwrap(),
            y: F::from(self.y).unwrap(),
        }
    }
}

impl<T: CoordinateType + Float> Vector<T>
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
    pub fn norm2(&self) -> T {
        (self.x * self.x + self.y * self.y).sqrt()
    }

    /// Return a vector with the same direction but length 1.
    ///
    /// # Panics
    /// Panics if the vector has length 0.
    pub fn normalized(&self) -> Self {
        *self / self.norm2()
    }

    /// Return the normal vector onto this vector.
    /// The normal has length 0.
    ///
    /// # Panics
    /// Panics if the vector has length 0.
    pub fn normal(&self) -> Self {
        self.normalized().rotate_ortho(Angle::R90)
    }
}

impl<T: CoordinateType + NumCast> Vector<T>
{
    /// Calculate length of vector.
    pub fn length<F: Float>(&self) -> F {
        let x = F::from(self.x).unwrap();
        let y = F::from(self.y).unwrap();

        x.hypot(y)
    }
}

/// Vector addition.
impl<T> Add for Vector<T>
    where T: CoordinateType + Add<Output=T>
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Vector {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

impl<T> AddAssign for Vector<T>
    where T: CoordinateType + AddAssign
{
    fn add_assign(&mut self, rhs: Self) {
        self.x += rhs.x;
        self.y += rhs.y;
    }
}

/// Vector subtraction.
impl<T> Sub for Vector<T>
    where T: CoordinateType + Sub<Output=T>
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Vector {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
        }
    }
}

impl<T> SubAssign for Vector<T>
    where T: CoordinateType + SubAssign
{
    fn sub_assign(&mut self, rhs: Self) {
        self.x -= rhs.x;
        self.y -= rhs.y;
    }
}

impl<T> Neg for Vector<T>
    where T: CoordinateType
{
    type Output = Self;

    fn neg(self) -> Self {
        Vector {
            x: T::zero() - self.x,
            y: T::zero() - self.y,
        }
    }
}

/// Scalar multiplication.
impl<T> Mul<T> for Vector<T>
    where T: CoordinateType + Mul<Output=T>
{
    type Output = Self;

    fn mul(self, rhs: T) -> Self {
        Vector {
            x: self.x * rhs,
            y: self.y * rhs,
        }
    }
}

/// In-place scalar multiplication.
impl<T> MulAssign<T> for Vector<T>
    where T: CoordinateType + MulAssign<T>
{
    fn mul_assign(&mut self, rhs: T) {
        self.x *= rhs;
        self.y *= rhs;
    }
}

/// Scalar division.
impl<T> Div<T> for Vector<T>
    where T: CoordinateType + Div<Output=T>
{
    type Output = Self;

    fn div(self, rhs: T) -> Self {
        Vector {
            x: self.x / rhs,
            y: self.y / rhs,
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use num_rational::Rational;
    use num_traits::Zero;

    #[test]
    fn test_rational_vector() {
        let _ = Vector::new(Rational::zero(), Rational::zero());
    }
}
