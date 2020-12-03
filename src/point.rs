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
use crate::vector::Vector;
use crate::traits::{MapPointwise, TryCastCoord};
use crate::CoordinateType;


use std::cmp::{Ord, Ordering};
use std::fmt;

use num_traits::{Float, NumCast};
pub use num_traits::Zero;
use std::ops::{Div, MulAssign, Mul, Neg, Sub, SubAssign, AddAssign, Add};
pub use std::ops::Deref;

/// Shorthand notation for creating a point.
///
/// # Example
/// ```
/// # #[macro_use]
/// # extern crate iron_shapes;
/// # fn main() {
/// use iron_shapes::prelude::*;
/// let p = point!(1, 2);
/// assert_eq!(p, Point::new(1, 2));
/// # }
/// ```
#[macro_export]
macro_rules! point {
 ($x:expr, $y:expr) => {Point::new($x, $y)}
}

/// A point is defined by a x and y coordinate in the euclidean plane.
#[derive(Copy, Clone, Hash, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Point<T>
    where T: CoordinateType {
    location: Vector<T>
}

impl<T> Deref for Point<T>
    where T: CoordinateType {
    type Target = Vector<T>;

    fn deref(&self) -> &Self::Target {
        &self.location
    }
}

impl<T: CoordinateType> From<Vector<T>> for Point<T> {
    fn from(v: Vector<T>) -> Self {
        Point {
            location: v
        }
    }
}

impl<T: CoordinateType> From<&Vector<T>> for Point<T> {
    fn from(v: &Vector<T>) -> Self {
        Self::from(*v)
    }
}

impl<T: CoordinateType> Into<Vector<T>> for Point<T> {
    fn into(self) -> Vector<T> {
        self.location
    }
}

///// Convert a type into a point by converting it into a vector first.
//impl<S, T: CoordinateType> From<S> for Point<T>
//    where S: Into<Vector<T>> {
//    fn from(s: S) -> Self {
//        let v: Vector<T> = s.into();
//        v.into()
//    }
//}

impl<T: CoordinateType> Point<T> {
    /// Create a new point with `x` and `y` coordinates.
    pub fn new(x: T, y: T) -> Self {
        Point { location: Vector::new(x, y) }
    }

    /// Get zero-Point.
    ///
    /// # Examples
    /// ```
    /// use iron_shapes::point::Point;
    ///
    /// let a = Point::zero();
    /// let b = Point::new(0, 0);
    ///
    /// assert_eq!(a, b);
    /// ```
    pub fn zero() -> Self { Vector::zero().into() }

    /// Check if this is the zero-Point.
    ///
    /// # Examples
    /// ```
    /// use iron_shapes::point::*;
    ///
    /// assert!(Point::<usize>::zero().is_zero());
    /// ```
    pub fn is_zero(&self) -> bool {
        self.x.is_zero() && self.y.is_zero()
    }

    /// Compute the squared distance to the `other` point.
    pub fn distance_sq(self, other: &Point<T>) -> T {
        let diff = self - *other;
        diff.norm2_squared()
    }

    /// Calculate the cross product of the two vectors defined by three points.
    ///
    /// A positive value implies that `self` → `a` → `b` is counter-clockwise, negative implies
    /// clockwise.
    ///
    /// (`b` - `self`) x (`c` - `b`)
    ///
    /// # Examples
    ///
    /// ```
    /// use iron_shapes::point::Point;
    ///
    /// let a = Point::new(1,0);
    /// let b = Point::new(1,1);
    /// let c = Point::new(0,1);
    ///
    /// let p = a.cross_prod3(b, c);
    ///
    /// assert_eq!(p, (b-a).cross_prod(c - b));
    /// ```
    pub fn cross_prod3(&self, b: Point<T>, c: Point<T>) -> T {
        (b.x - self.x) * (c.y - b.y) - (b.y - self.y) * (c.x - b.x)
    }

    /// Return the location of this point as a vector.
    pub fn v(&self) -> Vector<T> {
        self.location
    }
}

impl<T: CoordinateType + NumCast> Point<T> {
    pub fn distance<F: Float>(self, other: &Point<T>) -> F {
        let diff = self - *other;
        diff.length()
    }
}

/// Compare points.
///
/// The ordering is determined by the x-coordinates. If it is the same
/// for both points the y-coordinate is used.
///
/// Point `a` > Point `b` iff `a.x > b.x || (a.x == b.x && a.y > b.y)`.
impl<T: CoordinateType> PartialOrd for Point<T> {
    fn partial_cmp(&self, rhs: &Self) -> Option<Ordering> {
        match self.x.partial_cmp(&rhs.x) {
            Some(Ordering::Equal) => self.y.partial_cmp(&rhs.y),
            maybe_ordering => maybe_ordering
        }
    }
}

/// Compare points.
///
/// The ordering is determined by the x-coordinates. If it is the same
/// for both points the y-coordinate is used.
///
/// Point `a` > Point `b` iff `a.x > b.x || (a.x == b.x && a.y > b.y)`.
impl<T: CoordinateType + Ord> Ord for Point<T> {
    fn cmp(&self, rhs: &Self) -> Ordering {
        match self.x.cmp(&rhs.x) {
            Ordering::Equal => self.y.cmp(&rhs.y),
            ordering => ordering
        }
    }
}


/// Point wise transformation for a single point.
impl<T> MapPointwise<T> for Point<T>
    where T: CoordinateType
{
    /// Point wise transformation.
    fn transform<F>(&self, transformation: F) -> Self
        where F: Fn(Point<T>) -> Point<T> {
        transformation(*self)
    }
}


impl<T: CoordinateType> Into<(T, T)> for Point<T> {
    fn into(self) -> (T, T) {
        (self.x, self.y)
    }
}

impl<T: CoordinateType> Into<(T, T)> for &Point<T> {
    fn into(self) -> (T, T) {
        (self.x, self.y)
    }
}

impl<T: CoordinateType> From<(T, T)> for Point<T> {
    fn from(coords: (T, T)) -> Self {
        Point::new(coords.0, coords.1)
    }
}

impl<'a, T: CoordinateType> From<&'a (T, T)> for Point<T> {
    fn from(coords: &'a (T, T)) -> Self {
        Point::new(coords.0, coords.1)
    }
}

impl<'a, T: CoordinateType> From<&'a Point<T>> for Point<T> {
    fn from(v: &'a Point<T>) -> Self {
        Point::new(v.x, v.y)
    }
}

impl<T: CoordinateType> From<[T; 2]> for Point<T> {
    fn from(coords: [T; 2]) -> Self {
        Point::new(coords[0], coords[1])
    }
}

impl<T> fmt::Debug for Point<T>
    where T: fmt::Debug + CoordinateType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Point({:?}, {:?})", self.x, self.y)
    }
}

impl<T> fmt::Display for Point<T>
    where T: fmt::Display + CoordinateType
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({}, {})", self.x, self.y)
    }
}


impl<T: CoordinateType + NumCast> Point<T> {
    /// Convert Point into a Point with floating point data type.
    pub fn cast_to_float<F: CoordinateType + Float + NumCast>(&self) -> Point<F> {
        // TODO: find conversion that does not panic for sure.
        Point::new(
            F::from(self.x).unwrap(),
            F::from(self.y).unwrap(),
        )
    }
}

impl<T: CoordinateType + NumCast, Dst: CoordinateType + NumCast> TryCastCoord<T, Dst> for Point<T> {
    type Output = Point<Dst>;

    fn try_cast(&self) -> Option<Self::Output> {
        match (Dst::from(self.x), Dst::from(self.y)) {
            (Some(x), Some(y)) => Some(Point::new(x, y)),
            _ => None
        }
    }
}

/// Point addition.
impl<T, V> Add<V> for Point<T>
    where T: CoordinateType + Add<Output=T>,
          V: Into<Point<T>>
{
    type Output = Self;

    fn add(self, rhs: V) -> Self {
        let rhs = rhs.into();
        Point::new(
            self.x + rhs.x,
            self.y + rhs.y,
        )
    }
}

impl<T, V> AddAssign<V> for Point<T>
    where T: CoordinateType + AddAssign<T>,
          V: Into<Vector<T>>
{
    fn add_assign(&mut self, rhs: V) {
        let rhs = rhs.into();
        self.location += rhs;
    }
}

/// Subtract a point.
impl<T> Sub<Point<T>> for Point<T>
    where T: CoordinateType + Sub<Output=T>
{
    type Output = Vector<T>;

    fn sub(self, rhs: Point<T>) -> Self::Output {
        Vector::new(
            self.x - rhs.x,
            self.y - rhs.y,
        )
    }
}

/// Subtract a vector.
impl<T> Sub<Vector<T>> for Point<T>
    where T: CoordinateType + Sub<Output=T>
{
    type Output = Point<T>;

    fn sub(self, rhs: Vector<T>) -> Self::Output {
        Point::new(
            self.x - rhs.x,
            self.y - rhs.y,
        )
    }
}

impl<T, V> SubAssign<V> for Point<T>
    where T: CoordinateType + SubAssign<T>,
          V: Into<Vector<T>>
{
    fn sub_assign(&mut self, rhs: V) {
        let rhs = rhs.into();
        self.location -= rhs;
    }
}

impl<T> Neg for Point<T>
    where T: CoordinateType
{
    type Output = Self;

    fn neg(self) -> Self {
        Point::new(
            T::zero() - self.x,
            T::zero() - self.y,
        )
    }
}

/// Scalar multiplication.
impl<T> Mul<T> for Point<T>
    where T: CoordinateType + Mul<Output=T>
{
    type Output = Self;

    fn mul(self, rhs: T) -> Self {
        Point::new(
            self.x * rhs,
            self.y * rhs,
        )
    }
}

/// In-place scalar multiplication.
impl<T> MulAssign<T> for Point<T>
    where T: CoordinateType + MulAssign<T>
{
    fn mul_assign(&mut self, rhs: T) {
        self.location *= rhs;
    }
}

/// Scalar division.
impl<T> Div<T> for Point<T>
    where T: CoordinateType + Div<Output=T>
{
    type Output = Self;

    fn div(self, rhs: T) -> Self {
        Point::new(
            self.x / rhs,
            self.y / rhs,
        )
    }
}

