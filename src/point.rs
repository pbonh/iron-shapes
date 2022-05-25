// Copyright (c) 2018-2020 Thomas Kramer.
// SPDX-FileCopyrightText: 2018-2022 Thomas Kramer
//
// SPDX-License-Identifier: AGPL-3.0-or-later

//! Points represent a location in the two dimensional plane by an `x` and `y` coordinate.

use crate::vector::Vector;
use crate::traits::{MapPointwise, TryCastCoord, BoundingBox, TryBoundingBox};


use std::cmp::{Ord, Ordering};
use std::fmt;

use num_traits::{Float, NumCast};
pub use num_traits::Zero;
use std::ops::{Div, MulAssign, Mul, Neg, Sub, SubAssign, AddAssign, Add, DerefMut};
pub use std::ops::Deref;
use crate::prelude::Rect;

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
#[derive(Copy, Clone, Default, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Point<T> {
    /// Store the location as a translation from the origin.
    location: Vector<T>
}

impl<T: PartialEq> Eq for Point<T> {}

impl<T: PartialEq> PartialEq for Point<T> {
    fn eq(&self, other: &Self) -> bool {
        self.location == other.location
    }
}

impl<T> Deref for Point<T> {
    type Target = Vector<T>;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.location
    }
}

impl<T> DerefMut for Point<T> {
    #[inline]
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.location
    }
}

impl<T> From<Vector<T>> for Point<T> {
    #[inline]
    fn from(v: Vector<T>) -> Self {
        Point {
            location: v
        }
    }
}

impl<T: Copy> From<&Vector<T>> for Point<T> {
    #[inline]
    fn from(v: &Vector<T>) -> Self {
        Self::from(*v)
    }
}

impl<T> Into<Vector<T>> for Point<T> {
    #[inline]
    fn into(self) -> Vector<T> {
        self.location
    }
}

impl<T> Point<T> {
    /// Create a new point with `x` and `y` coordinates.
    #[inline]
    pub fn new(x: T, y: T) -> Self {
        Point { location: Vector::new(x, y) }
    }
}

impl<T: Copy> Point<T> {
    /// Return the location of this point as a vector.
    #[inline]
    pub fn v(&self) -> Vector<T> {
        self.location
    }
}

impl<T: Zero> Point<T> {
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
    #[inline]
    pub fn zero() -> Self { Vector::zero().into() }

    /// Check if this is the zero-Point.
    ///
    /// # Examples
    /// ```
    /// use iron_shapes::point::*;
    ///
    /// assert!(Point::<usize>::zero().is_zero());
    /// ```
    #[inline]
    pub fn is_zero(&self) -> bool {
        self.x.is_zero() && self.y.is_zero()
    }
}

impl<T: Copy + Mul<Output=T> + Add<Output=T> + Sub<Output=T>> Point<T> {
    /// Compute the squared distance to the `other` point.
    ///
    /// # Examples
    /// ```
    /// use iron_shapes::point::*;
    ///
    /// let a = Point::new(0, 0);
    /// let b = Point::new(2, 0);
    /// assert_eq!(a.distance_sq(&b), 2*2);
    /// ```
    #[inline]
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
    #[inline]
    pub fn cross_prod3(&self, b: Point<T>, c: Point<T>) -> T {
        (b.x - self.x) * (c.y - b.y) - (b.y - self.y) * (c.x - b.x)
    }
}

impl<T: Copy + Sub<Output=T> + NumCast> Point<T> {
    /// Compute the Euclidean distance betwen two points.
    #[inline]
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
impl<T: PartialOrd> PartialOrd for Point<T> {
    #[inline]
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
impl<T: Ord> Ord for Point<T> {
    #[inline]
    fn cmp(&self, rhs: &Self) -> Ordering {
        match self.x.cmp(&rhs.x) {
            Ordering::Equal => self.y.cmp(&rhs.y),
            ordering => ordering
        }
    }
}


/// Point wise transformation for a single point.
impl<T: Copy> MapPointwise<T> for Point<T> {
    /// Point wise transformation.
    #[inline]
    fn transform<F>(&self, transformation: F) -> Self
        where F: Fn(Point<T>) -> Point<T> {
        transformation(*self)
    }
}


impl<T: Copy> Into<(T, T)> for Point<T> {
    #[inline]
    fn into(self) -> (T, T) {
        (self.x, self.y)
    }
}

impl<T: Copy> Into<(T, T)> for &Point<T> {
    #[inline]
    fn into(self) -> (T, T) {
        (self.x, self.y)
    }
}

impl<T: Copy> From<(T, T)> for Point<T> {
    #[inline]
    fn from(coords: (T, T)) -> Self {
        Point::new(coords.0, coords.1)
    }
}

impl<'a, T: Copy> From<&'a (T, T)> for Point<T> {
    #[inline]
    fn from(coords: &'a (T, T)) -> Self {
        Point::new(coords.0, coords.1)
    }
}

impl<'a, T: Copy> From<&'a Point<T>> for Point<T> {
    #[inline]
    fn from(v: &'a Point<T>) -> Self {
        Point::new(v.x, v.y)
    }
}

impl<T: Copy> From<[T; 2]> for Point<T> {
    #[inline]
    fn from(coords: [T; 2]) -> Self {
        Point::new(coords[0], coords[1])
    }
}

impl<T> fmt::Debug for Point<T>
    where T: fmt::Debug {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Point({:?}, {:?})", self.x, self.y)
    }
}

impl<T> fmt::Display for Point<T>
    where T: fmt::Display
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({}, {})", self.x, self.y)
    }
}


impl<T: Copy + NumCast> Point<T> {
    /// Convert Point into a Point with floating point data type.
    #[inline]
    pub fn cast_to_float<F: Float + NumCast>(&self) -> Point<F> {
        // TODO: find conversion that does not panic for sure.
        Point::new(
            F::from(self.x).unwrap(),
            F::from(self.y).unwrap(),
        )
    }
}

impl<T: Copy + NumCast, Dst: Copy + NumCast> TryCastCoord<T, Dst> for Point<T> {
    type Output = Point<Dst>;

    #[inline]
    fn try_cast(&self) -> Option<Self::Output> {
        match (Dst::from(self.x), Dst::from(self.y)) {
            (Some(x), Some(y)) => Some(Point::new(x, y)),
            _ => None
        }
    }
}

/// Point addition.
impl<T, V> Add<V> for Point<T>
    where T: Copy + Add<Output=T>,
          V: Into<Point<T>>
{
    type Output = Self;

    #[inline]
    fn add(self, rhs: V) -> Self {
        let rhs = rhs.into();
        Point::new(
            self.x + rhs.x,
            self.y + rhs.y,
        )
    }
}

impl<T, V> AddAssign<V> for Point<T>
    where T: Copy + AddAssign<T>,
          V: Into<Vector<T>>
{
    #[inline]
    fn add_assign(&mut self, rhs: V) {
        let rhs = rhs.into();
        self.location += rhs;
    }
}

/// Subtract a point.
impl<T> Sub<Point<T>> for Point<T>
    where T: Copy + Sub<Output=T>
{
    type Output = Vector<T>;

    #[inline]
    fn sub(self, rhs: Point<T>) -> Self::Output {
        Vector::new(
            self.x - rhs.x,
            self.y - rhs.y,
        )
    }
}

/// Subtract a vector.
impl<T> Sub<Vector<T>> for Point<T>
    where T: Copy + Sub<Output=T>
{
    type Output = Point<T>;

    #[inline]
    fn sub(self, rhs: Vector<T>) -> Self::Output {
        Point::new(
            self.x - rhs.x,
            self.y - rhs.y,
        )
    }
}

impl<T, V> SubAssign<V> for Point<T>
    where T: Copy + SubAssign<T>,
          V: Into<Vector<T>>
{
    #[inline]
    fn sub_assign(&mut self, rhs: V) {
        let rhs = rhs.into();
        self.location -= rhs;
    }
}

impl<T> Neg for Point<T>
    where T: Copy + Neg<Output=T>
{
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        Point::new(
            -self.x,
            -self.y,
        )
    }
}

/// Scalar multiplication.
impl<T> Mul<T> for Point<T>
    where T: Copy + Mul<Output=T>
{
    type Output = Self;

    #[inline]
    fn mul(self, rhs: T) -> Self {
        Point::new(
            self.x * rhs,
            self.y * rhs,
        )
    }
}

/// In-place scalar multiplication.
impl<T> MulAssign<T> for Point<T>
    where T: Copy + MulAssign<T>
{
    #[inline]
    fn mul_assign(&mut self, rhs: T) {
        self.location *= rhs;
    }
}

/// Scalar division.
impl<T> Div<T> for Point<T>
    where T: Copy + Div<Output=T>
{
    type Output = Self;

    #[inline]
    fn div(self, rhs: T) -> Self {
        Point::new(
            self.x / rhs,
            self.y / rhs,
        )
    }
}

impl<T: Copy + Zero + Add<Output=T>> std::iter::Sum for Point<T> {
    /// Compute the sum of all points in the iterator.
    /// If the iterator is empty, (0, 0) is returned.
    fn sum<I: Iterator<Item=Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, p| acc + p)
    }
}

impl<T: Copy> BoundingBox<T> for Point<T> {
    fn bounding_box(&self) -> Rect<T> {
        Rect {
            lower_left: *self,
            upper_right: *self
        }
    }
}

impl<T: Copy> TryBoundingBox<T> for Point<T> {
    fn try_bounding_box(&self) -> Option<Rect<T>> {
        Some(self.bounding_box())
    }
}

