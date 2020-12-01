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
use crate::traits::Transform;
use crate::CoordinateType;

use std::cmp::{Ord, Ordering};

use num_traits::{Float, NumCast};

/// A point is defined by a x and y coordinate in the euclidean plane.
pub type Point<T> = Vector<T>;

/// Shorthand notation for creating a point.
macro_rules! point {
 ($x:expr, $y:expr) => {Point::new($x, $y)}
}

//#[derive(Copy, Clone, Hash, PartialEq, Eq, Debug)]
//pub struct Point<T>
//    where T: CoordinateType {
//    pub x: T,
//    pub y: T,
//}

//impl<T: CoordinateType> From<Vector<T>> for Point<T> {
//    fn from(v: Vector<T>) -> Self {
//        Point {
//            x: v.x,
//            y: v.y,
//        }
//    }
//}

//impl<T: CoordinateType> Into<Vector<T>> for Point<T> {
//    fn into(self) -> Vector<T> {
//        Vector {
//            x: self.x,
//            y: self.y,
//        }
//    }
//}

///// Convert a type into a point by converting it into a vector first.
//impl<S, T: CoordinateType> From<S> for Point<T>
//    where S: Into<Vector<T>> {
//    fn from(s: S) -> Self {
//        let v: Vector<T> = s.into();
//        v.into()
//    }
//}

impl<T: CoordinateType> Point<T> {
    pub fn distance_sq(self, other: &Point<T>) -> T {
        let diff: Vector<T> = self - *other;
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
impl<T> Transform<T> for Vector<T>
    where T: CoordinateType
{
    /// Point wise transformation.
    fn transform<F>(&self, transformation: F) -> Self
        where F: Fn(Point<T>) -> Point<T> {
        transformation(*self)
    }
}