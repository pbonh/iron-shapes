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

//! This module contains data types and functions for basic polygons without holes.

use crate::CoordinateType;

use crate::point::Point;
use crate::edge::Edge;
use crate::rect::Rect;

pub use crate::traits::{DoubledOrientedArea, TryBoundingBox, MapPointwise, WindingNumber};

use crate::types::*;

use std::iter::FromIterator;
use std::cmp::{Ord, PartialEq};
use std::slice::Iter;
use num_traits::NumCast;
use crate::traits::TryCastCoord;

/// A `SimplePolygon` is a polygon defined by vertices. It does not contain holes but can be
/// self-intersecting.
///
/// TODO: Implement `Deref` for accessing the vertices.
#[derive(Clone, Debug, Hash, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SimplePolygon<T>
    where T: CoordinateType {
    /// Vertices of the polygon.
    pub points: Vec<Point<T>>
}

/// Shorthand notation for creating a simple polygon.
/// # Example
/// ```
/// # #[macro_use]
/// # extern crate iron_shapes;
/// # fn main() {
/// use iron_shapes::prelude::*;
/// let p = simple_polygon!((0, 0), (1, 0), (1, 1));
/// assert_eq!(p, SimplePolygon::from(vec![(0, 0), (1, 0), (1, 1)]));
/// # }
/// ```
#[macro_export]
macro_rules! simple_polygon {
 ($($x:expr),*) => {SimplePolygon::new((vec![$($x.into()),*]))}
}

impl<T: CoordinateType> SimplePolygon<T> {
    /// Create a new polygon from a list of points.
    ///
    /// The orientation of the points is normalized to counter-clock-wise.
    pub fn new(points: Vec<Point<T>>) -> Self {
        let mut new = Self::new_raw(points);

        // Normalize orientation.
        if new.orientation() != Orientation::CounterClockWise {
            new.points.reverse();
        }

        new
    }

    /// Create a new polygon from a list of points.
    /// The points are taken as they are, without reordering
    /// or simplification.
    pub fn new_raw(points: Vec<Point<T>>) -> Self {
        Self {
            points
        }
    }

    /// Create a new simple polygon from a rectangle.
    pub fn from_rect(rect: &Rect<T>) -> Self {
        Self::new(
            vec![rect.lower_left(), rect.lower_right(),
                 rect.upper_right(), rect.upper_left()]
        )
    }

    /// Create empty polygon without any vertices.
    pub fn empty() -> Self {
        SimplePolygon {
            points: Vec::new()
        }
    }

    /// Get the number of vertices.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Shortcut for `self.points.iter()`.
    pub fn iter(&self) -> Iter<Point<T>> {
        self.points.iter()
    }

    /// Get the convex hull of the polygon.
    ///
    /// Implements Andrew's Monotone Chain algorithm.
    /// See: <http://geomalgorithms.com/a10-_hull-1.html>
    pub fn convex_hull(&self) -> SimplePolygon<T>
        where T: Ord {

        // Sort S by increasing x and then y-coordinate.
        let p: Vec<Point<T>> = {
            let mut p: Vec<Point<T>> = self.points.clone();
            p.sort();
            p
        };

        assert!(p.len() >= 2);

        // minmin = index of P with min x first and min y second
        let minmin = 0;
        let p_minmin = p[minmin];
        // maxmax = index of P with max x first and max y second
        let maxmax = p.len() - 1;
        let p_maxmax = p[maxmax];
        // minmax = index of P with min x first and max y second
        let minmax = p.iter().enumerate()
            .take_while(|(_, p)| p.x == p_minmin.x)
            .last().unwrap().0;
        let p_minmax = p[minmax];
        // maxmin = index of P with max x first and min y second
        let maxmin = p.iter().enumerate().rev()
            .take_while(|(_, p)| p.x == p_maxmax.x)
            .last().unwrap().0;
        let p_maxmin = p[maxmin];

        debug_assert!(minmin <= minmax);
        debug_assert!(minmax <= maxmin || p_minmax.x == p_maxmin.x);
        debug_assert!(maxmin <= maxmax);

        debug_assert!(p.iter().all(|&p| p_minmin <= p));
        debug_assert!(p.iter().all(|&p| p_maxmax >= p));
        debug_assert!(p.iter().all(|&p|
            p_minmax.x < p.x || (p_minmax.x == p.x && p_minmax.y >= p.y)
        ));
        debug_assert!(p.iter().all(|&p|
            p_maxmin.x > p.x || (p_maxmin.x == p.x && p_maxmin.y <= p.y)
        ));

        // Handle degenerate case where all x coordinates are equal.
        if p_minmin.x == p_maxmax.x {
            if p_minmin.y == p_maxmax.y {
                SimplePolygon::new(vec![p_minmin])
            } else {
                SimplePolygon::new(vec![p_minmin, p_maxmax])
            }
        } else {
            let build_half_hull =
                |l: Edge<T>, points: &[Point<T>]| {
                    let mut stack = Vec::new();
                    // Push starting point on stack.
                    stack.push(l.start);

                    for &pi in points {
                        // Skip all points that are not strictly right of `l`.
                        if l.side_of(pi) == Side::Right {
                            while stack.len() >= 2 {
                                let pt1 = stack[stack.len() - 1];
                                let pt2 = stack[stack.len() - 2];

                                if Edge::new(pt2, pt1).side_of(pi) == Side::Left {
                                    // `pi` is strictly left of the line defined by the top two elements
                                    // on the stack.
                                    break;
                                }
                                stack.pop();
                            }
                            stack.push(pi);
                        }
                    }

                    stack
                };

            // Compute the lower hull.
            let l_min = Edge::new(p_minmin, p_maxmin);
            let mut lower_half_hull = build_half_hull(l_min, &p[minmax + 1..maxmin]);

            // Push p_maxmin on stack if necessary.
            if p_maxmin != p_maxmax {
                lower_half_hull.push(p_maxmin);
            }

            // Compute the upper hull.
            let l_max = Edge::new(p_maxmax, p_minmax);
            let mut upper_half_hull = build_half_hull(l_max, &p[minmax + 1..maxmin]);

            // Push p_minmax on stack if necessary.
            if p_minmax != p_minmin {
                upper_half_hull.push(p_minmax);
            }

            // Join both hulls.
            lower_half_hull.append(&mut upper_half_hull);

            SimplePolygon::new(lower_half_hull)
        }
    }

    /// Get an iterator over the polygon points.
    /// Point 0 is appended to the end to close the cycle.
    fn iter_cycle(&self) -> impl Iterator<Item=&Point<T>> {
        self.points.iter()
            .cycle()
            .take(self.points.len() + 1)
    }

    /// Get all exterior edges of the polygon.
    /// # Examples
    ///
    /// ```
    /// use iron_shapes::simple_polygon::SimplePolygon;
    /// use iron_shapes::edge::Edge;
    /// let coords = vec![(0, 0), (1, 0)];
    ///
    /// let poly = SimplePolygon::from(coords);
    ///
    /// assert_eq!(poly.edges(), vec![Edge::new((0, 0), (1, 0)), Edge::new((1, 0), (0, 0))]);
    ///
    /// ```
    pub fn edges(&self) -> Vec<Edge<T>> {
        self.edges_iter().collect()
    }

    /// Iterate over all edges.
    fn edges_iter(&self) -> impl Iterator<Item=Edge<T>> + '_ {
        self.iter()
            .zip(self.iter_cycle().skip(1))
            .map(|(a, b)| Edge::new(a, b))
    }

    /// Test if all edges are parallel to the x or y axis.
    pub fn is_rectilinear(&self) -> bool {
        self.edges_iter().all(|e| e.is_rectilinear())
    }

    /// Get the vertex with lowest x-coordinate. Prefer lower y-coordinates to break ties.
    ///
    /// # Examples
    ///
    /// ```
    /// use iron_shapes::simple_polygon::SimplePolygon;
    /// use iron_shapes::point::Point;
    /// let coords = vec![(0, 0), (1, 0), (-1, 2), (-1, 1)];
    ///
    /// let poly = SimplePolygon::from(coords);
    ///
    /// assert_eq!(poly.lower_left_vertex(), Point::new(-1, 1));
    ///
    /// ```
    pub fn lower_left_vertex(&self) -> Point<T> {
        debug_assert!(self.points.len() > 0);

        self.lower_left_vertex_with_index().1
    }

    /// Get the vertex with lowest x-coordinate and its index.
    /// Prefer lower y-coordinates to break ties.
    fn lower_left_vertex_with_index(&self) -> (usize, Point<T>) {
        debug_assert!(self.points.len() > 0);

        // Find minimum.
        let min = self.points
            .iter()
            .enumerate()
            .min_by(|(_, &p1), (_, &p2)|
                p1.partial_cmp(&p2).unwrap());
        let (idx, point) = min.unwrap();

        (idx, point.clone().into())
    }

    /// Get the orientation of the polygon,
    /// i.e. check if it is wound clock-wise or counter-clock-wise.
    ///
    /// # Examples
    ///
    /// ```
    /// use iron_shapes::simple_polygon::SimplePolygon;
    /// use iron_shapes::point::Point;
    /// use iron_shapes::types::Orientation;
    /// let coords = vec![(0, 0), (3, 0), (3, 1)];
    ///
    /// let poly = SimplePolygon::from(coords);
    ///
    /// assert_eq!(poly.orientation(), Orientation::CounterClockWise);
    ///
    /// ```
    pub fn orientation(&self) -> Orientation {
        // Find the orientation by the polygon area.

        let area2 = self.area_doubled_oriented();

        if area2 > T::zero() {
            Orientation::CounterClockWise
        } else if area2 < T::zero() {
            Orientation::ClockWise
        } else {
            debug_assert!(area2 == T::zero());
            Orientation::Straight
        }
    }


    /// Get index of previous vertex.
    fn prev(&self, i: usize) -> usize {
        match i {
            0 => self.points.len() - 1,
            x => x - 1
        }
    }

    /// Get index of next vertex.
    fn next(&self, i: usize) -> usize {
        match i {
            _ if i == self.points.len() - 1 => 0,
            x => x + 1
        }
    }
}

impl<T> WindingNumber<T> for SimplePolygon<T>
    where T: CoordinateType {
    /// Calculate the winding number of the polygon around this point.
    ///
    /// TODO: Define how point on edges and vertices is handled.
    ///
    /// See: <http://geomalgorithms.com/a03-_inclusion.html>
    fn winding_number(&self, point: Point<T>) -> isize {
        let edges = self.edges();
        let mut winding_number = 0isize;

        // Edge Crossing Rules
        //
        // 1. an upward edge includes its starting endpoint, and excludes its final endpoint;
        // 2. a downward edge excludes its starting endpoint, and includes its final endpoint;
        // 3. horizontal edges are excluded
        // 4. the edge-ray intersection point must be strictly right of the point P.

        for e in edges {
            if e.start.y <= point.y { // Crosses upward?
                if e.end.y > point.y { // Crosses really upward?
                    // Yes, crosses upward.
                    if e.side_of(point) == Side::Left {
                        winding_number += 1;
                    }
                }
            } else if e.end.y <= point.y { // Crosses downward?
                // Yes, crosses downward.
                // `e.start.y > point.y` needs not to be checked anymore.
                if e.side_of(point) == Side::Right {
                    winding_number -= 1;
                }
            }
        }

        winding_number
    }
}

/// Create a polygon from a type that is convertible into an iterator of values convertible to `Point`s.
impl<I, T, P> From<I> for SimplePolygon<T>
    where T: CoordinateType,
          I: IntoIterator<Item=P>,
          Point<T>: From<P>
{
    fn from(iter: I) -> Self {
        let points: Vec<Point<T>> = iter.into_iter().map(
            |x| x.into()
        ).collect();

        SimplePolygon { points }
    }
}

// impl<T: CoordinateType> From<&Rect<T>> for SimplePolygon<T> {
//     fn from(rect: &Rect<T>) -> Self {
//         Self::new(
//             vec![rect.lower_left(), rect.lower_right(),
//                  rect.upper_right(), rect.upper_left()]
//         )
//     }
// }

//
// /// Create a polygon from a `Vec` of values convertible to `Point`s.
// impl<'a, T, P> From<&'a Vec<P>> for SimplePolygon<T>
//     where T: CoordinateType,
//           Point<T>: From<&'a P>
// {
//     fn from(vec: &'a Vec<P>) -> Self {
//         let points: Vec<Point<T>> = vec.into_iter().map(
//             |x| x.into()
//         ).collect();
//
//         SimplePolygon { points }
//     }
// }
//
// /// Create a polygon from a `Vec` of values convertible to `Point`s.
// impl<T, P> From<Vec<P>> for SimplePolygon<T>
//     where T: CoordinateType,
//           Point<T>: From<P>
// {
//     fn from(vec: Vec<P>) -> Self {
//         let points: Vec<Point<T>> = vec.into_iter().map(
//             |x| x.into()
//         ).collect();
//
//         SimplePolygon { points }
//     }
// }

/// Create a polygon from a iterator of values convertible to `Point`s.
impl<T, P> FromIterator<P> for SimplePolygon<T>
    where T: CoordinateType,
          P: Into<Point<T>>
{
    fn from_iter<I>(iter: I) -> Self
        where I: IntoIterator<Item=P>
    {
        let points: Vec<Point<T>> = iter.into_iter().map(
            |x| x.into()
        ).collect();

        assert!(points.len() >= 2, "A polygon needs to have at least two points.");

        SimplePolygon { points }
    }
}


impl<T> TryBoundingBox<T> for SimplePolygon<T>
    where T: CoordinateType {
    fn try_bounding_box(&self) -> Option<Rect<T>> {
        if self.len() > 0 {
            let mut x_min = self.points[0].x;
            let mut x_max = x_min;
            let mut y_min = self.points[0].y;
            let mut y_max = y_min;

            for p in self.iter().skip(1) {
                if p.x < x_min {
                    x_min = p.x;
                }
                if p.x > x_max {
                    x_max = p.x;
                }
                if p.y < y_min {
                    y_min = p.y;
                }
                if p.y > y_max {
                    y_max = p.y;
                }
            }

            Some(Rect::new((x_min, y_min), (x_max, y_max)))
        } else {
            None
        }
    }
}

impl<T> MapPointwise<T> for SimplePolygon<T>
    where T: CoordinateType {
    fn transform<F: Fn(Point<T>) -> Point<T>>(&self, tf: F) -> Self {
        let points = self.points.iter().map(|&p| tf(p)).collect();

        let mut new = SimplePolygon {
            points
        };

        // Make sure the polygon is oriented the same way as before.
        // TODO: Could be done more efficiently if the magnification/mirroring of the transformation is known.
        if new.orientation() != self.orientation() {
            new.points.reverse()
        }

        new
    }
}

impl<T: CoordinateType> DoubledOrientedArea<T> for SimplePolygon<T> {
    /// Calculates the doubled oriented area.
    ///
    /// Using doubled area allows to compute in the integers because the area
    /// of a polygon with integer coordinates is either integer or half-integer.
    ///
    /// The area will be positive if the vertices are listed counter-clockwise,
    /// negative otherwise.
    ///
    /// Complexity: O(n)
    ///
    /// # Examples
    ///
    /// ```
    /// use iron_shapes::traits::DoubledOrientedArea;
    /// use iron_shapes::simple_polygon::SimplePolygon;
    /// let coords = vec![(0, 0), (3, 0), (3, 1)];
    ///
    /// let poly = SimplePolygon::from(coords);
    ///
    /// assert_eq!(poly.area_doubled_oriented(), 3);
    ///
    /// ```
    fn area_doubled_oriented(&self) -> T {
        let mut sum = T::zero();
        for i in 0..self.points.len() {
            sum = sum + self.points[i].x * (self.points[self.next(i)].y -
                self.points[self.prev(i)].y);
        }
        sum
    }
}

impl<T> PartialEq for SimplePolygon<T>
    where T: CoordinateType {
    /// Equality test for simple polygons.
    ///
    /// Two polygons are equal iff a cyclic shift on their vertices can be applied
    /// such that the both lists of vertices match exactly.
    ///
    /// Complexity: O(n^2)
    ///
    /// TODO: Normalized ordering of vertices for faster comparison.
    fn eq(&self, rhs: &Self) -> bool {
        let n = self.len();
        if n == rhs.len() {
            for i in 0..n {
                let l = self.points.iter();
                let r = rhs.points.iter().cycle().skip(i).take(n);

                if l.eq(r) {
                    return true;
                }
            }
            false
        } else {
            false
        }
    }
}

impl<T: CoordinateType + NumCast, Dst: CoordinateType + NumCast> TryCastCoord<T, Dst> for SimplePolygon<T> {
    type Output = SimplePolygon<Dst>;

    fn try_cast(&self) -> Option<Self::Output> {
        let new_points: Option<Vec<_>> = self.points.iter()
            .map(|p| p.try_cast())
            .collect();

        new_points.map(|p| SimplePolygon::new_raw(p))
    }
}

/// Two simple polygons should be the same even if points are shifted cyclical.
#[test]
fn test_partial_eq() {
    let p1 = simple_polygon!((0, 0), (0, 1), (1, 1), (1, 0));
    let p2 = simple_polygon!((0, 0), (0, 1), (1, 1), (1, 0));
    assert_eq!(p1, p2);

    let p2 = simple_polygon!((0, 1), (1, 1), (1, 0), (0, 0));
    assert_eq!(p1, p2);
}

/// Simple sanity check for computation of bounding box.
#[test]
fn test_bounding_box() {
    let p = simple_polygon!((0, 0), (0, 1), (1, 1));
    assert_eq!(p.try_bounding_box(), Some(Rect::new((0, 0), (1, 1))));
}
