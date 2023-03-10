// Copyright (c) 2018-2020 Thomas Kramer.
// SPDX-FileCopyrightText: 2018-2022 Thomas Kramer
//
// SPDX-License-Identifier: AGPL-3.0-or-later

//! This module contains data types and functions for basic rectilinear polygons without holes.

use crate::CoordinateType;

use crate::point::Point;
use crate::rect::Rect;

pub use crate::traits::{DoubledOrientedArea, MapPointwise, TryBoundingBox, WindingNumber};

use crate::types::*;

use crate::prelude::SimpleTransform;
use crate::redge::{REdge, REdgeOrientation};
use crate::simple_polygon::SimplePolygon;
use crate::traits::TryCastCoord;
use num_traits::NumCast;
use std::cmp::{Ord, PartialEq};

/// A `SimpleRPolygon` is a rectilinear polygon. It does not contain holes but can be self-intersecting.
/// The vertices are stored in an implicit format (one coordinate of two neighbour vertices is always the same
/// for rectilinear polygons). This reduces memory usage but has the drawback that edges must
/// alternate between horizontal and vertical. Vertices between two edges of the same orientation will
/// be dropped.
///
#[derive(Clone, Debug, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SimpleRPolygon<T> {
    /// Vertices of the polygon.
    /// Begin with a y-coordinate. First edge is horizontal.
    half_points: Vec<T>,
}

impl<T: PartialEq> Eq for SimpleRPolygon<T> {}

/// Shorthand notation for creating a simple polygon.
/// # Example
/// ```
/// # #[macro_use]
/// # extern crate iron_shapes;
/// # fn main() {
/// use iron_shapes::prelude::*;
/// let p = simple_rpolygon!((0, 0), (1, 0), (1, 1), (0, 1));
/// assert_eq!(Some(p), SimpleRPolygon::try_new(vec![(0, 0), (1, 0), (1, 1), (0, 1)]));
/// # }
/// ```
#[macro_export]
macro_rules! simple_rpolygon {
 ($($x:expr),*) => {SimpleRPolygon::try_new((vec![$($x),*])).unwrap()}
}

impl<T> SimpleRPolygon<T> {
    /// Create empty polygon without any vertices.
    pub fn empty() -> Self {
        Self {
            half_points: Vec::new(),
        }
    }

    /// Get the number of vertices.
    #[deprecated(note = "use len() instead")]
    pub fn num_points(&self) -> usize {
        self.half_points.len()
    }

    /// Get the number of vertices.
    pub fn len(&self) -> usize {
        self.half_points.len()
    }

    /// Check if polygon has no vertices.
    pub fn is_empty(&self) -> bool {
        self.half_points.is_empty()
    }

    /// Reverse the order of the vertices in-place.
    pub fn reverse(&mut self) {
        self.half_points.reverse();
        self.half_points.rotate_right(1);
    }

    /// Get index of previous half-point.
    fn prev(&self, i: usize) -> usize {
        match i {
            0 => self.half_points.len() - 1,
            x => x - 1,
        }
    }

    /// Get index of next half-point.
    fn next(&self, i: usize) -> usize {
        match i {
            _ if i == self.half_points.len() - 1 => 0,
            x => x + 1,
        }
    }
}

impl<T: Clone> SimpleRPolygon<T> {
    /// Create a copy of this polygon whose vertices are ordered in reversed order.
    pub fn reversed(&self) -> Self {
        let mut result = self.clone();
        result.reverse();
        result
    }
}

impl<T: Copy> SimpleRPolygon<T> {
    /// Get `i`-th point of the polygon.
    fn get_point(&self, i: usize) -> Point<T> {
        if i % 2 == 0 {
            Point::new(self.half_points[self.prev(i)], self.half_points[i])
        } else {
            Point::new(self.half_points[i], self.half_points[self.prev(i)])
        }
    }

    /// Iterate over the points.
    pub fn points(&self) -> impl Iterator<Item = Point<T>> + '_ {
        (0..self.len()).map(move |i| self.get_point(i))
    }

    /// Get all exterior edges of the polygon.
    /// # Examples
    ///
    /// ```
    /// use iron_shapes::simple_rpolygon::SimpleRPolygon;
    /// use iron_shapes::redge::REdge;
    /// let coords = vec![(0, 0), (1, 0), (1, 1), (0, 1)];
    ///
    /// let poly = SimpleRPolygon::try_new(coords).unwrap();
    /// let edges: Vec<_> = poly.edges().collect();
    /// assert_eq!(edges, vec![
    ///     REdge::new((0, 0), (1, 0)),
    ///     REdge::new((1, 0), (1, 1)),
    ///     REdge::new((1, 1), (0, 1)),
    ///     REdge::new((0, 1), (0, 0)),
    /// ]);
    ///
    /// ```
    pub fn edges(&self) -> impl Iterator<Item = REdge<T>> + '_ {
        (0..self.len()).map(move |i| {
            let orientation = if i % 2 == 0 {
                REdgeOrientation::Horizontal
            } else {
                REdgeOrientation::Vertical
            };

            REdge::new_raw(
                self.half_points[self.prev(i)],
                self.half_points[self.next(i)],
                self.half_points[i],
                orientation,
            )
        })
    }
}

impl<T: CoordinateType> SimpleRPolygon<T> {
    /// Create new rectilinear polygon from points.
    /// Returns `None` if the polygon defined by the points is not rectilinear.
    /// ```
    /// use iron_shapes::simple_rpolygon::SimpleRPolygon;
    ///
    /// let poly1 = SimpleRPolygon::try_new(vec![(0, 0), (1, 0), (1, 1), (0, 1)]);
    /// assert!(poly1.is_some());
    ///
    /// // A triangle cannot be rectilinear.
    /// let poly1 = SimpleRPolygon::try_new(vec![(0, 0), (1, 0), (1, 1)]);
    /// assert!(poly1.is_none());
    ///
    /// ```
    pub fn try_new<P>(points: Vec<P>) -> Option<Self>
    where
        P: Copy + Into<Point<T>>,
    {
        if points.is_empty() {
            // Empty polygon.
            Some(Self {
                half_points: Vec::new(),
            })
        } else if let Some(last) = points.last() {
            let mut half_points = Vec::new();

            let mut last: Point<T> = (*last).into();
            #[derive(PartialEq, Eq, Debug)]
            enum Orientation {
                None,
                Vertical,
                Horizontal,
            }

            let mut orientation = Orientation::None;

            for p in points.iter().cycle().take(points.len() + 1) {
                let p: Point<T> = (*p).into();
                match (last.x == p.x, last.y == p.y) {
                    (true, true) => {
                        // Same point. Do nothing.
                    }
                    (false, false) => {
                        // Not rectilinear.
                        return None;
                    }
                    (false, true) => {
                        // Horizontal line. Store x.
                        if orientation != Orientation::Horizontal {
                            // Corner.
                            if orientation == Orientation::Vertical {
                                half_points.push(last.x);
                            }
                            orientation = Orientation::Horizontal;
                        }
                    }
                    (true, false) => {
                        // Vertical line.
                        if orientation != Orientation::Vertical {
                            // Corner.
                            if orientation == Orientation::Horizontal {
                                half_points.push(last.y);
                            }
                            orientation = Orientation::Vertical;
                        }
                    }
                }
                last = p;
            }

            debug_assert!(half_points.len() % 2 == 0);

            if orientation == Orientation::Vertical {
                // Last edge was horizontal.
                // Make sure to *start* with a horizontal edge.

                if !half_points.is_empty() {
                    half_points.rotate_left(1);
                }
            }
            Some(Self { half_points })
        } else {
            None
        }
    }

    /// Apply the transformation to this rectilinear polygon.
    pub fn transformed(&self, tf: &SimpleTransform<T>) -> Self {
        Self::try_new(self.points().map(|p| tf.transform_point(p)).collect()).unwrap()
        // Unwrap should be fine because the edges will remain axis-aligned under the Simple Transform.
    }

    /// Convert to a `SimplePolygon`.
    pub fn to_simple_polygon(&self) -> SimplePolygon<T> {
        SimplePolygon::new(self.points().collect())
    }

    /// Get the convex hull of the polygon.
    ///
    /// Implements Andrew's Monotone Chain algorithm.
    /// See: <http://geomalgorithms.com/a10-_hull-1.html>
    pub fn convex_hull(&self) -> SimplePolygon<T>
    where
        T: Ord,
    {
        crate::algorithms::convex_hull::convex_hull(self.points().collect())
    }

    /// Get the vertex with lowest x-coordinate. Prefer lower y-coordinates to break ties.
    ///
    /// # Examples
    ///
    /// ```
    /// use iron_shapes::simple_rpolygon::SimpleRPolygon;
    /// use iron_shapes::point::Point;
    /// let coords = vec![(0, 0), (1, 0), (1, 1), (0, 1)];
    ///
    /// let poly = SimpleRPolygon::try_new(coords).unwrap();
    ///
    /// assert_eq!(poly.lower_left_vertex(), Point::new(0, 0));
    ///
    /// ```
    pub fn lower_left_vertex(&self) -> Point<T> {
        debug_assert!(!self.is_empty());

        self.lower_left_vertex_with_index().1
    }

    /// Get the vertex with lowest x-coordinate and its index.
    /// Prefer lower y-coordinates to break ties.
    fn lower_left_vertex_with_index(&self) -> (usize, Point<T>) {
        debug_assert!(!self.is_empty());

        // Find minimum.
        let min = self
            .points()
            .enumerate()
            .min_by(|(_, p1), (_, p2)| p1.partial_cmp(p2).unwrap());
        let (idx, point) = min.unwrap();

        (idx, point)
    }

    /// Get the orientation of the polygon,
    /// i.e. check if it is wound clock-wise or counter-clock-wise.
    ///
    /// # Examples
    ///
    /// ```
    /// use iron_shapes::simple_rpolygon::SimpleRPolygon;
    /// use iron_shapes::point::Point;
    /// use iron_shapes::types::Orientation;
    /// let coords = vec![(0, 0), (1, 0), (1, 1), (0, 1)];
    ///
    /// let poly = SimpleRPolygon::try_new(coords).unwrap();
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

    // TODO:
    // /// Decompose into non-overlapping rectangles.
    // pub fn decompose_rectangles(&self) -> Vec<Rect<T>> {
    //     // Get vertical edges and order them by their lower end.
    //     struct Vertical {
    //         x: T,
    //         y_low: T,
    //         y_high: T,
    //     }
    //     impl Ord for Vertical {
    //         fn cmp(&self, other: &Self) -> Ordering {
    //
    //         }
    //     }
    //     let mut verticals = BinaryHeap::new();
    //     self.edges().filter(|e| e.is_vertical())
    //         .for_each(|e| {
    //             let (y_low, y_high) = if e.start < e.end {
    //                 (e.start, e.end)
    //             } else {
    //                 (e.end, e.start)
    //             };
    //             let v = Vertical {
    //                 x: e.offset,
    //                 y_low,
    //                 y_high,
    //             };
    //             verticals.push(v);
    //         });
    //     unimplemented!()
    // }
}

impl<T> WindingNumber<T> for SimpleRPolygon<T>
where
    T: CoordinateType,
{
    /// Calculate the winding number of the polygon around this point.
    ///
    /// TODO: Define how point on edges and vertices is handled.
    ///
    /// See: <http://geomalgorithms.com/a03-_inclusion.html>
    fn winding_number(&self, point: Point<T>) -> isize {
        let mut winding_number = 0isize;

        // Edge Crossing Rules
        //
        // 1. an upward edge includes its starting endpoint, and excludes its final endpoint;
        // 2. a downward edge excludes its starting endpoint, and includes its final endpoint;
        // 3. horizontal edges are excluded
        // 4. the edge-ray intersection point must be strictly right of the point P.

        for e in self.edges() {
            if e.start().y <= point.y {
                // Crosses upward?
                if e.end().y > point.y {
                    // Crosses really upward?
                    // Yes, crosses upward.
                    if e.side_of(point) == Side::Left {
                        winding_number += 1;
                    }
                }
            } else if e.end().y <= point.y {
                // Crosses downward?
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

impl<T: CoordinateType> From<Rect<T>> for SimpleRPolygon<T> {
    fn from(r: Rect<T>) -> Self {
        Self {
            half_points: vec![
                r.lower_left.y,
                r.upper_right.x,
                r.upper_right.y,
                r.lower_left.x,
            ],
        }
    }
}

#[test]
fn test_from_rect() {
    use super::rect::Rect;
    let r = Rect::new((0, 1), (2, 3));
    let p = SimpleRPolygon::from(r);
    assert_eq!(
        p.points().collect::<Vec<_>>(),
        [
            Point::new(0, 1),
            Point::new(2, 1),
            Point::new(2, 3),
            Point::new(0, 3)
        ]
    );
}

// /// Create a polygon from a type that is convertible into an iterator of values convertible to `Point`s.
// impl<I, T, P> TryFrom<I> for SimpleRPolygon<T>
//     where T: CoordinateType,
//           I: IntoIterator<Item=P>,
//           Point<T>: From<P>
// {
//     type Error = ();
//     /// Create a polygon from a type that is convertible into an iterator of values convertible to `Point`s.
//     /// Return `None` if the polygon is not rectilinear.
//     fn try_from(iter: I) -> Result<Self, ()> {
//         let points: Vec<Point<T>> = iter.into_iter().map(
//             |x| x.into()
//         ).collect();
//
//         match SimpleRPolygon::try_new(points) {
//             None => Err(()),
//             Some(p) => Ok(p)
//         }
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

impl<T> SimpleRPolygon<T>
where
    T: Copy + PartialOrd,
{
    /// Check if the polygon is an axis-aligned rectangle.
    pub fn is_rect(&self) -> bool {
        self.len() == 4
    }
}

impl<T> TryBoundingBox<T> for SimpleRPolygon<T>
where
    T: Copy + PartialOrd,
{
    fn try_bounding_box(&self) -> Option<Rect<T>> {
        if !self.is_empty() {
            let mut xmax = self.half_points[1];
            let mut ymax = self.half_points[0];
            let mut xmin = xmax;
            let mut ymin = ymax;
            self.half_points.chunks(2).for_each(|c| {
                let x = c[1];
                let y = c[0];
                if x > xmax {
                    xmax = x
                };
                if x < xmin {
                    xmin = x
                };
                if y > ymax {
                    ymax = y
                };
                if y < ymin {
                    ymin = y
                };
            });

            Some(Rect::new((xmin, ymin), (xmax, ymax)))
        } else {
            None
        }
    }
}

impl<T: CoordinateType> DoubledOrientedArea<T> for SimpleRPolygon<T> {
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
    /// use iron_shapes::simple_rpolygon::SimpleRPolygon;
    /// let coords = vec![(0, 0), (1, 0), (1, 1), (0, 1)];
    ///
    /// let poly = SimpleRPolygon::try_new(coords).unwrap();
    ///
    /// assert_eq!(poly.area_doubled_oriented(), 2);
    ///
    /// ```
    fn area_doubled_oriented(&self) -> T {
        debug_assert!(self.half_points.len() % 2 == 0);
        // Iterate over all horizontal edges. Compute the area
        // as the sum of (oriented edge length) * (edge distance to origin).
        let area: T = (0..self.len())
            .step_by(2)
            .map(move |i| {
                let start = self.half_points[self.prev(i)];
                let end = self.half_points[self.next(i)];
                let offset = self.half_points[i];

                (start - end) * offset
            })
            .fold(T::zero(), |acc, area| acc + area);

        area + area
    }
}

impl<T: PartialEq> PartialEq for SimpleRPolygon<T> {
    /// Equality test for simple polygons.
    ///
    /// Two polygons are equal iff a cyclic shift on their vertices can be applied
    /// such that the both lists of vertices match exactly.
    ///
    /// Complexity: O(n^2)
    ///
    /// TODO: Normalized ordering of vertices for faster comparison.
    fn eq(&self, rhs: &Self) -> bool {
        let n = self.half_points.len();
        debug_assert!(n % 2 == 0);
        if n == rhs.half_points.len() {
            for i in 0..n / 2 {
                let l = self.half_points.iter();
                let r = rhs.half_points.iter().cycle().skip(2 * i).take(n);

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

impl<T: CoordinateType + NumCast, Dst: CoordinateType + NumCast> TryCastCoord<T, Dst>
    for SimpleRPolygon<T>
{
    type Output = SimpleRPolygon<Dst>;

    fn try_cast(&self) -> Option<Self::Output> {
        let new_half_points: Option<Vec<_>> =
            self.half_points.iter().map(|&p| Dst::from(p)).collect();

        new_half_points.map(|p| SimpleRPolygon { half_points: p })
    }
}

#[test]
fn test_create_rpolygon() {
    let p = SimpleRPolygon::try_new(vec![(0, 0), (1, 0), (1, 2), (0, 2)]).unwrap();
    assert_eq!(p.half_points, vec![0, 1, 2, 0]);

    let p = SimpleRPolygon::try_new(vec![(1, 0), (1, 2), (0, 2), (0, 0)]).unwrap();
    assert_eq!(p.half_points, vec![0, 1, 2, 0]);

    // Zero-area polygon is converted to an empty polygon.
    let p = SimpleRPolygon::try_new(vec![(0, 1), (0, 0)]).unwrap();
    assert_eq!(p.half_points, vec![]);

    // Intermediate vertices on straight lines are removed.
    let p = SimpleRPolygon::try_new(vec![(0, 0), (1, 0), (1, 1), (1, 2), (0, 2), (0, 1)]).unwrap();

    assert_eq!(p.half_points, vec![0, 1, 2, 0]);
}

/// Two simple polygons should be the same even if points are shifted cyclical.
#[test]
fn test_partial_eq() {
    let p1 = simple_rpolygon!((0, 0), (0, 1), (1, 1), (1, 0));
    let p2 = simple_rpolygon!((0, 0), (0, 1), (1, 1), (1, 0));
    assert_eq!(p1, p2);

    let p2 = simple_rpolygon!((0, 1), (1, 1), (1, 0), (0, 0));
    assert_eq!(p1, p2);
}

/// Simple sanity check for computation of bounding box.
#[test]
fn test_bounding_box() {
    let p = simple_rpolygon!((0, 1), (2, 1), (2, 3), (0, 3));
    assert_eq!(p.try_bounding_box(), Some(Rect::new((0, 1), (2, 3))));
}

#[test]
fn test_reversed() {
    let p = simple_rpolygon!((0, 0), (1, 0), (1, 1), (2, 1), (2, 2), (0, 2));
    let p_rev_expected = simple_rpolygon!((0, 2), (2, 2), (2, 1), (1, 1), (1, 0), (0, 0));
    let p_rev_actual = p.reversed();
    assert_eq!(p_rev_actual, p_rev_expected);
    assert_eq!(p.reversed().reversed().half_points, p.half_points)
}
