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

//! An `REdge` is 'rectilinear' edge which is either horizontal or vertical.

use crate::vector::Vector;
use crate::point::Point;
use crate::rect::Rect;
use crate::cmp;

use crate::CoordinateType;

use num_traits::cast::NumCast;

pub use crate::traits::{BoundingBox, RotateOrtho, TryCastCoord};
pub use crate::types::Angle;

pub use crate::types::{Side, ContainsResult};
use crate::edge::Edge;
use std::convert::TryFrom;

/// Return type for the edge-edge intersection functions.
/// Stores all possible results of a edge to edge intersection.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum REdgeIntersection<T: CoordinateType> {
    /// No intersection.
    None,
    /// Intersection in a single point but not on an endpoint of an edge.
    Point(Point<T>),
    /// Intersection in an endpoint of an edge.
    EndPoint(Point<T>),
    /// Full or partial overlap.
    Overlap(REdge<T>),
}

/// Return type for the line-line intersection functions.
/// Stores all possible results of a line to line intersection.
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum RLineIntersection<T: CoordinateType> {
    /// No intersection at all.
    None,
    /// Intersection in a single point.
    /// Besides the intersection point also an other expression for the intersection point is given.
    /// The three values `(a, b, c)` describe the intersection point in terms of a starting point (the starting point
    /// of the edge which defines the line) and the direction of the edge multiplied by a fraction.
    ///
    /// `edge.start + edge.vector()*a/c == p` and
    /// `other_edge.start + other_edge.vector()*b/c == p`.
    Point(Point<T>),
    /// Lines are collinear.
    Collinear,
}

impl<T: CoordinateType> Into<(Point<T>, Point<T>)> for REdge<T> {
    fn into(self) -> (Point<T>, Point<T>) {
        (self.start(), self.end())
    }
}

impl<T: CoordinateType> Into<(Point<T>, Point<T>)> for &REdge<T> {
    fn into(self) -> (Point<T>, Point<T>) {
        (self.start(), self.end())
    }
}

impl<T: CoordinateType> Into<Edge<T>> for &REdge<T> {
    fn into(self) -> Edge<T> {
        Edge::new(self.start(), self.end())
    }
}

impl<T: CoordinateType> TryFrom<&Edge<T>> for REdge<T> {
    type Error = ();

    /// Try to convert an edge into a rectilinear edge.
    /// Returns none if the edge is not rectilinear.
    fn try_from(value: &Edge<T>) -> Result<Self, Self::Error> {
        match REdge::try_from_points(value.start, value.end) {
            None => Err(()),
            Some(e) => Ok(e)
        }
    }
}

// impl<T: CoordinateType> From<(Point<T>, Point<T>)> for REdge<T> {
//     fn from(points: (Point<T>, Point<T>)) -> Self {
//         REdge::new(points.0, points.1)
//     }
// }
//
// impl<T: CoordinateType> From<[Point<T>; 2]> for REdge<T> {
//     fn from(points: [Point<T>; 2]) -> Self {
//         REdge::new(points[0], points[1])
//     }
// }

/// Orientation of a rectilinear edge.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum REdgeOrientation {
    /// Horizontal edge.
    Horizontal,
    /// Vertical edge.
    Vertical,
}

/// An rectilinear edge (horizontal or vertical line segment) is represented by its starting point and end point.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct REdge<T: CoordinateType> {
    /// Start-coordinate of the edge.
    pub start: T,
    /// End-coordinate of the edge.
    pub end: T,
    /// Distance to to the origin (0, 0).
    pub offset: T,
    /// Orientation: Either horizontal or vertical.
    pub orientation: REdgeOrientation,
}

impl<T: CoordinateType> REdge<T> {
    /// Create a new `REdge` from two arguments that implement `Into<Point>`.
    /// The two points must lie either on a vertical or horizontal line, otherwise `None` is returned.
    ///
    /// # Panics
    /// Panics if the two points are not on the same horizontal or vertical line.
    pub fn new<C>(start: C, end: C) -> Self
        where C: Into<Point<T>>
    {
        Self::try_from_points(start, end).expect("Points must be on a vertical or horizontal line.")
    }

    /// Create a new rectilinear edge.
    ///
    /// # Parameters
    ///
    /// * `start`: Start-coordinate of the edge.
    /// * `end`: End-coordinate of the edge.
    /// * `offset`: Distance to to the origin (0, 0).
    /// * `orientation`: Orientation: Either horizontal or vertical.
    pub fn new_raw(start: T, end: T, offset: T, orientation: REdgeOrientation) -> Self {
        Self {
            start,
            end,
            offset,
            orientation,
        }
    }

    /// Create a new `REdge` from two arguments that implement `Into<Point>`.
    /// The two points must lie either on a vertical or horizontal line, otherwise `None` is returned.
    pub fn try_from_points<C>(start: C, end: C) -> Option<Self>
        where C: Into<Point<T>> {
        let s = start.into();
        let e = end.into();

        if s.x == e.x {
            // Vertical edge.
            Some(REdge {
                start: s.y,
                end: e.y,
                offset: s.x,
                orientation: REdgeOrientation::Vertical,
            })
        } else if s.y == e.y {
            // Horizontal edge.
            Some(REdge {
                start: s.x,
                end: e.x,
                offset: s.y,
                orientation: REdgeOrientation::Horizontal,
            })
        } else {
            // Edge is neither horizontal nor vertical.
            None
        }
    }

    /// Get the start point of the edge.
    pub fn start(&self) -> Point<T> {
        match self.orientation {
            REdgeOrientation::Horizontal => Point::new(self.start, self.offset),
            REdgeOrientation::Vertical => Point::new(self.offset, self.start)
        }
    }

    /// Get the end point of the edge.
    pub fn end(&self) -> Point<T> {
        match self.orientation {
            REdgeOrientation::Horizontal => Point::new(self.end, self.offset),
            REdgeOrientation::Vertical => Point::new(self.offset, self.end)
        }
    }

    /// Return the same edge but with the two points swapped.
    pub fn reversed(&self) -> Self {
        Self {
            start: self.end,
            end: self.start,
            offset: self.offset,
            orientation: self.orientation,
        }
    }

    /// Check if edge is degenerate.
    /// An edge is degenerate if start point and end point are equal.
    #[inline]
    pub fn is_degenerate(&self) -> bool {
        self.start == self.end
    }

    /// Test if this edge is either horizontal or vertical.
    #[inline]
    pub fn is_ortho(&self) -> bool {
        !self.is_degenerate() && true
    }

    /// Test if this edge is horizontal.
    #[inline]
    pub fn is_horizontal(&self) -> bool {
        !self.is_degenerate() &&
            self.orientation == REdgeOrientation::Horizontal
    }

    /// Test if this edge is vertical.
    #[inline]
    pub fn is_vertical(&self) -> bool {
        !self.is_degenerate() &&
            self.orientation == REdgeOrientation::Vertical
    }


    /// Returns the vector from `self.start()` to `self.end()`.
    pub fn vector(&self) -> Vector<T> {
        match self.orientation {
            REdgeOrientation::Horizontal => Vector::new(self.end - self.start, T::zero()),
            REdgeOrientation::Vertical => Vector::new(T::zero(), self.end - self.start)
        }
    }

    /// Tells on which side of the edge a point is.
    ///
    /// # Panics
    /// Panics if the edge is degenerate.
    ///
    /// Returns `Side::Left` if the point is on the left side,
    /// `Side::Right` if the point is on the right side
    /// or `Side::Center` if the point lies exactly on the line.
    pub fn side_of(&self, point: Point<T>) -> Side {
        assert!(!self.is_degenerate(), "Edge is degenerate.");

        let p_offset = match self.orientation {
            REdgeOrientation::Horizontal => T::zero() - point.y,
            REdgeOrientation::Vertical => point.x
        };

        if p_offset < self.offset {
            if self.start < self.end {
                Side::Left
            } else {
                Side::Right
            }
        } else if p_offset > self.offset {
            if self.start < self.end {
                Side::Right
            } else {
                Side::Left
            }
        } else {
            debug_assert!(p_offset == self.offset);
            Side::Center
        }
    }


    //    /// Find minimal distance between two edges.
//    pub fn distance_to_edge(self, other: Edge<T>) -> f64 {
//        let d1 = self.distance(other.start);
//        let d2 = self.distance(other.end);
//        let d3 = other.distance(self.start);
//        let d4 = other.distance(self.end);
//
//        let min1 = d1.min(d2);
//        let min2 = d3.min(d4);
//
//        // TODO: if intersects
//
//        min1.min(min2)
//    }


    /// Test if point lies on the edge.
    /// Includes start and end points of edge.
    pub fn contains_point(&self, point: Point<T>) -> ContainsResult {
        let (p_offset, p_projected) = match self.orientation {
            REdgeOrientation::Horizontal => (point.y, point.x),
            REdgeOrientation::Vertical => (point.x, point.y)
        };

        if p_offset != self.offset || self.is_degenerate() {
            ContainsResult::No
        } else if p_projected == self.start || p_projected == self.end {
            ContainsResult::OnBounds
        } else if (p_projected >= self.start && p_projected <= self.end)
            || (p_projected >= self.end && p_projected <= self.start) {
            ContainsResult::WithinBounds
        } else {
            ContainsResult::No
        }
    }

    /// Test if point lies on the line defined by the edge.
    pub fn line_contains_point(&self, point: Point<T>) -> bool {
        let p_offset = match self.orientation {
            REdgeOrientation::Horizontal => point.y,
            REdgeOrientation::Vertical => point.x
        };
        p_offset == self.offset
    }


    /// Test if two edges are parallel.
    pub fn is_parallel(&self, other: &REdge<T>) -> bool {
        self.orientation == other.orientation
    }


    /// Test if two edges are collinear, i.e. are on the same line.
    pub fn is_collinear(&self, other: &REdge<T>) -> bool
        where T: CoordinateType {
        self.is_parallel(other) && self.offset == other.offset
    }

    /// Test edges for coincidence.
    /// Two edges are coincident if they are oriented the same way
    /// and share more than one point (implies that they must be parallel).
    pub fn is_coincident(&self, other: &REdge<T>) -> bool {
        // self.is_collinear(other) &&
        //     !Interval::new(self.start, self.end)
        //         .intersection(&zInterval::new(other.start, other.end))
        //         .is_empty()

        self.is_collinear(other) && (
            self.start <= other.start && other.start <= self.end ||
                self.start <= other.end && other.end <= self.end ||
                self.end <= other.start && other.start <= self.start ||
                self.end <= other.end && other.end <= self.start ||
                other.start <= self.start && self.start <= other.end ||
                other.start <= self.end && self.end <= other.end ||
                other.end <= self.start && self.start <= other.start ||
                other.end <= self.end && self.end <= other.start
        )
    }


    /// Test if this edge is crossed by the line defined by the other edge.
    ///
    /// Returns `WithinBounds` if start and end point of this edge lie on different sides
    /// of the line defined by the `other` edge or `OnBounds` if at least one of the points
    /// lies on the line.
    pub fn crossed_by_line(&self, line: &REdge<T>) -> ContainsResult {
        // TODO: Handle degenerate cases.
        let side1 = line.side_of(self.start());

        if side1 == Side::Center {
            ContainsResult::OnBounds
        } else {
            let side2 = line.side_of(self.end());

            if side2 == Side::Center {
                ContainsResult::OnBounds
            } else {
                if side1 == side2 {
                    ContainsResult::No
                } else {
                    ContainsResult::WithinBounds
                }
            }
        }
    }

    /// Test if lines defined by the edges intersect.
    /// If the lines are collinear they are also considered intersecting.
    pub fn lines_intersect(&self, other: &REdge<T>) -> bool {
        !self.is_parallel(other) || self.is_collinear(other)
    }

    /// Test if two edges intersect.
    /// If the edges coincide, they also intersect.
    pub fn edges_intersect(&self, other: &REdge<T>) -> ContainsResult {

        // Two edges intersect if the start and end point of one edge
        // lie on opposite sides of the other edge.

        match self.edge_intersection(other) {
            REdgeIntersection::None => ContainsResult::No,
            REdgeIntersection::Overlap(_) => ContainsResult::WithinBounds,
            REdgeIntersection::Point(_) => ContainsResult::WithinBounds,
            REdgeIntersection::EndPoint(_) => ContainsResult::OnBounds,
        }
    }


    /// Calculate the distance from the point to the line given by the edge.
    ///
    /// Distance will be positive if the point lies on the right side of the edge and negative
    /// if the point is on the left side.
    pub fn oriented_distance_to_line(&self, point: Point<T>) -> T {
        assert!(!self.is_degenerate());

        let diff = match self.orientation {
            REdgeOrientation::Horizontal => self.offset - point.y,
            REdgeOrientation::Vertical => point.x - self.offset
        };

        if self.start > self.end {
            T::zero() - diff
        } else {
            diff
        }
    }

    /// Calculate the distance from the point to the line given by the edge.
    pub fn distance_to_line(&self, point: Point<T>) -> T {
        let diff = self.oriented_distance_to_line(point);
        if diff < T::zero() {
            T::zero() - diff
        } else {
            diff
        }
    }


    /// Find the perpendicular projection of a point onto the line of the edge.
    pub fn projection(&self, point: Point<T>) -> Point<T> {
        assert!(!self.is_degenerate());

        match self.orientation {
            REdgeOrientation::Horizontal => (point.x, self.offset),
            REdgeOrientation::Vertical => (self.offset, point.y)
        }.into()
    }

    /// Compute the intersection of the two lines defined by the edges.
    pub fn line_intersection(&self, other: &REdge<T>) -> RLineIntersection<T> {
        match (self.orientation, other.orientation) {
            (REdgeOrientation::Horizontal, REdgeOrientation::Horizontal) |
            (REdgeOrientation::Vertical, REdgeOrientation::Vertical) => {
                if self.offset == other.offset {
                    RLineIntersection::Collinear
                } else {
                    RLineIntersection::None
                }
            }
            (REdgeOrientation::Horizontal, REdgeOrientation::Vertical) => {
                RLineIntersection::Point(Point::new(other.offset, self.offset))
            }
            (REdgeOrientation::Vertical, REdgeOrientation::Horizontal) => {
                RLineIntersection::Point(Point::new(self.offset, other.offset))
            }
        }
    }

    /// Compute the intersection between two edges.
    pub fn edge_intersection(&self, other: &REdge<T>) -> REdgeIntersection<T> {
        match (self.orientation, other.orientation) {
            (REdgeOrientation::Horizontal, REdgeOrientation::Horizontal) |
            (REdgeOrientation::Vertical, REdgeOrientation::Vertical) => {
                if self.offset == other.offset {
                    // Sort start and end such that start comes first.
                    let (s1, e1) = if self.start < self.end {
                        (self.start, self.end)
                    } else {
                        (self.end, self.start)
                    };
                    let (s2, e2) = if other.start < other.end {
                        (other.start, other.end)
                    } else {
                        (other.end, other.start)
                    };
                    debug_assert!(s1 <= e1);
                    debug_assert!(s2 <= e2);

                    // Compute the intersection of the two intervals.
                    let s = cmp::max(s1, s2);
                    let e = cmp::min(e1, e2);

                    if s > e {
                        REdgeIntersection::None
                    } else if s < e {
                        // Make sure the orientation is the same as for `self`.
                        let (s, e) = if self.start < self.end {
                            (s, e)
                        } else {
                            (e, s)
                        };
                        REdgeIntersection::Overlap(REdge {
                            start: s,
                            end: e,
                            offset: self.offset,
                            orientation: self.orientation,
                        })
                    } else {
                        debug_assert!(s == e);
                        // Intersection in an endpoint.
                        let p = match self.orientation {
                            REdgeOrientation::Vertical => Point::new(self.offset, s),
                            REdgeOrientation::Horizontal => Point::new(s, self.offset),
                        };
                        REdgeIntersection::EndPoint(p)
                    }
                } else {
                    REdgeIntersection::None
                }
            }
            (o1, o2) => {
                let (horizontal, vertical) = match (o1, o2) {
                    (REdgeOrientation::Horizontal, REdgeOrientation::Vertical) => (self, other),
                    (REdgeOrientation::Vertical, REdgeOrientation::Horizontal) => (other, self),
                    _ => panic!()
                };
                let p = Point::new(vertical.offset, horizontal.offset);
                let is_on_horizontal = (horizontal.start <= p.x && p.x <= horizontal.end) ||
                    (horizontal.end <= p.x && p.x <= horizontal.start);
                let is_on_vertical = (vertical.start <= p.y && p.y <= vertical.end) ||
                    (vertical.end <= p.y && p.y <= vertical.start);
                if is_on_horizontal && is_on_vertical {
                    let is_endpoint_horizontal = p.x == horizontal.start || p.x == horizontal.end;
                    let is_endpoint_vertical = p.y == vertical.start || p.x == vertical.end;
                    if is_endpoint_horizontal || is_endpoint_vertical {
                        REdgeIntersection::EndPoint(p)
                    } else {
                        REdgeIntersection::Point(p)
                    }
                } else {
                    REdgeIntersection::None
                }
            }
        }
    }
}


impl<T: CoordinateType> BoundingBox<T> for REdge<T> {
    fn bounding_box(&self) -> Rect<T> {
        Rect::new(self.start(), self.end())
    }
}

impl<T: CoordinateType + NumCast, Dst: CoordinateType + NumCast> TryCastCoord<T, Dst> for REdge<T> {
    type Output = REdge<Dst>;

    fn try_cast(&self) -> Option<Self::Output> {
        match (Dst::from(self.start), Dst::from(self.end), Dst::from(self.offset)) {
            (Some(s), Some(e), Some(o)) => Some(
                REdge {
                    start: s,
                    end: e,
                    offset: o,
                    orientation: self.orientation,
                }
            ),
            _ => None
        }
    }
}


#[cfg(test)]
mod tests {
    use crate::redge::{REdge, RLineIntersection, REdgeIntersection};
    use crate::point::Point;
    use crate::types::*;

    #[test]
    fn test_is_parallel() {
        let e1 = REdge::new((0, 0), (1, 0));
        let e2 = REdge::new((100, 200), (101, 200));
        let e3 = REdge::new((1000, 2000), (1000, 2001));

        assert!(e1.is_parallel(&e2));
        assert!(!e1.is_parallel(&e3));
    }

    #[test]
    fn test_is_collinear() {
        let e1 = REdge::new((0, 0), (1, 0));
        let e2 = REdge::new((3, 0), (4, 0));
        assert!(e1.is_collinear(&e2));
        assert!(e2.is_collinear(&e1));

        // Not collinear.
        let e1 = REdge::new((0, 0), (1, 0));
        let e2 = REdge::new((3, 1), (4, 1));
        assert!(!e1.is_collinear(&e2));
        assert!(!e2.is_collinear(&e1));

        // Not collinear.
        let e1 = REdge::new((0, 0), (1, 0));
        let e2 = REdge::new((0, 0), (0, 1));
        assert!(!e1.is_collinear(&e2));
        assert!(!e2.is_collinear(&e1));
    }

    #[test]
    fn test_oriented_distance_to_line() {
        let xaxis = REdge::new((0, 0), (1, 0));
        let xaxis_neg = REdge::new((0, 0), (-1, 0));
        let yaxis = REdge::new((0, 0), (0, 1));
        let yaxis_neg = REdge::new((0, 0), (0, -1));

        let p = Point::new(2, 3);

        assert_eq!(xaxis.oriented_distance_to_line(p), -p.y);
        assert_eq!(xaxis_neg.oriented_distance_to_line(p), p.y);
        assert_eq!(yaxis.oriented_distance_to_line(p), p.x);
        assert_eq!(yaxis_neg.oriented_distance_to_line(p), -p.x);
    }

    #[test]
    fn test_distance_to_line() {
        let xaxis = REdge::new((0, 0), (1, 0));
        let yaxis = REdge::new((0, 0), (0, 1));
        let p = Point::new(2, 3);

        assert_eq!(xaxis.distance_to_line(p), p.y);
        assert_eq!(yaxis.distance_to_line(p), p.x);
    }

    #[test]
    fn test_line_contains_point() {
        let xaxis = REdge::new((0, 0), (1, 0));
        let yaxis = REdge::new((0, 0), (0, 1));
        let p0 = Point::new(2, 0);
        let p1 = Point::new(2, 3);
        let p2 = Point::new(0, 3);

        assert!(xaxis.line_contains_point(p0));
        assert!(!yaxis.line_contains_point(p0));
        assert!(!xaxis.line_contains_point(p1));
        assert!(yaxis.line_contains_point(p2));
        assert!(!xaxis.line_contains_point(p2));
    }

    #[test]
    fn test_contains() {
        let e1 = REdge::new((0, 0), (10, 0));
        let p0 = Point::new(0, 0);
        let p1 = Point::new(1, 0);
        let p2 = Point::new(0, 1);
        let p3 = Point::new(10, 0);
        let p4 = Point::new(11, 0);

        assert!(e1.contains_point(p0).inclusive_bounds());
        assert!(!e1.contains_point(p0).is_within_bounds());
        assert!(e1.contains_point(p1).inclusive_bounds());
        assert!(e1.contains_point(p2).is_no());
        assert!(e1.contains_point(p3).inclusive_bounds());
        assert!(e1.contains_point(p4).is_no());
    }

    #[test]
    fn test_projection() {
        let horizontal = REdge::new((0, 0), (1, 0));
        let vertical = REdge::new((0, 0), (0, 1));
        let p = Point::new(1, 2);

        assert_eq!(horizontal.projection(p), Point::new(p.x, 0));
        assert_eq!(vertical.projection(p), Point::new(0, p.y));
    }

    #[test]
    fn test_side_of() {
        let xaxis = REdge::new((0, 0), (1, 0));
        let yaxis = REdge::new((0, 0), (0, 1));

        let p1 = Point::new(-10, 0);
        let p2 = Point::new(10, -10);
        let p3 = Point::new(0, 10);

        assert_eq!(yaxis.side_of(p1), Side::Left);
        assert_eq!(yaxis.side_of(p2), Side::Right);
        assert_eq!(yaxis.side_of(p3), Side::Center);

        assert_eq!(xaxis.side_of(p1), Side::Center);
        assert_eq!(xaxis.side_of(p2), Side::Right);
        assert_eq!(xaxis.side_of(p3), Side::Left);
    }

    #[test]
    fn test_crossed_by() {
        let e1 = REdge::new((1, 0), (3, 0));
        let e2 = REdge::new((1, 1), (3, 1));

        let e3 = REdge::new((2, -1), (2, 1));
        let e4 = REdge::new((2, 100), (2, 101));

        // Coincident lines
        assert!(e1.crossed_by_line(&e1).inclusive_bounds());

        // Parallel but not coincident.
        assert!(e1.is_parallel(&e2));
        assert!(!e1.crossed_by_line(&e2).inclusive_bounds());

        // Crossing lines.
        assert!(e1.crossed_by_line(&e3).inclusive_bounds());
        assert!(e3.crossed_by_line(&e1).inclusive_bounds());

        // crossed_by is not commutative
        assert!(e1.crossed_by_line(&e4).inclusive_bounds());
        assert!(!e4.crossed_by_line(&e1).inclusive_bounds());
    }

    #[test]
    fn test_intersect() {
        let e1 = REdge::new((0, 0), (2, 0));
        let e2 = REdge::new((1, -1), (1, 1));
        let e3 = REdge::new((1, 100), (1, 101));

        // Self intersection.
        assert!(e1.edges_intersect(&e1).inclusive_bounds());

        assert!(e1.edges_intersect(&e2).inclusive_bounds());
        assert!(e2.edges_intersect(&e1).inclusive_bounds());

        assert_eq!(e3.side_of(e1.start()), Side::Left);
        assert_eq!(e3.side_of(e1.end()), Side::Right);
        assert!(e1.crossed_by_line(&e3).inclusive_bounds());
        assert!(!e1.edges_intersect(&e3).inclusive_bounds());

        // Intersection at an endpoint.
        let e1 = REdge::new((0, 0), (0, 1));
        let e2 = REdge::new((0, 1), (2, 1));
        assert_eq!(e1.edges_intersect(&e2), ContainsResult::OnBounds);
        assert_eq!(e2.edges_intersect(&e1), ContainsResult::OnBounds);
    }


    #[test]
    fn test_line_intersection() {
        let v = REdge::new((7, 1), (7, 2));
        let h = REdge::new((100, 77), (101, 77));

        assert_eq!(v.line_intersection(&h), RLineIntersection::Point(Point::new(7, 77)));

        assert_eq!(v.line_intersection(&v), RLineIntersection::Collinear);
        assert_eq!(h.line_intersection(&h), RLineIntersection::Collinear);

        let v1 = REdge::new((0, 1), (0, 2));
        let v2 = REdge::new((1, 1), (1, 2));
        assert_eq!(v1.line_intersection(&v2), RLineIntersection::None);
    }


    #[test]
    fn test_edge_intersection() {
        // Point intersection inside both edges.
        let e1 = REdge::new((0, 0), (0, 2));
        let e2 = REdge::new((-1, 1), (1, 1));
        assert_eq!(e1.edge_intersection(&e2), REdgeIntersection::Point(Point::new(0, 1)));
        assert_eq!(e2.edge_intersection(&e1), REdgeIntersection::Point(Point::new(0, 1)));

        // Point intersection on the end of one edge.
        let e1 = REdge::new((0, 0), (0, 1));
        let e2 = REdge::new((-1, 0), (1, 0));
        assert_eq!(e1.edge_intersection(&e2), REdgeIntersection::EndPoint(Point::new(0, 0)));
        assert_eq!(e2.edge_intersection(&e1), REdgeIntersection::EndPoint(Point::new(0, 0)));

        // No intersection
        let e1 = REdge::new((0, 0), (0, 1));
        let e2 = REdge::new((0, 2), (0, 3));
        assert_eq!(e1.edge_intersection(&e2), REdgeIntersection::None);
        assert_eq!(e2.edge_intersection(&e1), REdgeIntersection::None);

        let e1 = REdge::new((0, 0), (0, 1));
        let e2 = REdge::new((1, 0), (2, 0));
        assert_eq!(e1.edge_intersection(&e2), REdgeIntersection::None);
        assert_eq!(e2.edge_intersection(&e1), REdgeIntersection::None);
    }

    #[test]
    fn test_edge_intersection_overlap() {

        // Overlapping edges. Same orientations.
        let e1 = REdge::new((1, 1), (1, 3));
        let e2 = REdge::new((1, 1), (1, 2));
        assert!(e1.edges_intersect(&e2).inclusive_bounds());
        assert!(e1.is_collinear(&e2));

        assert_eq!(e1.edge_intersection(&e2),
                   REdgeIntersection::Overlap(e2));
        assert_eq!(e2.edge_intersection(&e1),
                   REdgeIntersection::Overlap(e2));


        // Overlapping edges. Opposing orientations.
        let e1 = REdge::new((1, 1), (1, 3));
        let e2 = REdge::new((1, 2), (1, 1));
        assert!(e1.edges_intersect(&e2).inclusive_bounds());
        assert!(e1.is_collinear(&e2));

        assert_eq!(e1.edge_intersection(&e2),
                   REdgeIntersection::Overlap(e2.reversed()));
        assert_eq!(e2.edge_intersection(&e1),
                   REdgeIntersection::Overlap(e2));


        // Overlapping edges. One is fully contained in the other. Same orientations.
        let e1 = REdge::new((1, 1), (1, 100));
        let e2 = REdge::new((1, 2), (1, 3));
        assert!(e1.edges_intersect(&e2).inclusive_bounds());
        assert!(e1.is_collinear(&e2));

        assert_eq!(e1.edge_intersection(&e2),
                   REdgeIntersection::Overlap(e2));
        assert_eq!(e2.edge_intersection(&e1),
                   REdgeIntersection::Overlap(e2));

        // Overlapping edges. One is fully contained in the other. Opposing orientations.
        let e1 = REdge::new((1, 1), (1, 100));
        let e2 = REdge::new((1, 3), (1, 2));
        assert!(e1.edges_intersect(&e2).inclusive_bounds());
        assert!(e1.is_collinear(&e2));

        assert_eq!(e1.edge_intersection(&e2),
                   REdgeIntersection::Overlap(e2.reversed()));
        assert_eq!(e2.edge_intersection(&e1),
                   REdgeIntersection::Overlap(e2));


        // Collinear edge, touch in exactly one point.
        let e1 = REdge::new((0, 0), (0, 1));
        let e2 = REdge::new((0, 1), (0, 2));
        assert!(e1.edges_intersect(&e2).inclusive_bounds());
        assert!(e1.is_collinear(&e2));

        assert_eq!(e1.edge_intersection(&e2),
                   REdgeIntersection::EndPoint(Point::new(0, 1)));
        assert_eq!(e2.edge_intersection(&e1),
                   REdgeIntersection::EndPoint(Point::new(0, 1)));
    }

    #[test]
    fn test_coincident() {
        let e1 = REdge::new((0, 0), (0, 2));
        let e2 = REdge::new((0, 1), (0, 3));
        assert!(e1.is_coincident(&e2));
        assert!(e2.is_coincident(&e1));

        let e1 = REdge::new((0, 0), (0, 1));
        let e2 = REdge::new((0, 1), (0, 3));
        assert!(e1.is_coincident(&e2));
        assert!(e2.is_coincident(&e1));

        let e1 = REdge::new((0, 0), (2, 0));
        let e2 = REdge::new((1, 0), (3, 0));
        assert!(e1.is_coincident(&e2));
        assert!(e2.is_coincident(&e1));


        let e1 = REdge::new((0, 0), (0, 1));
        let e2 = REdge::new((0, 0), (1, 0));
        assert!(!e1.is_coincident(&e2));
        assert!(!e2.is_coincident(&e1));
    }


    // #[test]
    // fn test_intersection_of_random_edges() {
    //     let tol = 1e-12;
    //     let seed1 = [1u8; 32];
    //     let seed2 = [2u8; 32];
    //
    //     let between = Uniform::from(-1.0..1.0);
    //     let mut rng = StdRng::from_sed(seed1);
    //
    //     let mut rand_edge = || -> Edge<f64> {
    //         let points: Vec<(f64, f64)> = (0..2).into_iter()
    //             .map(|_| (between.sample(&mut rng), between.sample(&mut rng)))
    //             .collect();
    //
    //         REdge::new(points[0], points[1])
    //     };
    //
    //
    //     let bernoulli_rare = Bernoulli::new(0.05).unwrap();
    //     let bernoulli_50 = Bernoulli::new(0.5).unwrap();
    //     let mut rng = StdRng::from_seed(seed2);
    //
    //     for _i in 0..1000 {
    //
    //         // Create a random pair of edges. 5% of the edge pairs share an endpoint.
    //         let (a, b) = {
    //             let a = rand_edge();
    //
    //             let b = {
    //                 let b = rand_edge();
    //
    //                 if bernoulli_rare.sample(&mut rng) {
    //                     // Share a point with the other edge.
    //                     let shared = if bernoulli_50.sample(&mut rng) {
    //                         a.start
    //                     } else {
    //                         a.end
    //                     };
    //
    //                     let result = REdge::new(b.start, shared);
    //                     if bernoulli_50.sample(&mut rng) {
    //                         result
    //                     } else {
    //                         result.reversed()
    //                     }
    //                 } else {
    //                     b
    //                 }
    //             };
    //             (a, b)
    //         };
    //
    //         let intersection_ab = a.edge_intersection_approx(&b, tol);
    //         assert_eq!(intersection_ab != EdgeIntersection::None, a.edges_intersect(&b).inclusive_bounds());
    //
    //         let intersection_ba = b.edge_intersection_approx(&a, tol);
    //         assert_eq!(intersection_ba != EdgeIntersection::None, b.edges_intersect(&a).inclusive_bounds());
    //     }
    // }
}