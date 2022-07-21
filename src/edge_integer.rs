// Copyright (c) 2018-2020 Thomas Kramer.
// SPDX-FileCopyrightText: 2018-2022 Thomas Kramer
//
// SPDX-License-Identifier: AGPL-3.0-or-later

//! Edge intersection functions for integer coordinates.

use std::cmp::Ordering;
use crate::point::Point;
pub use crate::edge::{Edge, EdgeIntersection, LineIntersection};
use crate::redge::{REdge, REdgeIntersection};

use crate::CoordinateType;

use num_traits::{PrimInt, Zero};
use crate::traits::BoundingBox;
use std::fmt::Debug;
use std::convert::TryFrom;

impl<T: CoordinateType + PrimInt + Debug> Edge<T> {
    /// Compute the intersection point of the lines defined by the two edges.
    /// Coordinates of intersection points are rounded towards zero.
    ///
    /// Degenerate lines don't intersect by definition.
    ///
    /// Returns `LineIntersection::None` iff the two lines don't intersect.
    /// Returns `LineIntersection::Collinear` iff both lines are equal.
    /// Returns `LineIntersection::Point(p,(a,b,c))` iff the lines intersect in exactly one point `p`.
    /// `f` is a value such that `self.start + self.vector()*a/c == p` and
    /// `other.start + other.vector()*b/c == p`.
    ///
    /// # Examples
    ///
    /// ```
    /// use iron_shapes::point::Point;
    /// use iron_shapes::edge::*;
    ///
    /// let e1 = Edge::new((0, 0), (2, 2));
    /// let e2 = Edge::new((0, 2), (2, 0));
    ///
    /// assert_eq!(e1.line_intersection_rounded(e2),
    ///     LineIntersection::Point(Point::new(1, 1), (4, 4, 8)));
    ///
    /// ```
    pub fn line_intersection_rounded(&self, other: Edge<T>) -> LineIntersection<T, T> {
        if self.is_degenerate() || other.is_degenerate() {
            LineIntersection::None
        } else {

            // TODO: faster implementation if both lines are orthogonal

            let ab = self.vector();
            let cd = other.vector();

            // Assert that the vectors have a non-zero length. This should already be the case
            // because the degenerate cases are handled before.
            debug_assert!(!ab.is_zero());
            debug_assert!(!cd.is_zero());

            // if ab.x.is_zero() {
            //     // Self is vertical.
            //     if cd.y.is_zero() {
            //         // Lines are orthogonal.
            //         // Get intersection point.
            //         let p = Point::new(self.start.x, other.start.y);
            //         unimplemented!()
            //         // TODO:
            //     } else if cd.x.is_zero() {
            //         // Lines are parallel.
            //         return if self.x == other.x {
            //             // Lines are collinear.
            //             LineIntersection::Collinear
            //         } else {
            //             LineIntersection::None
            //         }
            //     }
            // } else if ab.y.is_zero() {
            //     if cd.x.is_zero() {
            //         // Lines are orthogonal.
            //         // Get intersection point.
            //         let p = Point::new(other.start.x, start.start.y);
            //         unimplemented!()
            //         // TODO:
            //     } else if cd.y.is_zero() {
            //         // Lines are parallel.
            //         return if self.y == other.y {
            //             // Lines are collinear.
            //             LineIntersection::Collinear
            //         } else {
            //             LineIntersection::None
            //         }
            //     }
            // }

            let s = ab.cross_prod(cd);

            if s.is_zero() {
                // Lines are parallel
                debug_assert!(self.is_parallel(&other));

                // TODO: check more efficiently for collinear lines.
                if self.line_contains_point(other.start) {
                    // If the line defined by `self` contains at least one point of `other` then they are equal.
                    debug_assert!(self.is_collinear(&other));
                    LineIntersection::Collinear
                } else {
                    LineIntersection::None
                }
            } else {
                let ac = other.start - self.start;
                let ac_cross_cd = ac.cross_prod(cd);

                // let two = T::one() + T::one();
                // let one_vector = Vector::new(T::one(), T::one());

                // Compute exact solution but scaled by s.
                let exact_scaled_s = self.start * s + ab * ac_cross_cd;

                // Divide by s and round by truncating towards zero.
                // Where the result of an integer division is negative it is truncated towards zero.
                // let p: Point<T> = (exact_scaled_s * two + one_vector * s) / (s * two); // Round to next integer.
                let p: Point<T> = exact_scaled_s / s;

                // TODO: maybe remove computation of relative positions?
                let ca_cross_ab = ac.cross_prod(ab);
                debug_assert!({
                    let exact_scaled_s = other.start * s + cd * ca_cross_ab;
                    // let p2: Point<T> = (exact_scaled_s * two + one_vector * s) / (s * two);// Round to next integer.
                    let p2: Point<T> = exact_scaled_s / s;

                    p == p2
                }
                );

                let positions = if s < T::zero() {
                    (T::zero() - ac_cross_cd, T::zero() - ca_cross_ab, T::zero() - s)
                } else {
                    (ac_cross_cd, ca_cross_ab, s)
                };

                LineIntersection::Point(p, positions)
            }
        }
    }


    /// Compute the intersection with another edge.
    /// Coordinates of intersection points are rounded towards zero.
    ///
    /// `EdgeIntersection::EndPoint` is returned if and only if the intersection lies exactly on an end point.
    pub fn edge_intersection_rounded(&self, other: &Edge<T>) -> EdgeIntersection<T, T, Edge<T>> {
        // Swap direction of other edge such that both have the same direction.
        let other = if (self.start < self.end) != (other.start < other.end) {
            other.reversed()
        } else {
            *other
        };
        debug_assert_eq!(self.start < self.end, other.start < other.end,
                         "Edges should have the same orientation now.");

        // Try to convert the edges into rectilinear edges.
        if let Ok(a) = REdge::try_from(self) {
            if let Ok(b) = REdge::try_from(&other) {

                return match a.edge_intersection(&b) {
                    REdgeIntersection::None => EdgeIntersection::None,
                    REdgeIntersection::EndPoint(p) => {
                        debug_assert!(p == a.start() || p == a.end() || p == b.start() || p == b.end());
                        EdgeIntersection::EndPoint(p)
                    },
                    REdgeIntersection::Point(p) => EdgeIntersection::Point(p),
                    REdgeIntersection::Overlap(e) => EdgeIntersection::Overlap(e.into()),
                }
            }
        }

        // Check endpoints for coincidence.
        // This must be handled separately because equality of the intersection point and endpoints
        // will not necessarily be detected due to rounding errors.
        let same_start_start = self.start == other.start;
        let same_start_end = self.start == other.end;

        let same_end_start = self.end == other.start;
        let same_end_end = self.end == other.end;

        // Are the edges equal but not degenerate?
        let fully_coincident = (same_start_start & same_end_end) ^ (same_start_end & same_end_start);

        let result = if self.is_degenerate() {
            // First degenerate case
            if other.contains_point(self.start).inclusive_bounds() {
                EdgeIntersection::EndPoint(self.start)
            } else {
                EdgeIntersection::None
            }
        } else if other.is_degenerate() {
            // Second degenerate case
            if self.contains_point(other.start).inclusive_bounds() {
                EdgeIntersection::EndPoint(other.start)
            } else {
                EdgeIntersection::None
            }
        } else if fully_coincident {
            EdgeIntersection::Overlap(*self)
        } else if !self.bounding_box().touches(&other.bounding_box()) {
            // If bounding boxes do not touch, then intersection is impossible.
            EdgeIntersection::None
        } else {
            // Compute the intersection of the lines defined by the two edges.
            let line_intersection = self.line_intersection_rounded(other);

            // Then check if the intersection point is on both edges
            // or find the intersection if the edges overlap.
            match line_intersection {
                LineIntersection::None => EdgeIntersection::None,

                LineIntersection::Point(p, (pos1, pos2, len)) => {
                    if pos1 >= T::zero() && pos1 <= len
                        && pos2 >= T::zero() && pos2 <= len {
                        if pos1 == T::zero()
                            || pos1 == len
                            || pos2 == T::zero()
                            || pos2 == len {
                            EdgeIntersection::EndPoint(p)
                        } else {
                            EdgeIntersection::Point(p)
                        }
                    } else {
                        EdgeIntersection::None
                    }
                }
                LineIntersection::Collinear => {
                    debug_assert!(self.is_collinear(&other));

                    // Project all points of the two edges on the line defined by the first edge
                    // (scaled by the length of the first edge).
                    // This allows to calculate the interval of overlap in one dimension.

                    let (pa, pb) = self.into();
                    let (pc, pd) = other.into();

                    let b = pb - pa;
                    let c = pc - pa;
                    let d = pd - pa;

                    let dist_a = T::zero();
                    let dist_b = b.dot(b);

                    let dist_c = b.dot(c);
                    let dist_d = b.dot(d);

                    let start1 = (dist_a, pa);
                    let end1 = (dist_b, pb);

                    // Sort end points of other edge.
                    let (start2, end2) = if dist_c < dist_d {
                        ((dist_c, pc), (dist_d, pd))
                    } else {
                        ((dist_d, pd), (dist_c, pc))
                    };

                    // Find maximum by distance.
                    let start = if start1.0 < start2.0 {
                        start2
                    } else {
                        start1
                    };

                    // Find minimum by distance.
                    let end = if end1.0 < end2.0 {
                        end1
                    } else {
                        end2
                    };

                    // Check if the edges overlap in more than one point, in exactly one point or
                    // in zero points.
                    match start.0.cmp(&end.0) {
                        Ordering::Less => EdgeIntersection::Overlap(Edge::new(start.1, end.1)),
                        Ordering::Equal => EdgeIntersection::EndPoint(start.1),
                        Ordering::Greater => EdgeIntersection::None
                    }
                }
            }
        };

        // Check that the result is consistent with the edge intersection test.
        debug_assert_eq!(
            result == EdgeIntersection::None,
            self.edges_intersect(&other).is_no()
        );

        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_line_intersection_rounded() {
        let e1 = Edge::new((0, 0), (2, 2));
        let e2 = Edge::new((0, 2), (2, 0));
        let intersection = e1.line_intersection_rounded(e2);
        assert_eq!(intersection, LineIntersection::Point((1, 1).into(), (4, 4, 8)));


        let e1 = Edge::new((0, 1), (1, 1));
        let e2 = Edge::new((0, 0), (1, 2));
        let intersection = e1.line_intersection_rounded(e2);
        match intersection {
            LineIntersection::Point(p, _) => assert_eq!(p, Point::new(0, 1)),
            _ => assert!(false)
        }

        let e1 = Edge::new((0, 1), (1, 1));
        let e2 = Edge::new((0, 0), (1, 3));
        let intersection = e1.line_intersection_rounded(e2);
        match intersection {
            LineIntersection::Point(p, _) => assert_eq!(p, Point::new(0, 1)),
            _ => assert!(false)
        }

        // Test truncation towards zero.
        let e1 = Edge::new((0i32, 0), (1, 0));
        let e2 = Edge::new((0, 1), (12345677, -1));
        let intersection = e1.line_intersection_rounded(e2);
        match intersection {
            LineIntersection::Point(p, _) => assert_eq!(p, Point::new(12345677 / 2, 0)),
            _ => assert!(false)
        }

        // Test truncation towards zero.
        let e1 = Edge::new((0i32, 0), (-1, 0));
        let e2 = Edge::new((0, -1), (-12345677, 1));
        let intersection = e1.line_intersection_rounded(e2);
        match intersection {
            LineIntersection::Point(p, _) => assert_eq!(p, Point::new(-12345677 / 2, 0)),
            _ => assert!(false)
        }
    }

    #[test]
    fn test_edge_intersection_rounded() {

        // Intersection inside the edges.
        let e1 = Edge::new((-10, 0), (10, 0));
        let e2 = Edge::new((-1, -1), (2, 2));
        assert_eq!(e1.edge_intersection_rounded(&e2), EdgeIntersection::Point((0, 0).into()));
        assert_eq!(e2.edge_intersection_rounded(&e1), EdgeIntersection::Point((0, 0).into()));

        // Intersection on an endpoint.
        let e1 = Edge::new((-10, 0), (10, 0));
        let e2 = Edge::new((0, 0), (2, 2));
        assert_eq!(e1.edge_intersection_rounded(&e2), EdgeIntersection::EndPoint((0, 0).into()));
        assert_eq!(e2.edge_intersection_rounded(&e1), EdgeIntersection::EndPoint((0, 0).into()));


        // Intersection on both endpoint.
        let e1 = Edge::new((0, 0), (10, 0));
        let e2 = Edge::new((0, 0), (2, 2));
        assert_eq!(e1.edge_intersection_rounded(&e2), EdgeIntersection::EndPoint((0, 0).into()));
        assert_eq!(e2.edge_intersection_rounded(&e1), EdgeIntersection::EndPoint((0, 0).into()));


        // Intersection not on an endpoint but rounded down to an endpoint.
        // TODO: Rethink what should happen here. EndPoint or Point?
        let e1 = Edge::new((0, 0), (1, 0));
        let e2 = Edge::new((0, -1), (1, 10));
        assert_eq!(e1.edge_intersection_rounded(&e2), EdgeIntersection::Point((0, 0).into()));
        assert_eq!(e2.edge_intersection_rounded(&e1), EdgeIntersection::Point((0, 0).into()));

        // No intersection.
        let e1 = Edge::new((-10, 0), (10, 0));
        let e2 = Edge::new((1, 1), (2, 2));
        assert_eq!(e1.edge_intersection_rounded(&e2), EdgeIntersection::None);
        assert_eq!(e2.edge_intersection_rounded(&e1), EdgeIntersection::None);
    }

    #[test]
    fn test_end_point_intersection_at_negative_x() {
        let p = Point::new;

        // Negative coordinates.
        let e1 = Edge::new(p(-1, 2), p(0, 0));
        let e2 = Edge::new(p(-1, 2), p(0, 2));

        assert_eq!(e1.edge_intersection_rounded(&e2),
                   EdgeIntersection::EndPoint(p(-1, 2)));
    }
}