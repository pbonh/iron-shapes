// Copyright (c) 2018-2020 Thomas Kramer.
// SPDX-FileCopyrightText: 2018-2022 Thomas Kramer
//
// SPDX-License-Identifier: AGPL-3.0-or-later

//! Edge intersection functions for rational coordinates.
//!
//! Some computations on edges can be performed exactly when the coordinates
//! have rational type. For example the intersection of two edges with rational coordinates
//! can be expressed exactly as a point with rational coordinates.
//!
//! This modules contains implementations of such computations that take edges with rational coordinates
//! and produce an output in rational coordinates.
//!

use crate::point::Point;

use num_rational::Ratio;
use num_traits::Zero;

use std::cmp::Ordering;

pub use crate::edge::{Edge, EdgeIntersection, LineIntersection};
use crate::traits::BoundingBox;
use crate::CoordinateType;

impl<T: CoordinateType + num_integer::Integer> Edge<Ratio<T>> {
    /// Compute the intersection point of the lines defined by the two edges.
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
    /// extern crate num_rational;
    /// use num_rational::Ratio;
    /// use iron_shapes::point::Point;
    /// use iron_shapes::edge_rational::*;
    ///
    /// let r = |i| Ratio::from_integer(i);
    ///
    /// let e1 = Edge::new((r(0), r(0)), (r(2), r(2)));
    /// let e2 = Edge::new((r(0), r(2)), (r(2), r(0)));
    ///
    /// assert_eq!(e1.line_intersection_rational(e2),
    ///     LineIntersection::Point(Point::new(r(1), r(1)), (r(4), r(4), r(8))));
    ///
    /// ```
    pub fn line_intersection_rational(
        &self,
        other: Edge<Ratio<T>>,
    ) -> LineIntersection<Ratio<T>, Ratio<T>> {
        if self.is_degenerate() || other.is_degenerate() {
            LineIntersection::None
        } else {
            // TODO: faster implementation if both lines are orthogonal

            let ab = self.vector();
            let cd = other.vector();

            // Assert that the vectors have a non-zero length. This should already be the case
            // because the degenerate cases are handled before.
            debug_assert!(ab.norm2_squared() > Ratio::zero());
            debug_assert!(cd.norm2_squared() > Ratio::zero());

            let s = ab.cross_prod(cd);

            // TODO: What if approximate zero due to rounding error?
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
                let i = ac_cross_cd / s;

                let p: Point<Ratio<T>> = self.start + ab * i;

                let ca_cross_ab = ac.cross_prod(ab);

                // Check that the intersection point lies on the lines indeed.
                // TODO: uncomment this checks
                //                debug_assert!(self.cast().line_contains_point_approx(p, 1e-4));
                //                debug_assert!(other.cast().line_contains_point_approx(p, 1e-2));

                // debug_assert!({
                //     let j = ca_cross_ab / s;
                //     let p2: Point<Ratio<T>> = other.start + cd * j;
                //     (p - p2).norm2_squared() < Ratio::new(T::one(), 1_000_000)
                // });

                let positions = if s < Ratio::zero() {
                    (
                        Ratio::zero() - ac_cross_cd,
                        Ratio::zero() - ca_cross_ab,
                        Ratio::zero() - s,
                    )
                } else {
                    (ac_cross_cd, ca_cross_ab, s)
                };

                LineIntersection::Point(p, positions)
            }
        }
    }

    /// Compute the intersection with another edge.
    pub fn edge_intersection_rational(
        &self,
        other: &Edge<Ratio<T>>,
    ) -> EdgeIntersection<Ratio<T>, Ratio<T>, Edge<Ratio<T>>> {
        //        debug_assert!(tolerance >= Ratio::zero(), "Tolerance cannot be negative.");

        // Swap direction of other edge such that both have the same direction.
        let other = if (self.start < self.end) != (other.start < other.end) {
            other.reversed()
        } else {
            *other
        };

        // Check endpoints for coincidence.
        // This must be handled separately because equality of the intersection point and endpoints
        // will not necessarily be detected due to rounding errors.
        let same_start_start = self.start == other.start;
        let same_start_end = self.start == other.end;

        let same_end_start = self.end == other.start;
        let same_end_end = self.end == other.end;

        // TODO: optimize for chained edges (start1 == end2 ^ start2 == end1)

        // Are the edges equal but not degenerate?
        let fully_coincident =
            (same_start_start & same_end_end) ^ (same_start_end & same_end_start);

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
            let line_intersection = self.line_intersection_rational(other);

            // Then check if the intersection point is on both edges
            // or find the intersection if the edges overlap.
            match line_intersection {
                LineIntersection::None => EdgeIntersection::None,

                // Intersection in one point:
                LineIntersection::Point(p, (pos1, pos2, len)) => {
                    if pos1 >= Ratio::zero() && pos1 <= len && pos2 >= Ratio::zero() && pos2 <= len
                    {
                        if pos1 == Ratio::zero()
                            || pos1 == len
                            || pos2 == Ratio::zero()
                            || pos2 == len
                        {
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

                    let dist_a = Ratio::zero();
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
                    let start = if start1.0 < start2.0 { start2 } else { start1 };

                    // Find minimum by distance.
                    let end = if end1.0 < end2.0 { end1 } else { end2 };

                    // Check if the edges overlap in more than one point, in exactly one point or
                    // in zero points.
                    match start.0.cmp(&end.0) {
                        Ordering::Less => EdgeIntersection::Overlap(Edge::new(start.1, end.1)),
                        Ordering::Equal => EdgeIntersection::EndPoint(start.1),
                        Ordering::Greater => EdgeIntersection::None,
                    }
                }
            }
        };

        // Sanity check for the result.
        debug_assert!({
            match result {
                EdgeIntersection::Point(p) => {
                    self.contains_point(p).is_within_bounds()
                        && other.contains_point(p).is_within_bounds()
                }
                EdgeIntersection::EndPoint(p) => {
                    self.contains_point(p).on_bounds() || other.contains_point(p).on_bounds()
                }
                EdgeIntersection::None => self.edges_intersect(&other).is_no(),
                EdgeIntersection::Overlap(_) => true,
            }
        });

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
    use crate::traits::Scale;
    use num_rational::Rational64;

    #[test]
    fn test_rational_edge() {
        let e = Edge::new(
            (Rational64::from(0), Rational64::from(0)),
            (Rational64::from(1), Rational64::from(1)),
        );
        let v = e.vector();
        assert!(v.norm2_squared() == Rational64::from(2));
        assert!(!e.is_degenerate());
    }

    #[test]
    fn test_line_intersection_rational() {
        // Helper constructors.
        let rp = |a: i64, b: i64| Point::new(Rational64::from(a), Rational64::from(b));
        let r = Rational64::new;
        let i = Rational64::from;

        let e1 = Edge::new(rp(0, 0), rp(2, 2));
        let e2 = Edge::new(rp(1, 0), rp(0, 1));
        let e3 = Edge::new(rp(1, 0), rp(3, 2));

        assert_eq!(
            e1.line_intersection_rational(e2),
            LineIntersection::Point(Point::new(r(1, 2), r(1, 2)), (i(1), i(2), i(4)))
        );

        // Parallel lines should not intersect
        assert_eq!(e1.line_intersection_rational(e3), LineIntersection::None);

        let e4 = Edge::new(rp(-320, 2394), rp(94, -4482));
        let e5 = Edge::new(rp(71, 133), rp(-1373, 13847));

        if let LineIntersection::Point(intersection, _) = e4.line_intersection_rational(e5) {
            let a = e4.vector();
            let b = intersection - e4.start;
            let area = b.cross_prod(a);
            assert!(area.is_zero());

            let a = e5.vector();
            let b = intersection - e5.start;
            let area = b.cross_prod(a);
            assert!(area.is_zero());
        } else {
            assert!(false);
        }

        // Collinear lines.
        let scale = Rational64::new(1, 3);
        let e1 = Edge::new(rp(0, 0), rp(2, 2)).scale(scale);
        let e2 = Edge::new(rp(4, 4), rp(8, 8)).scale(scale);
        assert!(!e1.is_coincident(&e2));
        assert!(e1.is_parallel(&e2));
        assert_eq!(
            e1.line_intersection_rational(e2),
            LineIntersection::Collinear
        );
    }

    #[test]
    fn test_edge_intersection_rational() {}
}
