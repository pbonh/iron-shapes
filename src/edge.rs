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
use crate::point::Point;
use crate::rect::Rect;

use crate::CoordinateType;

use num_traits::Float;
use num_traits::cast::NumCast;

use crate::traits::Transform;
pub use crate::traits::{BoundingBox, RotateOrtho, Angle};

pub use crate::types::{Side, ContainsResult};


#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum EdgeIntersection<TP: CoordinateType, TO: CoordinateType> {
    /// No intersection.
    None,
    /// Intersection in a single point but not on an endpoint of an edge.
    Point(Point<TP>),
    /// Intersection in an endpoint of an edge.
    EndPoint(Point<TO>),
    /// Full or partial overlap.
    Overlap(Edge<TO>),
}

/// Store all possible results of a line to line intersection.
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum LineIntersection<TP: CoordinateType, TE: CoordinateType> {
    /// No intersection at all.
    None,
    /// Intersection in a single point.
    /// Besides the intersection point also an other expression for the intersection point is given.
    /// The three values `(a, b, c)` describe the intersection point in terms of a starting point (the starting point
    /// of the edge which defines the line) and the direction of the edge multiplied by a fraction.
    ///
    /// `edge.start + edge.vector()*a/c == p` and
    /// `other_edge.start + other_edge.vector()*b/c == p`.
    Point(Point<TP>, (TE, TE, TE)),
    /// Lines are collinear.
    Collinear,
}

impl<T: CoordinateType> Into<(Point<T>, Point<T>)> for Edge<T> {
    fn into(self) -> (Point<T>, Point<T>) {
        (self.start, self.end)
    }
}

impl<T: CoordinateType> Into<(Point<T>, Point<T>)> for &Edge<T> {
    fn into(self) -> (Point<T>, Point<T>) {
        (self.start, self.end)
    }
}

impl<T: CoordinateType> From<(Point<T>, Point<T>)> for Edge<T> {
    fn from(points: (Point<T>, Point<T>)) -> Self {
        Edge::new(points.0, points.1)
    }
}

impl<T: CoordinateType> From<[Point<T>; 2]> for Edge<T> {
    fn from(points: [Point<T>; 2]) -> Self {
        Edge::new(points[0], points[1])
    }
}

/// An edge (line segment) is represented by its starting point and end point.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct Edge<T: Copy> {
    pub start: Point<T>,
    pub end: Point<T>,
}

impl<T: CoordinateType> Edge<T> {
    /// Create a new `Edge` from two arguments that implement `Into<Point>`.
    pub fn new<C>(start: C, end: C) -> Self
        where C: Into<Point<T>> {
        Edge {
            start: start.into(),
            end: end.into(),
        }
    }

    /// Return the same edge but with the two points swapped.
    pub fn reversed(&self) -> Self {
        Edge {
            start: self.end,
            end: self.start,
        }
    }

    /// Check if edge is degenerate.
    /// An edge is degenerate if start point and end point are equal.
    pub fn is_degenerate(&self) -> bool {
        self.start == self.end
    }

    /// Test if this edge is either horizontal or vertical.
    pub fn is_ortho(&self) -> bool {
        !self.is_degenerate() &&
            (self.start.x == self.end.x || self.start.y == self.end.y)
    }

    /// Test if this edge is horizontal.
    pub fn is_horizontal(&self) -> bool {
        !self.is_degenerate() &&
            self.start.y == self.end.y
    }

    /// Test if this edge is vertical.
    pub fn is_vertical(&self) -> bool {
        !self.is_degenerate() &&
            self.start.x == self.end.x
    }


    /// Returns the vector from `self.start` to `self.end`.
    pub fn vector(&self) -> Vector<T> {
        self.end - self.start
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
        assert!(!self.is_degenerate());

        let a = self.vector();
        let b = point - self.start;

        let area = b.cross_prod(a);

        if area.is_zero() {
            Side::Center
        } else if area < T::zero() {
            Side::Left
        } else {
            Side::Right
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
        if self.start == point || self.end == point {
            ContainsResult::OnBounds
        } else if self.is_degenerate() {
            ContainsResult::No
        } else {
            let a = self.start - point;
            let b = self.end - point;
            // Check if the triangle point-start-end has zero area and that a and b have opposite directions.
            // If a and b have opposite directions then point is between start and end.
            if a.cross_prod(b).is_zero() && a.dot(b) <= T::zero() {
                ContainsResult::WithinBounds
            } else {
                ContainsResult::No
            }
        }
    }

    /// Test if point lies on the line defined by the edge.
    pub fn line_contains_point(&self, point: Point<T>) -> bool {
        if self.is_degenerate() {
            self.start == point
        } else {
            let l = self.vector();
            let b = point - self.start;

            l.cross_prod(b).is_zero()
        }
    }


    /// Test if two edges are parallel.
    pub fn is_parallel(&self, other: &Edge<T>) -> bool {
        if self.is_degenerate() || other.is_degenerate() {
            false
        } else {
            let a = self.vector();
            let b = other.vector();

            a.cross_prod(b).is_zero()
        }
    }


    /// Test if two edges are collinear, i.e. are on the same line.
    pub fn is_collinear(&self, other: &Edge<T>) -> bool
        where T: CoordinateType {
        if self.is_degenerate() || other.is_degenerate() {
            false
        } else {
            let v = self.vector();
            let a = other.start - self.start;
            let b = other.end - self.start;

            v.cross_prod(a).is_zero() &&
                v.cross_prod(b).is_zero()
        }
    }

    /// Test edges for coincidence.
    /// Two edges are coincident if they are oriented the same way
    /// and share more than one point (implies that they must be parallel).
    pub fn is_coincident(&self, other: &Edge<T>) -> bool {

        // TODO: use is_collinear
        // TODO: approximate version for floating point types

        let self_start_in_other = other.contains_point(self.start).inclusive_bounds();
        let self_end_in_other = other.contains_point(self.end).inclusive_bounds();
        let other_start_in_self = self.contains_point(other.start).inclusive_bounds();
        let other_end_in_self = self.contains_point(other.end).inclusive_bounds();

        let share_more_than_one_point =
            self.end != other.start && self.start != other.end && (
                other_start_in_self && self_end_in_other ||
                    other_end_in_self && self_start_in_other ||
                    other_start_in_self && other_end_in_self ||
                    self_start_in_other && self_end_in_other);

        // Sharing more than one point should imply that the edges are parallel.
        debug_assert!(if share_more_than_one_point { self.is_parallel(other) } else { true });

        let oriented_the_same_way = self.vector().dot(other.vector()) > T::zero();

        share_more_than_one_point && oriented_the_same_way
    }


    /// Test if two edges are approximately parallel.
    /// To be used for float coordinates.
    /// Inspired by algorithm on page 241 of "Geometric Tools for Computer Graphics".
    pub fn is_parallel_approx(&self, other: &Edge<T>, epsilon_squared: T) -> bool
    {
        if self.is_degenerate() || other.is_degenerate() {
            false
        } else {
            let d1 = self.vector();
            let d2 = other.vector();

            let len1_sqr = d1.norm2_squared();
            let len2_sqr = d2.norm2_squared();

            let cross = d1.cross_prod(d2);
            let cross_sqr = cross * cross;

            // TODO: require square tolerance form caller?
            cross_sqr <= len1_sqr * len2_sqr * epsilon_squared
        }
    }


    /// Test if two edges are approximately collinear, i.e. are on the same line.
    /// Inspired by algorithm on page 241 of "Geometric Tools for Computer Graphics".
    pub fn is_collinear_approx(&self, other: &Edge<T>, epsilon_squared: T) -> bool {
        if self.is_degenerate() || other.is_degenerate() {
            false
        } else {
            let d1 = self.vector();
            let d2 = other.vector();

            let len1_sqr = d1.norm2_squared();
            let len2_sqr = d2.norm2_squared();

            let cross = d1.cross_prod(d2);
            let cross_sqr = cross * cross;

            let approx_parallel = cross_sqr <= len1_sqr * len2_sqr * epsilon_squared;

            if approx_parallel {
                let e = other.start - self.start;
                let len_e_sqrt = e.norm2_squared();
                let cross = e.cross_prod(d1);
                let cross_sqr = cross * cross;
                cross_sqr <= len1_sqr * len_e_sqrt * epsilon_squared
            } else {
                false
            }
        }
    }


    /// Test if lines defined by the edges intersect.
    /// If the lines are collinear they are also considered intersecting.
    pub fn lines_intersect_approx(&self, other: &Edge<T>, epsilon_squared: T) -> bool {
        !self.is_parallel_approx(other, epsilon_squared) || self.is_collinear_approx(other, epsilon_squared)
    }

    /// Test if this edge is crossed by the line defined by the other edge.
    ///
    /// Returns `WithinBounds` if start and end point of this edge lie on different sides
    /// of the line defined by the `other` edge or `OnBounds` if at least one of the points
    /// lies on the line.
    pub fn crossed_by_line(&self, other: &Edge<T>) -> ContainsResult {
        // TODO: Handle degenerate cases.
        let side1 = other.side_of(self.start);

        if side1 == Side::Center {
            ContainsResult::OnBounds
        } else {
            let side2 = other.side_of(self.end);

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
    pub fn lines_intersect(&self, other: &Edge<T>) -> bool {
        !self.is_parallel(other) || self.is_collinear(other)
    }

    /// Test if two edges intersect.
    /// If the edges coincide, they also intersect.
    pub fn edges_intersect(&self, other: &Edge<T>) -> ContainsResult {

        // Two edges intersect if the start and end point of one edge
        // lie on opposite sides of the other edge.

        if self.is_degenerate() {
            other.contains_point(self.start)
        } else if other.is_degenerate() {
            self.contains_point(other.start)
        } else if !self.bounding_box().touches(&other.bounding_box()) {
            ContainsResult::No
        } else {
            // TODO:
//        else if self.is_ortho() && other.is_ortho() {
//            // We know now that the bounding boxes touch each other.
//            // For rectilinear edges this implies that they touch somewhere or overlap.
//            true
//        } else {
            match self.crossed_by_line(&other) {
                ContainsResult::No => ContainsResult::No,
                r => r.min(other.crossed_by_line(self))
            }
        }
    }
}


impl<T: CoordinateType + NumCast> Edge<T> {
    /// Test if point lies on the line defined by the edge.
    pub fn line_contains_point_approx<F: Float + NumCast>(&self, point: Point<T>, tolerance: F) -> bool {
        if self.is_degenerate() {
            self.start == point
        } else {
            self.distance_to_line_abs_approx::<F>(point) <= tolerance
        }
    }

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
    /// use iron_shapes::point::Point;
    /// use iron_shapes::edge::*;
    ///
    /// let e1 = Edge::new((0, 0), (2, 2));
    /// let e2 = Edge::new((0, 2), (2, 0));
    ///
    /// assert_eq!(e1.line_intersection_approx(&e2, 1e-6),
    ///     LineIntersection::Point(Point::new(1., 1.), (4, 4, 8)));
    ///
    /// assert_eq!(e1.vector().cast() * 0.5, Point::new(1., 1.));
    /// ```
    ///
    pub fn line_intersection_approx<F: Float>(&self, other: &Edge<T>, tolerance: F) -> LineIntersection<F, T> {
        // TODO: implement algorithm on page 241 of Geometric Tools for Computer Graphics.

        debug_assert!(tolerance >= F::zero(), "Tolerance cannot be negative.");

        if self.is_degenerate() {
            LineIntersection::None
        } else if other.is_degenerate() {
            LineIntersection::None
        } else {

            // TODO: faster implementation if both lines are orthogonal

            let ab = self.vector();
            let cd = other.vector();

            // Assert that the vectors have a non-zero length. This should already be the case
            // because the degenerate cases are handled before.
            debug_assert!(ab.norm2_squared() > T::zero());
            debug_assert!(cd.norm2_squared() > T::zero());

            let s = ab.cross_prod(cd);
            let s_float = F::from(s).unwrap();

            // TODO: What if approximate zero due to rounding error?
            if s_float.abs() <= tolerance * ab.length() * cd.length() {
                // Lines are parallel
                // TODO use assertion
//                debug_assert!(self.is_parallel_approx(other, tolerance));

                // TODO: check more efficiently for collinear lines.
                if self.line_contains_point_approx(other.start, tolerance) {
                    // If the line defined by `self` contains at least one point of `other` then they are equal.
                    // TODO use assertion
//                    debug_assert!(self.is_collinear_approx(other, tolerance));
                    LineIntersection::Collinear
                } else {
                    LineIntersection::None
                }
            } else {
                let ac = other.start - self.start;
                let ac_cross_cd = ac.cross_prod(cd);
                let i = F::from(ac_cross_cd).unwrap() / s_float;

                let p: Point<F> = self.start.cast() + ab.cast() * i;

                let ca_cross_ab = ac.cross_prod(ab);

                // Check that the intersection point lies on the lines indeed.
                debug_assert!(self.cast().line_contains_point_approx(p, tolerance + tolerance));
                debug_assert!(other.cast().line_contains_point_approx(p, tolerance + tolerance));

                debug_assert!({
                    let j = F::from(ca_cross_ab).unwrap() / s_float;
                    let p2: Point<F> = other.start.cast() + cd.cast() * j;
                    (p - p2).norm2_squared() < tolerance + tolerance
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
    pub fn edge_intersection_approx<F: Float>(&self, other: &Edge<T>, tolerance: F) -> EdgeIntersection<F, T> {
        debug_assert!(tolerance >= F::zero(), "Tolerance cannot be negative.");

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

        // Are the edges equal but not degenerate?
        let fully_coincident = (same_start_start & same_end_end) ^ (same_start_end & same_end_start);

        let result = if self.is_degenerate() {
            // First degenerate case
            if other.contains_point_approx(self.start, tolerance) {
                EdgeIntersection::EndPoint(self.start)
            } else {
                EdgeIntersection::None
            }
        } else if other.is_degenerate() {
            // Second degenerate case
            if self.contains_point_approx(other.start, tolerance) {
                EdgeIntersection::EndPoint(other.start)
            } else {
                EdgeIntersection::None
            }
        } else if fully_coincident {
            EdgeIntersection::Overlap(*self)
        } else if !self.bounding_box().touches(&other.bounding_box()) { // TODO: Add tolerance here.
            // If bounding boxes do not touch, then intersection is impossible.
            EdgeIntersection::None
        } else {
            // Compute the intersection of the lines defined by the two edges.
            let line_intersection = self.line_intersection_approx(&other, tolerance);

            // Then check if the intersection point is on both edges
            // or find the intersection if the edges overlap.

            match line_intersection {
                LineIntersection::None => EdgeIntersection::None,

                // Coincident at an endpoint:
                LineIntersection::Point(_, _) if same_start_start || same_start_end =>
                    EdgeIntersection::EndPoint(self.start),

                // Coincident at an endpoint:
                LineIntersection::Point(_, _) if same_end_start || same_end_end =>
                    EdgeIntersection::EndPoint(self.end),

                // Intersection in one point:
                LineIntersection::Point(p, (pos1, pos2, len)) => {
                    let len_f = F::from(len).unwrap();
                    let zero_tol = -tolerance;
                    let len_tol = len_f + tolerance;

                    let pos1 = F::from(pos1).unwrap();
                    let pos2 = F::from(pos2).unwrap();

                    // Check if the intersection point `p` lies on the edge.
                    if pos1 >= zero_tol && pos1 <= len_tol
                        && pos2 >= zero_tol && pos2 <= len_tol {
                        // Intersection

                        // TODO: rework

                        let zero_tol = tolerance;
                        let len_tol = len_f - tolerance;

                        if pos1 <= zero_tol {
                            EdgeIntersection::EndPoint(self.start)
                        } else if pos1 >= len_tol {
                            EdgeIntersection::EndPoint(self.end)
                        } else if pos2 <= zero_tol {
                            EdgeIntersection::EndPoint(other.start)
                        } else if pos2 >= len_tol {
                            EdgeIntersection::EndPoint(other.end)
                        } else {
                            // Intersection in-between the endpoints.
                            debug_assert!({
                                              let e1 = self.cast();
                                              let e2 = other.cast();
                                              p != e1.start
                                                  && p != e1.end
                                                  && p != e2.start
                                                  && p != e2.end
                                          },
                                          "Intersection should not be an endpoint.");
                            EdgeIntersection::Point(p)
                        }
                    } else {
                        // No intersection.
                        EdgeIntersection::None
                    }
                }
                LineIntersection::Collinear => {
                    // TODO use assertion
//                    debug_assert!(self.is_collinear_approx(other, tolerance));

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
                    if start.0 < end.0 {
                        EdgeIntersection::Overlap(Edge::new(start.1, end.1))
                    } else if start.0 == end.0 {
                        EdgeIntersection::EndPoint(start.1)
                    } else {
                        EdgeIntersection::None
                    }
                }
            }
        };

        // Check that the result is consistent with the edge intersection test.
//        debug_assert_eq!(
//            result != EdgeIntersection::None,
//            self.edges_intersect(other)
//        );

        debug_assert!(
            if self.edges_intersect(&other).inclusive_bounds() {
                result != EdgeIntersection::None
            } else {
                true
            }
        );

        result
    }
}


impl<T: CoordinateType + NumCast> Edge<T> {
    /// Try to cast into other data type.
    /// When the conversion fails `None` is returned.
    pub fn try_cast<Target: NumCast>(&self) -> Option<Edge<Target>>
        where Target: CoordinateType {
        match (self.start.try_cast(), self.end.try_cast()) {
            (Some(start), Some(end)) => Some(Edge { start, end }),
            _ => None
        }
    }

    /// Cast to other data type.
    ///
    /// # Panics
    /// Panics when the conversion fails.
    pub fn cast<Target>(&self) -> Edge<Target>
        where Target: CoordinateType + NumCast {
        Edge {
            start: self.start.cast(),
            end: self.end.cast(),
        }
    }

    /// Cast to float.
    ///
    /// # Panics
    /// Panics when the conversion fails.
    pub fn cast_to_float<Target>(&self) -> Edge<Target>
        where Target: CoordinateType + NumCast + Float {
        Edge {
            start: self.start.cast_to_float(),
            end: self.end.cast_to_float(),
        }
    }

    /// Calculate the distance from the point to the line given by the edge.
    ///
    /// Distance will be positive if the point lies on the right side of the edge and negative
    /// if the point is on the left side.
    pub fn distance_to_line<F: Float>(&self, point: Point<T>) -> F {
        assert!(!self.is_degenerate());

        let a = self.vector();
        let b = point - self.start;

        let area = b.cross_prod(a);

        let distance = F::from(area).unwrap() / a.length();
        distance
    }

    /// Calculate distance from point to the edge.
    pub fn distance<F: Float>(&self, point: Point<T>) -> F {
        let dist_ortho: F = self.distance_to_line_abs_approx(point);
        let dist_a: F = (point - self.start).length();
        let dist_b = (point - self.end).length();

        dist_ortho.max(dist_a.min(dist_b))
    }

    /// Find the perpendicular projection of a point onto the line of the edge.
    pub fn projection_approx<F: Float>(&self, point: Point<T>) -> Point<F> {
        assert!(!self.is_degenerate());

        let p = (point - self.start).cast();

        let d = (self.vector()).cast();

        // Orthogonal projection of point onto the line.
        let projection = self.start.cast() + d * d.dot(p) / d.norm2_squared();
        projection
    }

    /// Find the mirror image of `point`.
    pub fn reflection_approx<F: Float>(&self, point: Point<T>) -> Point<F> {
        let proj = self.projection_approx(point);
        proj + proj - point.cast()
    }

    /// Calculate the absolute distance from the point onto the unbounded line coincident with this edge.
    pub fn distance_to_line_abs_approx<F: Float>(&self, point: Point<T>) -> F {
        let dist: F = self.distance_to_line(point);
        dist.abs()
    }


    /// Test if point lies approximately on the edge.
    /// Returns true if `point` is up to `tolerance` away from the edge
    /// and lies between start and end points (inclusive).
    pub fn contains_point_approx<F: Float>(&self, point: Point<T>, tolerance: F) -> bool {
        debug_assert!(tolerance >= F::zero());
        if self.is_degenerate() {
            let l = (self.start - point).norm2_squared();
            F::from(l).unwrap() <= tolerance
        } else {
            let p = point - self.start;

            let v = self.vector();
            // Rotate by -pi/4 to get the normal.
            let n: Vector<F> = v.rotate_ortho(Angle::R270).cast() / v.length();
            // Project on to normal
            let o = n.dot(p.cast());

            if o.abs() >= tolerance {
                // Point is not close enough to line.
                false
            } else {
                // Is point between start and end?

                // Project on to line
                let l = F::from(v.dot(p)).unwrap();

                l >= -tolerance && l <= F::from(v.norm2_squared()).unwrap() + tolerance
            }
        }
    }
}

impl<T: CoordinateType> Transform<T> for Edge<T> {
    fn transform<F: Fn(Point<T>) -> Point<T>>(&self, tf: F) -> Self {
        Edge {
            start: tf(self.start),
            end: tf(self.end),
        }
    }
}

impl<T: CoordinateType> BoundingBox<T> for Edge<T> {
    fn bounding_box(&self) -> Rect<T> {
        Rect::new(self.start, self.end)
    }
}


#[cfg(test)]
mod tests {
    extern crate rand;

    use std::f64;
    use crate::edge::Edge;
    use crate::point::Point;
    use crate::types::*;
    use super::*;

    use self::rand::distributions::{Uniform, Distribution, Bernoulli};
    use self::rand::rngs::StdRng;
    use self::rand::SeedableRng;

    #[test]
    fn test_is_parallel() {
        let e1 = Edge::new((1, 2), (3, 4));
        let e2 = Edge::new((1 - 10, 2 - 20), (5 - 10, 6 - 20));
        let e3 = Edge::new((1 - 10, 2 - 20), (5 - 10, 6 - 20 + 1));

        assert!(e1.is_parallel(&e2));
        assert!(!e1.is_parallel(&e3));

        assert!(e1.is_parallel_approx(&e2, 0));
        assert!(!e1.is_parallel_approx(&e3, 0));

        assert!(e1.is_parallel_approx(&e3, 1));
    }

    #[test]
    fn test_is_collinear() {
        let e1 = Edge::new((0, 0), (1, 2));
        let e2 = Edge::new((10, 20), (100, 200));
        assert!(e1.is_collinear(&e2));
        assert!(e2.is_collinear(&e1));
        assert!(e1.is_collinear_approx(&e2, 0));
        assert!(e2.is_collinear_approx(&e1, 0));

        // Not collinear.
        let e1 = Edge::new((0i64, 0), (1, 2));
        let e2 = Edge::new((10, 20), (1000, 2001));
        assert!(!e1.is_collinear(&e2));
        assert!(!e2.is_collinear(&e1));
        assert!(!e1.is_collinear_approx(&e2, 0));
        assert!(!e2.is_collinear_approx(&e1, 0));

        assert!(e1.is_collinear_approx(&e2, 1));
        assert!(e2.is_collinear_approx(&e1, 1));
    }

    #[test]
    fn test_distance_to_line() {
        let e1 = Edge::new((1, 0), (2, 1));
        let p0 = Point::new(2, 0);

        let d: f64 = e1.distance_to_line(p0);
        let diff = (d - f64::sqrt(2.) / 2.).abs();

        assert!(diff < 1e-9);
    }

    #[test]
    fn test_line_contains_point() {
        let e1 = Edge::new((1, 2), (5, 6));
        let p0 = Point::new(1, 2);
        let p1 = Point::new(3, 4);
        let p2 = Point::new(5, 6);
        let p3 = Point::new(6, 7);
        let p4 = Point::new(0, 1);
        let p5 = Point::new(0, 0);

        assert!(e1.line_contains_point(p0));
        assert!(e1.line_contains_point(p1));
        assert!(e1.line_contains_point(p2));
        assert!(e1.line_contains_point(p3));
        assert!(e1.line_contains_point(p4));
        assert!(!e1.line_contains_point(p5));

        let e1 = Edge::new((0, 0), (1, 1));
        assert!(e1.line_contains_point(Point::new(2, 2)));
        assert!(e1.line_contains_point(Point::new(-1, -1)));
    }

    #[test]
    fn test_contains() {
        let e1 = Edge::new((1, 2), (5, 6));
        let p0 = Point::new(1, 2);
        let p1 = Point::new(3, 4);
        let p2 = Point::new(5, 6);
        let p3 = Point::new(0, 0);
        let p4 = Point::new(6, 7);

        assert!(e1.contains_point(p0).inclusive_bounds());
        assert!(e1.contains_point(p1).inclusive_bounds());
        assert!(e1.contains_point(p2).inclusive_bounds());
        assert!(!e1.contains_point(p3).inclusive_bounds());
        assert!(!e1.contains_point(p4).inclusive_bounds());

        let tol = 1e-6;
        assert!(e1.contains_point_approx(p0, tol));
        assert!(e1.contains_point_approx(p1, tol));
        assert!(e1.contains_point_approx(p2, tol));
        assert!(!e1.contains_point_approx(p3, tol));
        assert!(!e1.contains_point_approx(p4, tol));

        let e1 = Edge::new((0, 0), (1, 1));
        let p0 = Point::new(2, 2);
        assert!(!e1.contains_point(p0).inclusive_bounds());
    }

    #[test]
    fn test_projection() {
        let e1 = Edge::new((-6., -5.), (4., 7.));
        let p1 = Point::new(1., 2.);
        let p2 = Point::new(-10., 10.);

        let proj1 = e1.projection_approx(p1);
        let proj2 = e1.projection_approx(p2);

        assert!(e1.contains_point_approx(proj1, PREC_DISTANCE));
        assert!(e1.contains_point_approx(proj2, PREC_DISTANCE));
    }

    #[test]
    fn test_side_of() {
        let e1 = Edge::new((1, 0), (4, 4));
        let p1 = Point::new(-10, 0);
        let p2 = Point::new(10, -10);
        let p3 = Point::new(1, 0);

        assert_eq!(e1.side_of(p1), Side::Left);
        assert_eq!(e1.side_of(p2), Side::Right);
        assert_eq!(e1.side_of(p3), Side::Center);
    }

    #[test]
    fn test_crossed_by() {
        let e1 = Edge::new((1, 0), (4, 4));
        let e2 = Edge::new((1 + 1, 0), (4 + 1, 4));
        let e3 = Edge::new((1, 0), (4, 5));
        let e4 = Edge::new((2, -2), (0, 0));

        // Coincident lines
        assert!(e1.crossed_by_line(&e1).inclusive_bounds());

        // Parallel but not coincident.
        assert!(e1.is_parallel(&e2));
        assert!(!e1.crossed_by_line(&e2).inclusive_bounds());

        // Crossing lines.
        assert!(e1.crossed_by_line(&e3).inclusive_bounds());

        // crossed_by is not commutative
        assert!(!e1.crossed_by_line(&e4).inclusive_bounds());
        assert!(e4.crossed_by_line(&e1).inclusive_bounds());
    }

    #[test]
    fn test_intersect() {
        let e1 = Edge::new((0, 0), (4, 4));
        let e2 = Edge::new((1, 0), (0, 1));
        let e3 = Edge::new((0, -1), (-1, 1));

        assert!(e1.edges_intersect(&e1).inclusive_bounds());
        assert!(e1.edges_intersect(&e2).inclusive_bounds());
        assert!(!e1.crossed_by_line(&e3).inclusive_bounds());
        assert!(e3.crossed_by_line(&e1).inclusive_bounds());
        assert!(!e1.edges_intersect(&e3).inclusive_bounds());

        // Bounding boxes overlap but edges are not touching.
        let e1 = Edge::new((0, 0), (8, 8));
        let e2 = Edge::new((4, 0), (3, 1));
        assert!(!e1.is_ortho());
        assert!(!e2.is_ortho());
        assert!(e1.bounding_box().touches(&e2.bounding_box()));
        assert!(!e1.edges_intersect(&e2).inclusive_bounds());
        assert!(e1.lines_intersect(&e2));

        // Intersection at an endpoint.
        let e1 = Edge::new((0, 0), (4, 4));
        let e2 = Edge::new((1, 1), (2, 0));
        assert_eq!(e1.edges_intersect(&e2), ContainsResult::OnBounds);
        assert_eq!(e2.edges_intersect(&e1), ContainsResult::OnBounds);
    }


    #[test]
    fn test_line_intersection() {
        let tol = 1e-12;
        let e1 = Edge::new((0, 0), (2, 2));
        let e2 = Edge::new((1, 0), (0, 1));
        let e3 = Edge::new((1, 0), (3, 2));

        assert_eq!(e1.line_intersection_approx(&e2, tol),
                   LineIntersection::Point(Point::new(0.5, 0.5), (1, 2, 4)));
        // Parallel lines should not intersect
        assert_eq!(e1.line_intersection_approx(&e3, tol), LineIntersection::None);


        let e4 = Edge::new((-320., 2394.), (94., -448382.));
        let e5 = Edge::new((71., 133.), (-13733., 1384.));

        if let LineIntersection::Point(intersection, _) = e4.line_intersection_approx(&e5, tol) {
            assert!(e4.distance_to_line::<f64>(intersection).abs() < tol);
            assert!(e5.distance_to_line::<f64>(intersection).abs() < tol);
        } else {
            assert!(false);
        }

        // Collinear lines.
        let e1 = Edge::new((0., 0.), (2., 2.));
        let e2 = Edge::new((4., 4.), (8., 8.));
        assert!(!e1.is_coincident(&e2));
        assert!(e1.is_parallel_approx(&e2, tol));
        assert_eq!(e1.line_intersection_approx(&e2, tol), LineIntersection::Collinear);
    }


    #[test]
    fn test_edge_intersection() {
        let tol = 1e-20;
        // Point intersection inside both edges.
        let e1 = Edge::new((0, 0), (2, 2));
        let e2 = Edge::new((2, 0), (0, 2));
        assert!(e1.edges_intersect(&e2).inclusive_bounds());
        assert_eq!(e1.edge_intersection_approx(&e2, tol),
                   EdgeIntersection::Point(Point::new(1f64, 1f64)));
        assert_eq!(e2.edge_intersection_approx(&e1, tol),
                   EdgeIntersection::Point(Point::new(1f64, 1f64)));

        // Point intersection on the end of one edge.
        let e1 = Edge::new((0, 0), (2, 2));
        let e2 = Edge::new((2, 0), (1, 1));
        assert!(e1.edges_intersect(&e2).inclusive_bounds());
        assert_eq!(e1.edge_intersection_approx(&e2, tol),
                   EdgeIntersection::EndPoint(Point::new(1, 1)));
        assert_eq!(e2.edge_intersection_approx(&e1, tol),
                   EdgeIntersection::EndPoint(Point::new(1, 1)));

        // No intersection
        let e1 = Edge::new((0, 0), (4, 4));
        let e2 = Edge::new((3, 0), (2, 1));
        assert!(!e1.edges_intersect(&e2).inclusive_bounds());
        assert_eq!(e1.edge_intersection_approx(&e2, tol),
                   EdgeIntersection::None);
        assert_eq!(e2.edge_intersection_approx(&e1, tol),
                   EdgeIntersection::None);
    }

    #[test]
    fn test_edge_intersection_overlap() {
        let tol = 1e-12;
        // Overlapping edges. Same orientations.
        let e1 = Edge::new((0, 0), (2, 0));
        let e2 = Edge::new((1, 0), (3, 0));
        assert!(e1.edges_intersect(&e2).inclusive_bounds());
        assert!(e1.is_collinear(&e2));

        assert_eq!(e1.edge_intersection_approx(&e2, tol),
                   EdgeIntersection::Overlap(Edge::new((1, 0), (2, 0))));
        assert_eq!(e2.edge_intersection_approx(&e1, tol),
                   EdgeIntersection::Overlap(Edge::new((1, 0), (2, 0))));


        // Overlapping edges. Opposing orientations.
        let e1 = Edge::new((0, 0), (2, 2));
        let e2 = Edge::new((3, 3), (1, 1));
        assert!(e1.edges_intersect(&e2).inclusive_bounds());
        assert!(e1.is_collinear(&e2));

        assert_eq!(e1.edge_intersection_approx(&e2, tol),
                   EdgeIntersection::Overlap(Edge::new((1, 1), (2, 2))));
        assert_eq!(e2.edge_intersection_approx(&e1, tol),
                   EdgeIntersection::Overlap(Edge::new((2, 2), (1, 1))));


        // Overlapping edges. One is fully contained in the other. Same orientations.
        let e1 = Edge::new((0, 0), (4, 4));
        let e2 = Edge::new((1, 1), (2, 2));
        assert!(e1.edges_intersect(&e2).inclusive_bounds());
        assert!(e1.is_collinear(&e2));

        assert_eq!(e1.edge_intersection_approx(&e2, tol),
                   EdgeIntersection::Overlap(e2));
        assert_eq!(e2.edge_intersection_approx(&e1, tol),
                   EdgeIntersection::Overlap(e2));

        // Overlapping edges. One is fully contained in the other. Opposing orientations.
        let e1 = Edge::new((0, 0), (4, 4));
        let e2 = Edge::new((2, 2), (1, 1));
        assert!(e1.edges_intersect(&e2).inclusive_bounds());
        assert!(e1.is_collinear(&e2));

        assert_eq!(e1.edge_intersection_approx(&e2, tol),
                   EdgeIntersection::Overlap(e2.reversed()));
        assert_eq!(e2.edge_intersection_approx(&e1, tol),
                   EdgeIntersection::Overlap(e2));


        // Collinear edge, touch in exactly one point.
        let e1 = Edge::new((0, 0), (1, 1));
        let e2 = Edge::new((1, 1), (2, 2));
        assert!(e1.edges_intersect(&e2).inclusive_bounds());
        assert!(e1.is_collinear(&e2));

        assert_eq!(e1.edge_intersection_approx(&e2, tol),
                   EdgeIntersection::EndPoint(Point::new(1, 1)));
        assert_eq!(e2.edge_intersection_approx(&e1, tol),
                   EdgeIntersection::EndPoint(Point::new(1, 1)));

        // Edges touch in exactly one point. Floating point.
        let e1 = Edge::new((1.1, 1.2), (0., 0.));
        let e2 = Edge::new((1.1, 1.2), (2., 2.));
        assert!(e1.edges_intersect(&e2).inclusive_bounds());
        assert!(e2.edges_intersect(&e1).inclusive_bounds());

        assert_eq!(e1.edge_intersection_approx(&e2, tol),
                   EdgeIntersection::EndPoint(Point::new(1.1, 1.2)));
        assert_eq!(e2.edge_intersection_approx(&e1, tol),
                   EdgeIntersection::EndPoint(Point::new(1.1, 1.2)));

        // Intersection on endpoint with floats.
        let e1 = Edge {
            start: Vector { x: 20.90725737794763, y: 65.33301386746126 },
            end: Vector { x: 32.799556599584776, y: 63.13131890182373 },
        };
        let e2 = Edge {
            start: Vector { x: 22.217533978705163, y: 70.84660296990562 },
            end: Vector { x: 32.799556599584776, y: 63.13131890182373 },
        };
        assert!(e1.edges_intersect(&e2).inclusive_bounds());
        assert!(e2.edges_intersect(&e1).inclusive_bounds());

        assert_eq!(e1.edge_intersection_approx(&e2, tol),
                   EdgeIntersection::EndPoint(Vector { x: 32.799556599584776, y: 63.13131890182373 }));
        assert_eq!(e2.edge_intersection_approx(&e1, tol),
                   EdgeIntersection::EndPoint(Vector { x: 32.799556599584776, y: 63.13131890182373 }));
    }

    #[test]
    fn test_coincident() {
        let e1 = Edge::new((0, 0), (2, 2));
        let e2 = Edge::new((1, 1), (3, 3));
        assert!(e1.is_coincident(&e2));
        assert!(e2.is_coincident(&e1));

        let e1 = Edge::new((0, 0), (3, 3));
        let e2 = Edge::new((1, 1), (2, 2));
        assert!(e1.is_coincident(&e2));
        assert!(e2.is_coincident(&e1));


        let e1 = Edge::new((0, 0), (1, 1));
        let e2 = Edge::new((1, 1), (2, 2));
        assert!(!e1.is_coincident(&e2));
        assert!(!e2.is_coincident(&e1));
    }


    #[test]
    fn test_intersection_of_random_edges() {
        let tol = 1e-12;
        let seed1 = [1u8; 32];
        let seed2 = [2u8; 32];

        let between = Uniform::from(-1.0..1.0);
        let mut rng = StdRng::from_seed(seed1);

        let mut rand_edge = || -> Edge<f64> {
            let points: Vec<(f64, f64)> = (0..2).into_iter()
                .map(|_| (between.sample(&mut rng), between.sample(&mut rng)))
                .collect();

            Edge::new(points[0], points[1])
        };


        let bernoulli_rare = Bernoulli::new(0.05).unwrap();
        let bernoulli_50 = Bernoulli::new(0.5).unwrap();
        let mut rng = StdRng::from_seed(seed2);

        for _i in 0..1000 {

            // Create a random pair of edges. 5% of the edge pairs share an endpoint.
            let (a, b) = {
                let a = rand_edge();

                let b = {
                    let b = rand_edge();

                    if bernoulli_rare.sample(&mut rng) {
                        // Share a point with the other edge.
                        let shared = if bernoulli_50.sample(&mut rng) {
                            a.start
                        } else {
                            a.end
                        };

                        let result = Edge::new(b.start, shared);
                        if bernoulli_50.sample(&mut rng) {
                            result
                        } else {
                            result.reversed()
                        }
                    } else {
                        b
                    }
                };
                (a, b)
            };

            let intersection_ab = a.edge_intersection_approx(&b, tol);
            assert_eq!(intersection_ab != EdgeIntersection::None, a.edges_intersect(&b).inclusive_bounds());

            let intersection_ba = b.edge_intersection_approx(&a, tol);
            assert_eq!(intersection_ba != EdgeIntersection::None, b.edges_intersect(&a).inclusive_bounds());
        }
    }
}