// Copyright (c) 2018-2020 Thomas Kramer.
// SPDX-FileCopyrightText: 2018-2022 Thomas Kramer
//
// SPDX-License-Identifier: AGPL-3.0-or-later

//! Data structures and functions for dealing with rectangles which consist of
//! vertical and horizontal edges.

use crate::cmp::{max, min};
use crate::edge::IntoEdges;
use crate::polygon::{Polygon, ToPolygon};
use crate::prelude::{Point, REdge, Vector};
use crate::traits::*;
use crate::CoordinateType;
use num_traits::{NumCast, One, Zero};
use std::ops::{Add, Div, Mul, Sub};

/// A rectangle which is oriented along the x an y axis and
/// represented by its lower left and upper right corner.
#[derive(Clone, Copy, Hash, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Rect<T> {
    /// Lower left corner of the rectangle.
    pub lower_left: Point<T>,
    /// Upper right corner of the rectangle.
    pub upper_right: Point<T>,
}

impl<T: PartialOrd + Copy> Rect<T> {
    /// Construct the bounding box of the two points. Order does not matter.
    ///
    /// # Examples
    ///
    /// ```
    /// use iron_shapes::prelude::*;
    ///
    /// // Create a rectangle based on two corner points.
    /// let rect1 = Rect::new(Point::new(0, 0), Point::new(1, 2));
    /// // Any type that implements `Into<Point<T>>` can be used for the corner points.
    /// let rect2 = Rect::new((1, 2), (0, 0));
    /// // Ordering of the corner points does not matter.
    /// assert_eq!(rect1, rect2);
    /// // Even though `(0, 0)` was passed as second argument it is recognized as lower left corner.
    /// assert_eq!(rect2.lower_left(), Point::new(0, 0));
    /// ```
    pub fn new<C>(c1: C, c2: C) -> Self
    where
        C: Into<Point<T>>,
    {
        let p1 = c1.into();
        let p2 = c2.into();

        let (x1, x2) = if p1.x < p2.x {
            (p1.x, p2.x)
        } else {
            (p2.x, p1.x)
        };

        let (y1, y2) = if p1.y < p2.y {
            (p1.y, p2.y)
        } else {
            (p2.y, p1.y)
        };

        Rect {
            lower_left: Point::new(x1, y1),
            upper_right: Point::new(x2, y2),
        }
    }
}

impl<T: Copy> Rect<T> {
    /// Get the lower left corner.
    #[inline]
    pub fn lower_left(&self) -> Point<T> {
        self.lower_left
    }

    /// Get the upper left corner.
    #[inline]
    pub fn upper_left(&self) -> Point<T> {
        Point::new(self.lower_left.x, self.upper_right.y)
    }

    /// Get the upper right corner.
    #[inline]
    pub fn upper_right(&self) -> Point<T> {
        self.upper_right
    }

    /// Get the lower right corner.
    #[inline]
    pub fn lower_right(&self) -> Point<T> {
        Point::new(self.upper_right.x, self.lower_left.y)
    }
}

impl<T: PartialOrd + Copy> Rect<T> {
    /// Check if rectangle contains the point.
    /// Inclusive boundaries.
    ///
    /// # Example
    /// ```
    /// use iron_shapes::prelude::*;
    /// let rect = Rect::new((0, 0), (10, 20));
    /// // Contains point somewhere in the center.
    /// assert!(rect.contains_point(Point::new(5, 5)));
    /// // Also contains point on the boundaries.
    /// assert!(rect.contains_point(Point::new(0, 0)));
    /// // Does not contain point outside of the rectangle.
    /// assert!(!rect.contains_point(Point::new(10, 21)));
    /// ```
    pub fn contains_point(&self, p: Point<T>) -> bool {
        self.lower_left.x <= p.x
            && p.x <= self.upper_right.x
            && self.lower_left.y <= p.y
            && p.y <= self.upper_right.y
    }

    /// Check if rectangle contains the point.
    /// Exclusive boundaries.
    ///
    /// # Example
    /// ```
    /// use iron_shapes::prelude::*;
    /// let rect = Rect::new((0, 0), (10, 20));
    /// // Contains point somewhere in the center.
    /// assert!(rect.contains_point_exclusive(Point::new(5, 5)));
    /// // Does not contain points on boundaries.
    /// assert!(!rect.contains_point_exclusive(Point::new(0, 0)));
    /// // Does not contain point outside of the rectangle.
    /// assert!(!rect.contains_point_exclusive(Point::new(10, 21)));
    /// ```
    pub fn contains_point_exclusive(&self, p: Point<T>) -> bool {
        self.lower_left.x < p.x
            && p.x < self.upper_right.x
            && self.lower_left.y < p.y
            && p.y < self.upper_right.y
    }

    /// Check if rectangle contains other rectangle.
    /// Inclusive boundaries.
    ///
    /// # Example
    ///
    /// ```
    /// use iron_shapes::prelude::*;
    ///
    /// let outer = Rect::new((0, 0), (2, 2));
    /// let inner = Rect::new((0, 0), (1, 1));
    /// assert!(outer.contains_rectangle(&inner));
    /// assert!(!inner.contains_rectangle(&outer));
    /// ```
    pub fn contains_rectangle(&self, other: &Self) -> bool {
        self.contains_point(other.lower_left) && self.contains_point(other.upper_right)
    }

    /// Check if rectangle contains other rectangle.
    /// Exclusive boundaries.
    ///
    /// # Example
    ///
    /// ```
    /// use iron_shapes::prelude::*;
    ///
    /// let outer = Rect::new((0, 0), (3, 3));
    /// let inner = Rect::new((1, 1), (2, 2));
    /// assert!(outer.contains_rectangle_exclusive(&inner));
    /// assert!(!inner.contains_rectangle_exclusive(&outer));
    ///
    /// let not_inner = Rect::new((0, 0), (1, 1)); // This shares the boundary with `outer`.
    /// assert!(!outer.contains_rectangle_exclusive(&not_inner));
    /// ```
    pub fn contains_rectangle_exclusive(&self, other: &Self) -> bool {
        self.contains_point_exclusive(other.lower_left)
            && self.contains_point_exclusive(other.upper_right)
    }

    /// Test if the both rectangles touch each other, i.e. if they either share a boundary or are overlapping.
    pub fn touches(&self, other: &Self) -> bool {
        !(self.lower_left.x > other.upper_right.x
            || self.lower_left.y > other.upper_right.y
            || self.upper_right.x < other.lower_left.x
            || self.upper_right.y < other.lower_left.y)
    }

    /// Compute the boolean intersection of two rectangles.
    /// This function excludes the boundaries, hence a zero-area intersection is considered `None`.
    /// See `intersection_inclusive_bounds()` zero-area intersections should be returned as `Some(rectangle)`.
    ///
    /// # Example
    ///
    /// ```
    /// use iron_shapes::prelude::*;
    ///
    /// // Create two overlapping rectangles.
    /// let a = Rect::new((0, 0), (2, 2));
    /// let b = Rect::new((1, 1), (3, 3));
    ///
    /// // Compute the intersection.
    /// assert_eq!(a.intersection(&b), Some(Rect::new((1, 1), (2, 2))));
    ///
    /// // Create a non-overlapping rectangle.
    /// let c = Rect::new((100, 100), (200, 200));
    /// // The intersection with a non-overlapping rectangle is `None`.
    /// assert_eq!(a.intersection(&c), None);
    /// ```
    pub fn intersection(&self, other: &Self) -> Option<Self> {
        let llx = max(self.lower_left.x, other.lower_left.x);
        let lly = max(self.lower_left.y, other.lower_left.y);

        let urx = min(self.upper_right.x, other.upper_right.x);
        let ury = min(self.upper_right.y, other.upper_right.y);

        if llx < urx && lly < ury {
            Some(Rect::new((llx, lly), (urx, ury)))
        } else {
            None
        }
    }

    /// Compute the boolean intersection of two rectangles and include the boundaries.
    /// This allows to get zero-area intersection results for example if the two
    /// rectangles touch on a boundary or one of the rectangle is already zero-area.
    ///
    /// # Example
    ///
    /// ```
    /// use iron_shapes::prelude::*;
    ///
    /// // Create two rectangles which intersect in a single point.
    /// let a = Rect::new((0, 0), (2, 2));
    /// let b = Rect::new((2, 2), (3, 3));
    ///
    /// // Compute the intersection.
    /// assert_eq!(a.intersection_inclusive_bounds(&b), Some(Rect::new((2, 2), (2, 2))));
    ///
    /// ```
    pub fn intersection_inclusive_bounds(&self, other: &Self) -> Option<Self> {
        let llx = max(self.lower_left.x, other.lower_left.x);
        let lly = max(self.lower_left.y, other.lower_left.y);

        let urx = min(self.upper_right.x, other.upper_right.x);
        let ury = min(self.upper_right.y, other.upper_right.y);

        if llx <= urx && lly <= ury {
            Some(Rect::new((llx, lly), (urx, ury)))
        } else {
            None
        }
    }

    /// Create the smallest `Rect` that contains the original `Rect` and the `point`.
    ///
    /// # Example
    ///
    /// ```
    /// use iron_shapes::prelude::*;
    ///
    /// let r1 = Rect::new((0,0), (1,2));
    ///
    /// let r2 = r1.add_point(Point::new(10, 11));
    ///
    /// assert_eq!(r2, Rect::new((0,0), (10,11)));
    ///
    /// ```
    pub fn add_point(&self, point: Point<T>) -> Self {
        Rect::new(
            Point::new(
                min(self.lower_left.x, point.x),
                min(self.lower_left.y, point.y),
            ),
            Point::new(
                max(self.upper_right.x, point.x),
                max(self.upper_right.y, point.y),
            ),
        )
    }

    /// Get the smallest `Rect` that contains both rectangles `self` and `rect`.
    ///
    /// # Example
    ///
    /// ```
    /// use iron_shapes::prelude::*;
    ///
    /// let r1 = Rect::new((0,0), (1,2));
    /// let r2 = Rect::new((4,5), (6,7));
    ///
    /// let r3 = r1.add_rect(&r2);
    ///
    /// assert_eq!(r3, Rect::new((0,0), (6,7)));
    ///
    /// ```
    pub fn add_rect(&self, rect: &Self) -> Self {
        self.add_point(rect.lower_left).add_point(rect.upper_right)
    }
}

impl<T: Sub<Output = T> + Copy + Ord + Zero> Rect<T> {
    /// Compute the shortest from the rectangle to the point `p`.
    /// The distance is zero if the point is inside the rectangle.
    ///
    /// # Example
    ///
    /// ```
    /// use iron_shapes::prelude::*;
    ///
    /// let r = Rect::new((0,0), (10, 10));
    ///
    /// assert_eq!(r.distance_to_point((5, 15).into()), Vector::new(0, 5));
    ///
    /// // Distance to point inside the rectangle is zero.
    /// assert_eq!(r.distance_to_point((5, 5).into()), Vector::new(0, 0));
    ///
    /// ```
    pub fn distance_to_point(&self, p: Point<T>) -> Vector<T> {
        let ll = self.lower_left();
        let ul = self.upper_right();

        // Compute x component of distance.
        let dx_neg = (p.x - ll.x).min(Zero::zero());
        let dx_pos = (p.x - ul.x).max(Zero::zero());
        let dx = dx_neg + dx_pos;

        // Compute y component of distance.
        let dy_neg = (p.y - ll.y).min(Zero::zero());
        let dy_pos = (p.y - ul.y).max(Zero::zero());
        let dy = dy_neg + dy_pos;

        Vector::new(dx, dy)
    }
}

#[test]
fn test_distance_to_point() {
    let rect = Rect::new((10, 10), (20, 20));

    assert_eq!(rect.distance_to_point((15, 15).into()).norm1(), 0);
    assert_eq!(rect.distance_to_point((10, 10).into()).norm1(), 0);
    assert_eq!(rect.distance_to_point((10, 20).into()).norm1(), 0);

    assert_eq!(rect.distance_to_point((0, 0).into()).norm1(), 20);
    assert_eq!(rect.distance_to_point((0, 5).into()).norm1(), 15);
    assert_eq!(rect.distance_to_point((0, 10).into()).norm1(), 10);
    assert_eq!(rect.distance_to_point((0, 15).into()).norm1(), 10);
    assert_eq!(rect.distance_to_point((0, 20).into()).norm1(), 10);
    assert_eq!(rect.distance_to_point((0, 25).into()).norm1(), 15);
}

impl<T: Copy + Sub<Output = T>> Rect<T> {
    /// Compute the width of the rectangle.
    #[inline]
    pub fn width(&self) -> T {
        self.upper_right.x - self.lower_left.x
    }

    /// Compute the height of the rectangle.
    #[inline]
    pub fn height(&self) -> T {
        self.upper_right.y - self.lower_left.y
    }
}

impl<T: Copy + Add<Output = T> + Div<Output = T> + One> Rect<T> {
    /// Get the center point of the rectangle.
    /// When using integer coordinates the resulting
    /// coordinates will be truncated to the next integers.
    pub fn center(&self) -> Point<T> {
        let _2 = T::one() + T::one();
        (self.lower_left() + self.upper_right()) / _2
    }
}

impl<T: Copy + Add<Output = T> + Sub<Output = T> + PartialOrd> Rect<T> {
    /// Create an enlarged copy of this rectangle.
    /// The vertical boundaries will be shifted towards the outside by `add_x`.
    /// The horizontal boundaries will be shifted towards the outside by `add_y`.
    pub fn sized(&self, add_x: T, add_y: T) -> Self {
        Rect::new(
            (self.lower_left.x - add_x, self.lower_left.y - add_y),
            (self.upper_right.x + add_x, self.upper_right.y + add_y),
        )
    }

    /// Create an enlarged copy of this rectangle.
    pub fn sized_isotropic(&self, add: T) -> Self {
        self.sized(add, add)
    }
}

impl<T: Copy + Add<Output = T> + Sub<Output = T> + Mul<Output = T>> DoubledOrientedArea<T>
    for Rect<T>
{
    /// Calculate doubled oriented area of rectangle.
    fn area_doubled_oriented(&self) -> T {
        let diff = self.upper_right - self.lower_left;
        let area = diff.x * diff.y;
        area + area
    }
}

impl<T: Copy> BoundingBox<T> for Rect<T> {
    /// Get bounding box of rectangle (which is equal to the rectangle itself).
    fn bounding_box(&self) -> Rect<T> {
        *self
    }
}

impl<T: Copy> TryBoundingBox<T> for Rect<T> {
    /// Get bounding box of rectangle (always exists).
    fn try_bounding_box(&self) -> Option<Rect<T>> {
        Some(*self)
    }
}

/// Point wise transformation of the two corner points.
impl<T: Copy + PartialOrd> MapPointwise<T> for Rect<T> {
    /// Point wise transformation.
    fn transform<F>(&self, transformation: F) -> Self
    where
        F: Fn(Point<T>) -> Point<T>,
    {
        Self::new(
            transformation(self.lower_left),
            transformation(self.upper_right),
        )
    }
}

/// Iterate over all points of the rectangle.
/// Starts with the lower left corner and iterates counter clock-wise.
impl<'a, T> IntoIterator for &'a Rect<T>
where
    T: Copy,
{
    type Item = Point<T>;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        vec![
            self.lower_left(),
            self.lower_right(),
            self.upper_right(),
            self.upper_left(),
        ]
        .into_iter()
    }
}

/// Iterate over all points of the rectangle.
/// Starts with the lower left corner and iterates counter clock-wise.
impl<T> IntoIterator for Rect<T>
where
    T: Copy,
{
    type Item = Point<T>;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        (&self).into_iter()
    }
}

impl<T: Copy> ToPolygon<T> for Rect<T> {
    fn to_polygon(&self) -> Polygon<T> {
        Polygon::from(self)
    }
}

impl<T: CoordinateType + NumCast, Dst: CoordinateType + NumCast> TryCastCoord<T, Dst> for Rect<T> {
    type Output = Rect<Dst>;

    fn try_cast(&self) -> Option<Self::Output> {
        match (self.lower_left.try_cast(), self.upper_right.try_cast()) {
            (Some(ll), Some(ur)) => Some(Rect::new(ll, ur)),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rect_intersection() {
        let a = Rect::new((0, 0), (2, 4));
        let b = Rect::new((1, 2), (3, 5));
        assert_eq!(a.intersection(&b), Some(Rect::new((1, 2), (2, 4))));

        let a = Rect::new((0, 0), (2, 2));
        let b = Rect::new((1, 1), (3, 3));
        assert_eq!(a.intersection(&b), Some(Rect::new((1, 1), (2, 2))));

        let a = Rect::new((0, 0), (1, 1));
        let b = Rect::new((2, 2), (3, 3));
        assert_eq!(a.intersection(&b), None);

        let a = Rect::new((0, 0), (2, 2));
        let b = Rect::new((1, 2), (5, 5));
        assert_eq!(a.intersection(&b), None);
    }
}

/// Iterator over edges of a rectangle.
#[derive(Clone)]
pub struct RectEdgeIterator<T> {
    rect: Rect<T>,
    pos: u8,
}

impl<T> RectEdgeIterator<T> {
    fn new(rect: Rect<T>) -> Self {
        Self { rect, pos: 0 }
    }
}

impl<T: CoordinateType> Iterator for RectEdgeIterator<T> {
    type Item = REdge<T>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= 4 {
            None
        } else {
            let point = |idx: u8| -> Point<T> {
                match idx {
                    0 => self.rect.lower_right(),
                    1 => self.rect.upper_right(),
                    2 => self.rect.upper_left(),
                    3 => self.rect.lower_left(),
                    _ => unreachable!(),
                }
            };
            let edge = REdge::new(point(self.pos), point((self.pos + 1) % 4));
            self.pos += 1;
            Some(edge)
        }
    }
}

impl<T: CoordinateType> IntoEdges<T> for &Rect<T> {
    type Edge = REdge<T>;
    type EdgeIter = RectEdgeIterator<T>;

    fn into_edges(self) -> Self::EdgeIter {
        RectEdgeIterator::new(*self)
    }
}

#[test]
fn test_edges_iterator() {
    let rect = Rect::new((1, 2), (3, 4));
    let edges: Vec<_> = rect.into_edges().collect();
    assert_eq!(
        edges,
        vec![
            REdge::new((3, 2), (3, 4)),
            REdge::new((3, 4), (1, 4)),
            REdge::new((1, 4), (1, 2)),
            REdge::new((1, 2), (3, 2)),
        ]
    );
}
