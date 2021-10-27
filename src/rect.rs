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

//! Data structures and functions for dealing with rectangles which consist of
//! vertical and horizontal edges.

use crate::point::Point;
use crate::traits::*;
use crate::cmp::{min, max};
use crate::CoordinateType;
use num_traits::NumCast;
use crate::polygon::{ToPolygon, Polygon};

/// A rectangle which is oriented along the x an y axis and
/// represented by its lower left and upper right corner.
#[derive(Clone, Copy, Hash, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Rect<T> {
    /// Lower left corner of the rectangle.
    pub lower_left: Point<T>,
    /// Upper right corner of the rectangle.
    pub upper_right: Point<T>,
}

impl<T: PartialEq> Eq for Rect<T> {}

impl<T: PartialEq> PartialEq for Rect<T> {
    fn eq(&self, other: &Self) -> bool {
        self.lower_left == other.lower_left && self.upper_right == other.upper_right
    }
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
        where C: Into<Point<T>> {
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

impl<T: CoordinateType> Rect<T> {
    /// Get the center point of the rectangle.
    /// When using integer coordinates the resulting
    /// coordinates will be truncated to the next integers.
    pub fn center(&self) -> Point<T> {
        let _2 = T::one() + T::one();
        (self.lower_left() + self.upper_right()) / _2
    }

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
        self.lower_left.x <= p.x && p.x <= self.upper_right.x &&
            self.lower_left.y <= p.y && p.y <= self.upper_right.y
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
        self.lower_left.x < p.x && p.x < self.upper_right.x &&
            self.lower_left.y < p.y && p.y < self.upper_right.y
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
        self.contains_point_exclusive(other.lower_left) && self.contains_point_exclusive(other.upper_right)
    }

    /// Test if the both rectangles touch each other, i.e. if they either share a boundary or are overlapping.
    pub fn touches(&self, other: &Self) -> bool {
        !(
            self.lower_left.x > other.upper_right.x ||
                self.lower_left.y > other.upper_right.y ||
                self.upper_right.x < other.lower_left.x ||
                self.upper_right.y < other.lower_left.y
        )
    }

    /// Compute the boolean intersection of two rectangles.
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
            Point::new(min(self.lower_left.x, point.x),
                       min(self.lower_left.y, point.y)),
            Point::new(max(self.upper_right.x, point.x),
                       max(self.upper_right.y, point.y)),
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
        self.add_point(rect.lower_left)
            .add_point(rect.upper_right)
    }

    /// Create an enlarged copy of this rectangle.
    /// The vertical boundaries will be shifted towards the outside by `add_x`.
    /// The horizontal boundaries will be shifted towards the outside by `add_y`.
    pub fn sized(&self, add_x: T, add_y: T) -> Self {
        Rect::new((self.lower_left.x - add_x, self.lower_left.y - add_y),
                  (self.upper_right.x + add_x, self.upper_right.y + add_y))
    }
}

impl<T: CoordinateType> DoubledOrientedArea<T> for Rect<T> {
    /// Calculate doubled oriented area of rectangle.
    fn area_doubled_oriented(&self) -> T {
        let diff = self.upper_right - self.lower_left;
        diff.x * diff.y * (T::one() + T::one())
    }
}

impl<T: CoordinateType> BoundingBox<T> for Rect<T> {
    /// Get bounding box of rectangle (which is equal to the rectangle itself).
    fn bounding_box(&self) -> Rect<T> {
        self.clone()
    }
}

impl<T: CoordinateType> TryBoundingBox<T> for Rect<T> {
    /// Get bounding box of rectangle (always exists).
    fn try_bounding_box(&self) -> Option<Rect<T>> {
        Some(self.clone())
    }
}

/// Point wise transformation of the two corner points.
impl<T> MapPointwise<T> for Rect<T>
    where T: CoordinateType
{
    /// Point wise transformation.
    fn transform<F>(&self, transformation: F) -> Self
        where F: Fn(Point<T>) -> Point<T> {
        Self::new(
            transformation(self.lower_left),
            transformation(self.upper_right),
        )
    }
}

/// Iterate over all points of the rectangle.
/// Starts with the lower left corner and iterates counter clock-wise.
impl<'a, T> IntoIterator for &'a Rect<T>
    where T: CoordinateType {
    type Item = Point<T>;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        vec![self.lower_left(), self.lower_right(),
             self.upper_right(), self.upper_left()].into_iter()
    }
}

/// Iterate over all points of the rectangle.
/// Starts with the lower left corner and iterates counter clock-wise.
impl<T> IntoIterator for Rect<T>
    where T: CoordinateType {
    type Item = Point<T>;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        (&self).into_iter()
    }
}

impl<T: CoordinateType> ToPolygon<T> for Rect<T> {
    fn to_polygon(&self) -> Polygon<T> {
        Polygon::from(self)
    }
}

impl<T: CoordinateType + NumCast, Dst: CoordinateType + NumCast> TryCastCoord<T, Dst> for Rect<T> {
    type Output = Rect<Dst>;

    fn try_cast(&self) -> Option<Self::Output> {
        match (self.lower_left.try_cast(), self.upper_right.try_cast()) {
            (Some(ll), Some(ur)) => Some(Rect::new(ll, ur)),
            _ => None
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