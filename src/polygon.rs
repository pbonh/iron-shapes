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

//! This module contains data types and functions for polygons with holes.

use crate::CoordinateType;

use crate::point::Point;
use crate::edge::Edge;
use crate::rect::Rect;

pub use crate::traits::{DoubledOrientedArea, BoundingBox, MapPointwise, WindingNumber};

use crate::types::*;
pub use crate::simple_polygon::*;

use std::iter::FromIterator;
use std::cmp::{Ord, PartialEq};
use crate::traits::TryCastCoord;
use num_traits::NumCast;
use itertools::Itertools;

/// A polygon possibly with holes. The polygon is defined by a hull and a list of holes
/// which are both `SimplePolygon`s.
#[derive(Clone, Hash, Debug, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Polygon<T>
    where T: CoordinateType {
    /// The outer hull of the polygon.
    pub exterior: SimplePolygon<T>,
    /// A list of holes in the polygon.
    pub interiors: Vec<SimplePolygon<T>>,
}

/// Shorthand notation for creating a polygon.
///
/// # Example
/// ```
/// # #[macro_use]
/// # extern crate iron_shapes;
/// # fn main() {
/// use iron_shapes::prelude::*;
/// let p = polygon!((0, 0), (1, 0), (1, 1));
/// assert_eq!(p, Polygon::new(vec![(0, 0), (1, 0), (1, 1)]));
/// # }
/// ```
#[macro_export]
macro_rules! polygon {
 ($($x:expr),*) => {Polygon::new((vec![$($x),*]))}
}

/// Create a polygon from a `Vec` of values convertible to `Point`s.
impl<'a, T, P> From<&'a Vec<P>> for Polygon<T>
    where T: CoordinateType,
          Point<T>: From<&'a P>
{
    fn from(vec: &'a Vec<P>) -> Self {
        Polygon {
            exterior: vec.into(),
            interiors: Vec::new(),
        }
    }
}


/// Create a polygon from a `Vec` of values convertible to `Point`s.
impl<T, P> From<Vec<P>> for Polygon<T>
    where T: CoordinateType,
          Point<T>: From<P>
{
    fn from(vec: Vec<P>) -> Self {
        Polygon {
            exterior: vec.into(),
            interiors: Vec::new(),
        }
    }
}


/// Create a polygon from a iterator of values convertible to `Point`s.
impl<T, P> FromIterator<P> for Polygon<T>
    where T: CoordinateType,
          P: Into<Point<T>>
{
    fn from_iter<I>(iter: I) -> Self
        where I: IntoIterator<Item=P>
    {
        let exterior: SimplePolygon<T> = SimplePolygon::from_iter(iter);
        Polygon {
            exterior,
            interiors: Vec::new(),
        }
    }
}

/// Create a polygon from a simple polygon.
impl<T> From<SimplePolygon<T>> for Polygon<T>
    where T: CoordinateType
{
    fn from(simple_polygon: SimplePolygon<T>) -> Self {
        Polygon {
            exterior: simple_polygon,
            interiors: Vec::new(),
        }
    }
}

/// Create a polygon from a simple polygon.
impl<T> From<&SimplePolygon<T>> for Polygon<T>
    where T: CoordinateType
{
    fn from(simple_polygon: &SimplePolygon<T>) -> Self {
        Polygon {
            exterior: simple_polygon.clone(),
            interiors: Vec::new(),
        }
    }
}

/// Create a polygon from a rectangle.
impl<T> From<Rect<T>> for Polygon<T>
    where T: CoordinateType
{
    fn from(rect: Rect<T>) -> Self {
        Polygon::from(
            vec![rect.lower_left(), rect.lower_right(),
                 rect.upper_right(), rect.upper_left()]
        )
    }
}

/// Create a polygon from a rectangle.
impl<T> From<&Rect<T>> for Polygon<T>
    where T: CoordinateType
{
    fn from(rect: &Rect<T>) -> Self {
        Polygon::from(
            vec![rect.lower_left(), rect.lower_right(),
                 rect.upper_right(), rect.upper_left()]
        )
    }
}

/// Trait for the conversion of a geometric shape to a polygon.
pub trait ToPolygon<T>
    where T: CoordinateType {
    /// Convert the geometric object into a polygon.
    fn to_polygon(&self) -> Polygon<T>;
}


impl<T: CoordinateType> Polygon<T> {
    /// Create a new polygon from a sequence of points.
    pub fn new<I>(i: I) -> Self
        where I: Into<Self> {
        i.into()
    }

    /// Create empty polygon without any vertices.
    pub fn empty() -> Self {
        Polygon {
            exterior: SimplePolygon::empty(),
            interiors: Vec::new(),
        }
    }

    /// Create a new polygon from a hull and a list of holes.
    pub fn new_with_holes<E, I>(exterior: E, holes: Vec<I>) -> Self
        where E: Into<SimplePolygon<T>>,
              I: Into<SimplePolygon<T>> {
        Polygon {
            exterior: exterior.into(),
            interiors: holes.into_iter().map(|i| i.into()).collect(),
        }
    }

    /// Get the number of vertices.
    pub fn len(&self) -> usize {
        self.exterior.len()
    }

    /// Get all exterior edges of the polygon.
    pub fn edges(&self) -> Vec<Edge<T>> {
        self.exterior.edges()
    }

    /// Get the convex hull of the polygon.
    ///
    /// Implements Andrew's Monotone Chain algorithm.
    /// See: http://geomalgorithms.com/a10-_hull-1.html
    pub fn convex_hull(&self) -> Polygon<T>
        where T: Ord {
        self.exterior.convex_hull().into()
    }

    /// Get the vertex with lowest x-coordinate of the exterior polygon.
    /// Prefer lower y-coordinates to break ties.
    ///
    /// # Examples
    ///
    /// ```
    /// use iron_shapes::polygon::Polygon;
    /// use iron_shapes::point::Point;
    /// let coords = vec![(0, 0), (1, 0), (-1, 2), (-1, 1)];
    ///
    /// let poly = Polygon::new(coords);
    ///
    /// assert_eq!(poly.lower_left_vertex(), Point::new(-1, 1));
    ///
    /// ```
    pub fn lower_left_vertex(&self) -> Point<T> {
        self.exterior.lower_left_vertex()
    }


    /// Get the orientation of the exterior polygon.
    ///
    /// # Examples
    ///
    /// ```
    /// use iron_shapes::polygon::Polygon;
    /// use iron_shapes::point::Point;
    /// use iron_shapes::types::Orientation;
    /// let coords = vec![(0, 0), (3, 0), (3, 1)];
    ///
    /// let poly = Polygon::new(coords);
    ///
    /// assert_eq!(poly.orientation(), Orientation::CounterClockWise);
    ///
    /// ```
    pub fn orientation(&self) -> Orientation {
        self.exterior.orientation()
    }
}

impl<T> BoundingBox<T> for Polygon<T>
    where T: CoordinateType {
    fn bounding_box(&self) -> Rect<T> {
        // TODO: What if holes exceed the exterior boundary?
        let bbox = self.exterior.bounding_box();

        debug_assert!(
            self.interiors.iter()
                .all(|p| bbox.contains_rectangle(p.bounding_box())),
            "Bounding boxes of interior polygons exceed the bounding box of the exterior polygon."
        );

        bbox
    }
}

impl<T> WindingNumber<T> for Polygon<T>
    where T: CoordinateType {
    /// Calculate the winding number of the polygon around this point.
    ///
    /// TODO: Define how point on edges and vertices is handled.
    ///
    /// See: http://geomalgorithms.com/a03-_inclusion.html
    fn winding_number(&self, point: Point<T>) -> isize {
        let ext = self.exterior.winding_number(point);
        let int: isize = self.interiors.iter()
            .map(|p| p.winding_number(point))
            .sum();
        ext + int
    }
}

impl<T> MapPointwise<T> for Polygon<T>
    where T: CoordinateType {
    fn transform<F: Fn(Point<T>) -> Point<T>>(&self, tf: F) -> Self {
        Polygon {
            exterior: self.exterior.transform(&tf),
            interiors: self.interiors.iter()
                .map(|p| p.transform(&tf))
                .collect(),
        }
    }
}

impl<T: CoordinateType> DoubledOrientedArea<T> for Polygon<T> {
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
    /// use iron_shapes::polygon::{Polygon, DoubledOrientedArea};
    /// let coords = vec![(0, 0), (3, 0), (3, 1)];
    ///
    /// let poly = Polygon::new(coords);
    ///
    /// assert_eq!(poly.area_doubled_oriented(), 3);
    ///
    /// ```
    fn area_doubled_oriented(&self) -> T {
        let ext = self.exterior.area_doubled_oriented();
        let int = self.interiors.iter()
            .map(|p| p.area_doubled_oriented())
            .fold(T::zero(), |a, b| a + b);
        ext + int
    }
}

impl<T> PartialEq for Polygon<T>
    where T: CoordinateType {
    /// Equality test for polygons.
    ///
    /// Two polygons are equal iff a cyclic shift on their vertices can be applied
    /// such that the both lists of vertices match exactly.
    fn eq(&self, rhs: &Self) -> bool {

        // TODO: Equality check for polygons with holes.
        assert!(self.interiors.is_empty() && rhs.interiors.is_empty(),
                "Equality check for polygons with holes not yet implemented.");

        self.exterior == rhs.exterior
    }
}

impl<T: CoordinateType + NumCast, Dst: CoordinateType + NumCast> TryCastCoord<T, Dst> for Polygon<T> {
    type Output = Polygon<Dst>;

    fn try_cast(&self) -> Option<Self::Output> {

        if let Some(new_hull) = self.exterior.try_cast() {
            let new_holes : Vec<_> = self.interiors.iter()
                .map(|hole| hole.try_cast())
                .while_some()
                .collect();
            if new_holes.len() == self.interiors.len() {
                Some(Polygon::new_with_holes(new_hull, new_holes))
            } else {
                // Some wholes could not be casted.
                None
            }
        } else {
            // The hull could not be casted.
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::prelude::*;

    #[test]
    fn test_create_polygon() {
        let coords = vec![(0, 0), (1, 0), (1, 1), (0, 1)];

        let poly: Polygon<i32> = (&coords).into();

        assert_eq!(poly.exterior.len(), coords.len());
    }

    #[test]
    fn test_area() {
        let coords = vec![(0, 0), (1, 0), (1, 1), (0, 1)];

        let poly: Polygon<i32> = (&coords).into();

        assert_eq!(poly.area_doubled_oriented(), 2);
    }

    #[test]
    fn test_orientation() {
        use crate::types::Orientation;
        let coords = vec![(0, 0), (1, 0), (1, 1)];

        let poly: Polygon<i32> = (&coords).into();

        assert_eq!(poly.orientation(), Orientation::CounterClockWise);
    }

    #[test]
    fn test_bounding_box() {
        use crate::traits::BoundingBox;
        use crate::rect::Rect;
        let coords = vec![(1, 0), (-1, -2), (1, 0), (42, 37)];

        let poly: Polygon<i32> = (&coords).into();

        assert_eq!(poly.bounding_box(), Rect::new((-1, -2), (42, 37)))
    }

    #[test]
    fn test_winding_number() {
        let coords = vec![(0, 0), (2, 0), (2, 2), (0, 2)];

        let poly: Polygon<i32> = (&coords).into();

        assert_eq!(poly.winding_number(Point::new(1, 1)), 1);
        assert_eq!(poly.winding_number(Point::new(-1, -1)), 0);
        assert_eq!(poly.winding_number(Point::new(10, 10)), 0);

        // Test point on edges
        assert_eq!(poly.winding_number(Point::new(1, 0)), 1); // Bottom edge
        assert_eq!(poly.winding_number(Point::new(2, 1)), 0); // Right edge
        assert_eq!(poly.winding_number(Point::new(1, 2)), 0); // Top edge
        assert_eq!(poly.winding_number(Point::new(0, 1)), 1); // Left edge

        // Test point on vertex.
        assert_eq!(poly.winding_number(Point::new(0, 0)), 1);
        assert_eq!(poly.winding_number(Point::new(2, 0)), 0);
        assert_eq!(poly.winding_number(Point::new(2, 2)), 0);
        assert_eq!(poly.winding_number(Point::new(0, 2)), 0);
    }

    #[test]
    fn test_convex_hull() {
        let poly = Polygon::new(vec![(0, 0), (1, 1), (2, 0), (2, 2), (0, 2)]);
        let exp_hull = Polygon::new(vec![(0, 0), (2, 0), (2, 2), (0, 2)]);
        assert_eq!(poly.convex_hull(), exp_hull);

        let poly = Polygon::new(vec![(1, 0), (2, 1), (1, 2), (1, 1), (0, 1)]);
        let exp_hull = Polygon::new(vec![(1, 0), (2, 1), (1, 2), (0, 1)]);
        assert_eq!(poly.convex_hull(), exp_hull);

        // Degenerate case. All x coordinates are the same.
        let poly = Polygon::new(vec![(0, 0), (0, 1), (0, 7)]);
        let exp_hull = Polygon::new(vec![(0, 0), (0, 7)]);
        assert_eq!(poly.convex_hull(), exp_hull);

        // Degenerate case. All y coordinates are the same.
        let poly = Polygon::new(vec![(0, 0), (1, 0), (7, 0)]);
        let exp_hull = Polygon::new(vec![(0, 0), (7, 0)]);
        assert_eq!(poly.convex_hull(), exp_hull);


        // Degenerate case. All points are equal.
        let poly4 = Polygon::new(vec![(0, 0), (0, 0), (0, 0)]);
        let exp_hull4 = Polygon::new(vec![(0, 0)]);
        assert_eq!(poly4.convex_hull(), exp_hull4);
    }


    #[test]
    fn test_partial_eq() {
        let poly1 = Polygon::new(vec![(0, 0), (1, 0), (1, 1)]);
        let poly2 = Polygon::new(vec![(1, 1), (0, 0), (1, 0)]);
        let poly3 = Polygon::new(vec![(0, 0), (1, 0), (1, 2)]);

        assert_eq!(poly1, poly2);
        assert_eq!(poly2, poly1);

        assert_ne!(poly1, poly3)
    }
}
