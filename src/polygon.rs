// Copyright (c) 2018-2020 Thomas Kramer.
// SPDX-FileCopyrightText: 2018-2022 Thomas Kramer
//
// SPDX-License-Identifier: AGPL-3.0-or-later

//! This module contains data types and functions for polygons with holes.

use crate::CoordinateType;

use crate::edge::Edge;
use crate::point::Point;
use crate::rect::Rect;

pub use crate::traits::{BoundingBox, DoubledOrientedArea, MapPointwise, WindingNumber};

pub use crate::simple_polygon::*;
use crate::types::*;

use crate::traits::TryCastCoord;
use itertools::Itertools;
use num_traits::{Num, NumCast};
use std::cmp::{Ord, PartialEq};
use std::iter::FromIterator;

/// A polygon possibly with holes. The polygon is defined by a hull and a list of holes
/// which are both `SimplePolygon`s.
#[derive(Clone, Hash, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Polygon<T> {
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
where
    T: CoordinateType,
    Point<T>: From<&'a P>,
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
where
    T: Copy,
    Point<T>: From<P>,
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
where
    T: Copy,
    P: Into<Point<T>>,
{
    fn from_iter<I>(iter: I) -> Self
    where
        I: IntoIterator<Item = P>,
    {
        let exterior: SimplePolygon<T> = SimplePolygon::from_iter(iter);
        Polygon {
            exterior,
            interiors: Vec::new(),
        }
    }
}

/// Create a polygon from a simple polygon.
impl<T> From<SimplePolygon<T>> for Polygon<T> {
    fn from(simple_polygon: SimplePolygon<T>) -> Self {
        Polygon {
            exterior: simple_polygon,
            interiors: Vec::new(),
        }
    }
}

/// Create a polygon from a simple polygon.
impl<T> From<&SimplePolygon<T>> for Polygon<T>
where
    T: Copy,
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
where
    T: Copy,
{
    fn from(rect: Rect<T>) -> Self {
        Polygon::from(vec![
            rect.lower_left(),
            rect.lower_right(),
            rect.upper_right(),
            rect.upper_left(),
        ])
    }
}

/// Create a polygon from a rectangle.
impl<T> From<&Rect<T>> for Polygon<T>
where
    T: Copy,
{
    fn from(rect: &Rect<T>) -> Self {
        Polygon::from(vec![
            rect.lower_left(),
            rect.lower_right(),
            rect.upper_right(),
            rect.upper_left(),
        ])
    }
}

/// Trait for the conversion of a geometric shape to a polygon.
pub trait ToPolygon<T> {
    /// Convert the geometric object into a polygon.
    fn to_polygon(&self) -> Polygon<T>;
}

impl<T> Polygon<T> {
    /// Create empty polygon without any vertices.
    pub fn empty() -> Self {
        Polygon {
            exterior: SimplePolygon::empty(),
            interiors: Vec::new(),
        }
    }
}

impl<T: Copy> Polygon<T> {
    /// Get the number of vertices.
    pub fn len(&self) -> usize {
        self.exterior.len()
    }

    /// Check if polygon has no vertices.
    pub fn is_empty(&self) -> bool {
        self.exterior.is_empty()
    }

    /// Get all exterior edges of the polygon.
    pub fn edges(&self) -> Vec<Edge<T>> {
        self.exterior.edges()
    }

    /// Iterate over all edges of the polygon, including interior edges.
    pub fn all_edges_iter(&self) -> impl Iterator<Item = Edge<T>> + '_ {
        self.exterior.edges_iter().chain(
            self.interiors
                .iter()
                .flat_map(|interior| interior.edges_iter()),
        )
    }
}

impl<T: CoordinateType> Polygon<T> {
    /// Create a new polygon from a sequence of points.
    pub fn new<I>(i: I) -> Self
    where
        I: Into<Self>,
    {
        i.into()
    }

    /// Create a new polygon from a hull and a list of holes.
    pub fn new_with_holes<E, I>(exterior: E, holes: Vec<I>) -> Self
    where
        E: Into<SimplePolygon<T>>,
        I: Into<SimplePolygon<T>>,
    {
        Polygon {
            exterior: exterior.into(),
            interiors: holes.into_iter().map(|i| i.into()).collect(),
        }
    }

    /// Get the convex hull of the polygon.
    ///
    /// Implements Andrew's Monotone Chain algorithm.
    /// See: <http://geomalgorithms.com/a10-_hull-1.html>
    pub fn convex_hull(&self) -> Polygon<T>
    where
        T: Ord,
    {
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
    /// assert_eq!(poly.orientation::<i64>(), Orientation::CounterClockWise);
    ///
    /// ```
    pub fn orientation<Area>(&self) -> Orientation
    where
        Area: Num + From<T> + PartialOrd,
    {
        self.exterior.orientation::<Area>()
    }
}

impl<T> TryBoundingBox<T> for Polygon<T>
where
    T: Copy + PartialOrd,
{
    fn try_bounding_box(&self) -> Option<Rect<T>> {
        // TODO: What if holes exceed the exterior boundary?
        let bbox = self.exterior.try_bounding_box();

        if let Some(bbox) = &bbox {
            debug_assert!(
                self.interiors.iter()
                    .filter_map(|p| p.try_bounding_box())
                    .all(|internal_bbox| bbox.contains_rectangle(&internal_bbox)),
                "Bounding boxes of interior polygons exceed the bounding box of the exterior polygon."
            );
        } else {
            // If the bounding box of the hull is not defined there should also be no
            // defined bounding boxes for the holes.
            let num_internal_bboxes = self
                .interiors
                .iter()
                .filter_map(|p| p.try_bounding_box())
                .count();
            debug_assert_eq!(
                num_internal_bboxes, 0,
                "Polygon with empty zero-vertex hull cannot contain holes."
            );
        }

        bbox
    }
}

impl<T> WindingNumber<T> for Polygon<T>
where
    T: CoordinateType,
{
    /// Calculate the winding number of the polygon around this point.
    ///
    /// TODO: Define how point on edges and vertices is handled.
    ///
    /// See: <http://geomalgorithms.com/a03-_inclusion.html>
    fn winding_number(&self, point: Point<T>) -> isize {
        let ext = self.exterior.winding_number(point);
        let int: isize = self.interiors.iter().map(|p| p.winding_number(point)).sum();
        ext + int
    }
}

impl<T> MapPointwise<T> for Polygon<T>
where
    T: CoordinateType,
{
    fn transform<F: Fn(Point<T>) -> Point<T>>(&self, tf: F) -> Self {
        Polygon {
            exterior: self.exterior.transform(&tf),
            interiors: self.interiors.iter().map(|p| p.transform(&tf)).collect(),
        }
    }
}

impl<T, A> DoubledOrientedArea<A> for Polygon<T>
where
    T: CoordinateType,
    A: Num + From<T>,
{
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
    /// let area: i64 = poly.area_doubled_oriented();
    /// assert_eq!(area, 3);
    ///
    /// ```
    fn area_doubled_oriented(&self) -> A {
        let ext: A = self.exterior.area_doubled_oriented();
        let int = self
            .interiors
            .iter()
            .map(|p| p.area_doubled_oriented())
            .fold(A::zero(), |a, b| a + b);
        ext + int
    }
}

impl<T: PartialEq> Eq for Polygon<T> {}

impl<T> PartialEq for Polygon<T>
where
    T: PartialEq,
{
    /// Equality test for polygons.
    ///
    /// Two polygons are equal iff a cyclic shift on their vertices can be applied
    /// such that the both lists of vertices match exactly.
    fn eq(&self, rhs: &Self) -> bool {
        // TODO: Equality check for polygons with holes.
        assert!(
            self.interiors.is_empty() && rhs.interiors.is_empty(),
            "Equality check for polygons with holes not yet implemented."
        );

        self.exterior == rhs.exterior
    }
}

impl<T: CoordinateType + NumCast, Dst: CoordinateType + NumCast> TryCastCoord<T, Dst>
    for Polygon<T>
{
    type Output = Polygon<Dst>;

    fn try_cast(&self) -> Option<Self::Output> {
        if let Some(new_hull) = self.exterior.try_cast() {
            let new_holes: Vec<_> = self
                .interiors
                .iter()
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

        let doubled_area: i64 = poly.area_doubled_oriented();
        assert_eq!(doubled_area, 2);
    }

    #[test]
    fn test_orientation() {
        use crate::types::Orientation;
        let coords = vec![(0, 0), (1, 0), (1, 1)];

        let poly: Polygon<i32> = (&coords).into();

        assert_eq!(poly.orientation::<i64>(), Orientation::CounterClockWise);
    }

    #[test]
    fn test_bounding_box() {
        use crate::rect::Rect;
        use crate::traits::TryBoundingBox;
        let coords = vec![(1, 0), (-1, -2), (1, 0), (42, 37)];

        let poly: Polygon<i32> = (&coords).into();

        assert_eq!(poly.try_bounding_box(), Some(Rect::new((-1, -2), (42, 37))))
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
