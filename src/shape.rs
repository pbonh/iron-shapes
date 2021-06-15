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

//! Abstractions for geometrical shapes.

use crate::prelude::*;
use crate::traits::{TryBoundingBox, MapPointwise};
use num_traits::NumCast;

/// Abstracted geometrical shape.
#[derive(PartialEq, Eq, Clone, Debug, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum Geometry<T: CoordinateType> {
    /// Point.
    Point(Point<T>),
    /// Edge.
    Edge(Edge<T>),
    /// Rect.
    Rect(Rect<T>),
    /// SimplePolygon.
    SimplePolygon(SimplePolygon<T>),
    /// Polygon.
    Polygon(Polygon<T>),
    /// Path.
    Path(Path<T>),
    /// Text.
    Text(Text<T>),
}

impl<T: CoordinateType> Geometry<T> {
    /// Create a transformed copy of the geometric object.
    pub fn transformed(&self, tf: &SimpleTransform<T>) -> Self {
        let trans = |p| tf.transform_point(p);
        match self {
            Geometry::Point(p) => tf.transform_point(*p).into(),
            Geometry::Edge(g) => g.transform(trans).into(),
            Geometry::Rect(g) => g.transform(trans).into(),
            Geometry::SimplePolygon(g) => g.transform(trans).into(),
            Geometry::Polygon(g) => g.transform(trans).into(),
            Geometry::Path(p) => p.transform(tf).into(),
            Geometry::Text(g) => g.transform(trans).into(),
        }
    }
}

/// Implement `From` for `Geometry`.
macro_rules! geometry_from {
 ( $t:tt ) => {
       impl<T: CoordinateType> From<$t<T>> for Geometry<T> {
            fn from(x: $t<T>) -> Geometry<T> {
                Geometry::$t(x)
            }
        }
 };
}

// Implement `From<_<T>> for Geometry<T>` for all shapes.
geometry_from!(Point);
geometry_from!(Edge);
geometry_from!(Rect);
geometry_from!(SimplePolygon);
geometry_from!(Polygon);
geometry_from!(Path);
geometry_from!(Text);

impl<T: CoordinateType> TryBoundingBox<T> for Geometry<T> {
    /// Calculate the bounding box of this geometrical shape by calling the bounding box method of the concrete type.
    fn try_bounding_box(&self) -> Option<Rect<T>> {
        match self {
            Geometry::Point(p) => p.try_bounding_box(),
            Geometry::Edge(e) => e.try_bounding_box(),
            Geometry::Rect(e) => e.try_bounding_box(),
            Geometry::SimplePolygon(e) => e.try_bounding_box(),
            Geometry::Polygon(e) => e.try_bounding_box(),
            Geometry::Path(p) => p.try_bounding_box(),
            Geometry::Text(t) => t.try_bounding_box()
        }
    }
}


// impl<T: CoordinateType> MapPointwise<T> for Geometry<T> {
//     /// Point wise transformation.
//     fn transform<F>(&self, transformation: F) -> Self
//         where F: Fn(Point<T>) -> Point<T> {
//         match self {
//             Geometry::Point(e) => e.transform(transformation).into(),
//             Geometry::Edge(e) => e.transform(transformation).into(),
//             Geometry::Rect(e) => e.transform(transformation).into(),
//             Geometry::SimplePolygon(e) => e.transform(transformation).into(),
//             Geometry::Polygon(e) => e.transform(transformation).into(),
//             Geometry::Path(p) => unimplemented!(),
//             Geometry::Text(t) => t.transform(transformation).into(),
//         }
//     }
// }

impl<T: CoordinateType + NumCast> DoubledOrientedArea<T> for Geometry<T> {
    /// Area calculation.
    fn area_doubled_oriented(&self) -> T {
        match self {
            Geometry::Point(_) => T::zero(),
            Geometry::Edge(_) => T::zero(),
            Geometry::Rect(e) => e.area_doubled_oriented(),
            Geometry::SimplePolygon(e) => e.area_doubled_oriented(),
            Geometry::Polygon(e) => e.area_doubled_oriented(),
            Geometry::Path(p) => {
                // TODO: Find a way without type conversions.
                T::from(FloatType::round(
                    p.area_approx::<FloatType>() * (2.0 as FloatType)
                )).unwrap()
            }
            Geometry::Text(_) => T::zero(),
        }
    }
}

impl<T: CoordinateType + NumCast> ToPolygon<T> for Geometry<T> {
    /// Convert a geometry into a polygon.
    ///
    /// The coordinate type must implement `NumCast` because there is currently
    /// no way to convert a `Path` into a polygon without converting it to a float type first.
    ///
    /// # Examples
    /// ```
    /// use iron_shapes::prelude::*;
    /// let rect = Rect::new((0, 0), (1, 2));
    /// // Convert the rectangle to a `Geometry`.
    /// let g: Geometry<_> = rect.into();
    /// assert_eq!(g.to_polygon(), rect.to_polygon())
    /// ```
    fn to_polygon(&self) -> Polygon<T> {
        match self {
            Geometry::Point(_) => Polygon::empty(),
            Geometry::Edge(_) => Polygon::empty(),
            Geometry::Rect(e) => e.to_polygon(),
            Geometry::SimplePolygon(e) => Polygon::from(e),
            Geometry::Polygon(e) => e.clone(),
            Geometry::Path(p) => p.to_polygon_approx().cast().into(),
            Geometry::Text(_) => Polygon::empty(),
        }
    }
}

impl<T: CoordinateType + NumCast> Into<Polygon<T>> for Geometry<T> {
    /// Convert a geometry into a polygon.
    fn into(self) -> Polygon<T> {
        self.to_polygon()
    }
}

impl<T: CoordinateType + NumCast, Dst: CoordinateType + NumCast> TryCastCoord<T, Dst> for Geometry<T> {
    type Output = Geometry<Dst>;

    fn try_cast(&self) -> Option<Self::Output> {
        match self {
            Geometry::Point(p) => p.try_cast().map(|s| s.into()),
            Geometry::Edge(e) => e.try_cast().map(|s| s.into()),
            Geometry::Rect(r) => r.try_cast().map(|s| s.into()),
            Geometry::SimplePolygon(p) => p.try_cast().map(|s| s.into()),
            Geometry::Polygon(p) => p.try_cast().map(|s| s.into()),
            Geometry::Path(p) => p.try_cast().map(|s| s.into()),
            Geometry::Text(t) => t.try_cast().map(|s| s.into()),
        }
    }
}

#[test]
/// Sanity check to make sure that area computation is consistent with conversion to polygons.
fn test_convert_to_polygon() {
    let geometries: Vec<Geometry<_>> = vec![
        Point::new(0, 0).into(),
        Edge::new((0, 0), (1, 1)).into(),
        Rect::new((0, 0), (1, 1)).into(),
        SimplePolygon::from(vec![(0, 0), (1, 0), (1, 1)]).into(),
        Polygon::new(vec![(0, 0), (1, 0), (1, 1)]).into()
    ];

    for g in geometries {
        assert_eq!(g.area_doubled_oriented(), g.to_polygon().area_doubled_oriented());
    }
}
