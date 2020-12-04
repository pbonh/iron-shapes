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
use num_traits::NumCast;

/// Abstracted geometrical shape.
#[derive(PartialEq, Eq, Clone, Debug, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum Geometry<T: CoordinateType> {
    Point(Point<T>),
    Edge(Edge<T>),
    Rect(Rect<T>),
    SimplePolygon(SimplePolygon<T>),
    Polygon(Polygon<T>),
    Path(Path<T>),
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

impl<T: CoordinateType + NumCast> BoundingBox<T> for Geometry<T> {
    /// Calculate the bounding box of this geometrical shape by calling the bounding box method of the concrete type.
    fn bounding_box(&self) -> Rect<T> {
        match self {
            Geometry::Point(p) => Rect::new(p, p),
            Geometry::Edge(e) => Rect::new(e.start, e.end),
            Geometry::Rect(e) => e.bounding_box(),
            Geometry::SimplePolygon(e) => e.bounding_box(),
            Geometry::Polygon(e) => e.bounding_box(),
            Geometry::Path(p) => p.bounding_box()
        }
    }
}

impl<T: CoordinateType> MapPointwise<T> for Geometry<T> {
    /// Point wise transformation.
    fn transform<F>(&self, transformation: F) -> Self
        where F: Fn(Point<T>) -> Point<T> {
        match self {
            Geometry::Point(e) => e.transform(transformation).into(),
            Geometry::Edge(e) => e.transform(transformation).into(),
            Geometry::Rect(e) => e.transform(transformation).into(),
            Geometry::SimplePolygon(e) => e.transform(transformation).into(),
            Geometry::Polygon(e) => e.transform(transformation).into(),
            Geometry::Path(p) => unimplemented!()
        }
    }
}

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
        }
    }
}

impl<T: CoordinateType> ToPolygon<T> for Geometry<T> {
    /// Convert a geometry into a polygon.
    ///
    /// # Examples
    /// ```
    /// use iron_shapes::prelude::*;
    /// ```
    fn to_polygon(&self) -> Polygon<T> {
        match self {
            Geometry::Point(_) => Polygon::empty(),
            Geometry::Edge(_) => Polygon::empty(),
            Geometry::Rect(e) => e.into(),
            Geometry::SimplePolygon(e) => Polygon::from(e),
            Geometry::Polygon(e) => e.clone(),
            Geometry::Path(p) => unimplemented!()
        }
    }
}

impl<T: CoordinateType> Into<Polygon<T>> for Geometry<T> {
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
        SimplePolygon::new(vec![(0, 0), (1, 0), (1, 1)]).into(),
        Polygon::new(vec![(0, 0), (1, 0), (1, 1)]).into()
    ];

    for g in geometries {
        assert_eq!(g.area_doubled_oriented(), g.to_polygon().area_doubled_oriented());
    }
}
