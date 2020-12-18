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

//! Multi-polygons are a set of multiple polygons.

use crate::CoordinateType;

use crate::point::Point;
use crate::polygon::Polygon;

pub use crate::traits::{DoubledOrientedArea, BoundingBox, MapPointwise, WindingNumber};

use std::iter::FromIterator;

/// A `MultiPolygon` is a list of polygons. There is no restrictions on the polygons (they can be
/// intersecting, empty, etc.).
#[derive(Clone, Debug, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct MultiPolygon<T>
    where T: CoordinateType {
    /// The list of polygons which defines the content of this multi polygon.
    pub polygons: Vec<Polygon<T>>
}

impl<T: CoordinateType> MultiPolygon<T> {

    /// Create a `MultiPolygon` from a vector of `Polygon`s.
    pub fn new(polygons: Vec<Polygon<T>>) -> Self {
        MultiPolygon {
            polygons
        }
    }

    /// Return the number of polygons.
    pub fn len(&self) -> usize {
        self.polygons.len()
    }
}


impl<T> WindingNumber<T> for MultiPolygon<T>
    where T: CoordinateType {
    fn winding_number(&self, point: Point<T>) -> isize {
        self.polygons.iter()
            .map(|p| p.winding_number(point))
            .sum()
    }
}

impl<T> MapPointwise<T> for MultiPolygon<T>
    where T: CoordinateType {
    fn transform<F: Fn(Point<T>) -> Point<T>>(&self, tf: F) -> Self {
        MultiPolygon::new(
            self.polygons.iter()
                .map(|p| p.transform(&tf))
                .collect()
        )
    }
}

impl<T: CoordinateType, IP: Into<Polygon<T>>> From<IP> for MultiPolygon<T> {
    fn from(x: IP) -> Self {
        MultiPolygon::new(vec![x.into()])
    }
}


impl<T: CoordinateType, IP: Into<Polygon<T>>> FromIterator<IP> for MultiPolygon<T> {
    fn from_iter<I: IntoIterator<Item=IP>>(iter: I) -> Self {
        MultiPolygon::new(
            iter.into_iter()
                .map(|p| p.into()).collect()
        )
    }
}

impl<T: CoordinateType> IntoIterator for MultiPolygon<T> {
    type Item = Polygon<T>;
    type IntoIter = ::std::vec::IntoIter<Polygon<T>>;

    fn into_iter(self) -> Self::IntoIter {
        self.polygons.into_iter()
    }
}

//impl<'a, T: CoordinateType> IntoIterator for &'a MultiPolygon<T> {
//    type Item = &'a Polygon<T>;
//    type IntoIter = ::std::vec::IntoIter<&'a Polygon<T>>;
//
//    fn into_iter(self) -> Self::IntoIter {
//        let p = &self.polygons;
//        p.into_iter()
//    }
//}


