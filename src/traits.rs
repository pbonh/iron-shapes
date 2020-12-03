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
use crate::CoordinateType;

use crate::rect::Rect;
use crate::vector::Vector;
use crate::point::Point;

use num_traits::Float;
use num_traits::cast::NumCast;
use crate::types::Angle;

pub trait BoundingBox<T>
    where T: CoordinateType {
    /// Return the bounding box of this geometry.
    fn bounding_box(&self) -> Rect<T>;
}

pub trait TryBoundingBox<T>
    where T: CoordinateType {
    /// Return the bounding box of this geometry if a bounding box is defined..
    fn try_bounding_box(&self) -> Option<Rect<T>>;
}

impl<T: CoordinateType + NumCast> TryBoundingBox<T> for dyn BoundingBox<T> {
    fn try_bounding_box(&self) -> Option<Rect<T>> {
        Some(self.bounding_box())
    }
}

pub trait DoubledOrientedArea<T>
    where T: CoordinateType {
    /// Calculate doubled oriented area of this geometry.
    /// Using the doubled area allows to compute the area without using fractions.
    fn area_doubled_oriented(&self) -> T;
}

impl<T: CoordinateType + NumCast> dyn DoubledOrientedArea<T> {
    pub fn area<F: Float + NumCast>(&self) -> F {
        F::abs(F::from(self.area_doubled_oriented()).unwrap() / (F::one() + F::one()))
    }
}

pub trait Translate<T>
    where T: CoordinateType {
    fn translate(&self, v: Vector<T>) -> Self;
}

pub trait Scale<T>
    where T: CoordinateType {
    fn scale(&self, factor: T) -> Self;
}

pub trait MapPointwise<T>
    where T: CoordinateType {
    /// Point wise transformation.
    fn transform<F>(&self, transformation: F) -> Self
        where F: Fn(Point<T>) -> Point<T>;
}

impl<S, T> Scale<T> for S
    where T: CoordinateType, S: MapPointwise<T> {
    fn scale(&self, factor: T) -> S {
        self.transform(|p: Point<T>| p * factor)
    }
}

impl<S, T> Translate<T> for S
    where T: CoordinateType, S: MapPointwise<T> {
    fn translate(&self, v: Vector<T>) -> S {
        self.transform(|p: Point<T>| p + v)
    }
}


pub trait WindingNumber<T>
    where T: CoordinateType {
    /// Calculate the winding number of the polygon around this point.
    ///
    /// TODO: Define how point on edges and vertices is handled.
    ///
    /// See: http://geomalgorithms.com/a03-_inclusion.html
    fn winding_number(&self, point: Point<T>) -> isize;

    /// Check if `point` is inside the polygon, i.e. the polygons winds around the point
    /// a non-zero number of times.
    ///
    /// For points on edges the following convention is used:
    /// Points on left or bottom edges are inside, points on right or top edges outside.
    fn contains_point_non_oriented(&self, point: Point<T>) -> bool {
        self.winding_number(point) != 0
    }

    /// Check if `point` is inside the polygon, i.e. the polygon winds around the point
    /// an odd number of times.
    ///
    /// For points on edges the following convention is used:
    /// Points on left or bottom edges are inside, points on right or top edges outside.
    fn contains_point(&self, point: Point<T>) -> bool {
        (self.winding_number(point) % 2 + 2) % 2 == 1
    }
}


/// Rotate by a integer multiple of 90 degrees.
pub trait RotateOrtho<T: CoordinateType> {
    fn rotate_ortho(&self, a: Angle) -> Self;
}

impl<S, T> RotateOrtho<T> for S
    where T: CoordinateType, S: MapPointwise<T> {
    fn rotate_ortho(&self, a: Angle) -> S {
        self.transform(|p: Point<T>|
            match a {
                Angle::R0 => p,
                Angle::R90 => Point::new(
                    T::zero() - p.y,
                    p.x,
                ),
                Angle::R180 => -p,
                Angle::R270 => Point::new(
                    p.y,
                    T::zero() - p.x,
                )
            }
        )
    }
}


/// Mirror at the x or y axis.
pub trait Mirror<T: CoordinateType> {
    fn mirror_x(&self) -> Self;
    fn mirror_y(&self) -> Self;
}

impl<S, T> Mirror<T> for S
    where T: CoordinateType, S: MapPointwise<T> {
    fn mirror_x(&self) -> S {
        self.transform(|p| Point::new(T::zero() - p.x, p.y))
    }

    fn mirror_y(&self) -> S {
        self.transform(|p| Point::new(p.x, T::zero() - p.y))
    }
}

//
//pub trait EachPoint<'a, T: 'a>
//    where T: CoordinateType {
//    type EachPoints: Iterator<Item=&'a Point<T>>;
//    fn each_point(&self) -> Self::EachPoints;
//}
