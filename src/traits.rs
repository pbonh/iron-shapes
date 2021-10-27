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

//! Common traits for geometrical objects.

use crate::rect::Rect;
use crate::vector::Vector;
use crate::point::Point;

use num_traits::{Float, Zero};
use num_traits::cast::NumCast;
use crate::types::Angle;
use std::ops::{Sub, Add, Mul};

/// Calculation of the 'bounding box', i.e. the smallest rectangle that contains the geometrical object.
pub trait BoundingBox<T>: TryBoundingBox<T> {
    /// Return the bounding box of this geometry.
    fn bounding_box(&self) -> Rect<T>;
}

/// Try the calculation of the 'bounding box', i.e. the smallest rectangle that contains the geometrical object.
/// In some cases this is not always possible, so the try might fail.
/// For instance a set of polygons does not have a bounding box if the set is empty.
pub trait TryBoundingBox<T> {
    /// Return the bounding box of this geometry if a bounding box is defined.
    fn try_bounding_box(&self) -> Option<Rect<T>>;
}


/// Calculate the doubled oriented area of a geometry.
/// Using the doubled area allows to compute the area without using fractions. This is especially
/// helpful when computing in integer coordinates.
pub trait DoubledOrientedArea<T> {
    /// Calculate doubled oriented area of this geometry.
    /// Using the doubled area allows to compute the area without using fractions.
    fn area_doubled_oriented(&self) -> T;
}

impl<T: NumCast> dyn DoubledOrientedArea<T> {
    /// Compute the area of a geometrical shape by first computing
    /// doubled oriented area and then taking absolute value of the half.
    pub fn area<F: Float + NumCast>(&self) -> F {
        F::abs(F::from(self.area_doubled_oriented()).unwrap() / (F::one() + F::one()))
    }
}

/// Translate the geometrical object by a vector.
pub trait Translate<T> {
    /// Translate the geometrical object by a vector `v`.
    fn translate(&self, v: Vector<T>) -> Self;
}

/// Scale the geometrical shape. Scaling center is the origin `(0, 0)`.
pub trait Scale<T> {
    /// Scale the geometrical shape. Scaling center is the origin `(0, 0)`.
    fn scale(&self, factor: T) -> Self;
}

/// Transform the geometrical object by transforming each point of it.
pub trait MapPointwise<T> {
    /// Point wise transformation.
    fn transform<F>(&self, transformation: F) -> Self
        where F: Fn(Point<T>) -> Point<T>;
}

impl<S, T> Scale<T> for S
    where T: Copy + Mul<Output=T>, S: MapPointwise<T> {
    fn scale(&self, factor: T) -> S {
        self.transform(|p: Point<T>| p * factor)
    }
}

impl<S, T> Translate<T> for S
    where T: Copy + Add<Output=T>, S: MapPointwise<T> {
    fn translate(&self, v: Vector<T>) -> S {
        self.transform(|p: Point<T>| p + v)
    }
}

/// Compute the winding number of a geometrical object around a point.
/// The winding number is used to check if a point is contained in a shape.
pub trait WindingNumber<T> {
    /// Calculate the winding number of the polygon around this point.
    ///
    /// TODO: Define how point on edges and vertices is handled.
    ///
    /// See: <http://geomalgorithms.com/a03-_inclusion.html>
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
pub trait RotateOrtho<T> {
    /// Rotate the geometrical shape by a multiple of 90 degrees.
    fn rotate_ortho(&self, a: Angle) -> Self;
}

impl<S, T> RotateOrtho<T> for S
    where T: Copy + Zero + Sub<Output=T>, S: MapPointwise<T> {
    fn rotate_ortho(&self, a: Angle) -> S {
        self.transform(|p: Point<T>|
            match a {
                Angle::R0 => p,
                Angle::R90 => Point::new(
                    T::zero() - p.y,
                    p.x,
                ),
                Angle::R180 => Point::zero()-p.v(),
                Angle::R270 => Point::new(
                    p.y,
                    T::zero() - p.x,
                )
            }
        )
    }
}


/// Mirror at the x or y axis.
pub trait Mirror<T> {
    /// Mirror this shape at the `x`-axis.
    fn mirror_x(&self) -> Self;
    /// Mirror this shape at the `y`-axis.
    fn mirror_y(&self) -> Self;
}

impl<S, T> Mirror<T> for S
    where T: Copy + Zero + Sub<Output=T>, S: MapPointwise<T> {
    /// Return the geometrical object mirrored at the `x` axis.
    fn mirror_x(&self) -> S {
        self.transform(|p| Point::new(T::zero() - p.x, p.y))
    }

    /// Return the geometrical object mirrored at the `y` axis.
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

/// This trait defines the type-casting of the coordinate types for geometrical objects.
pub trait TryCastCoord<Src, Dst>
    where Src: NumCast,
          Dst: NumCast {
    /// Output type of the cast. This is likely the same geometrical type just with other
    /// coordinate types.
    type Output;

    /// Try to cast to target data type.
    ///
    /// Conversion from float to int can fail and will return `None`.
    /// Float values like infinity or non-a-number
    /// have no integer representation.
    fn try_cast(&self) -> Option<Self::Output>;

    /// Cast to target data type.
    ///
    /// Conversion from float to int can fail and panic because float values
    /// like infinity or non-a-number have no integer representation.
    ///
    /// # Panics
    /// Panics if casting of the coordinate values fails. For instance when trying to cast a 'NaN' float into a integer.
    /// Use `try_cast` to detect and handle this failure.
    fn cast(&self) -> Self::Output {
        self.try_cast().unwrap()
    }
}
