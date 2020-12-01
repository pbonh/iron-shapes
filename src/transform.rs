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

use crate::vector::Vector;
use crate::point::Point;
use crate::matrix2d::*;
use crate::traits::{Angle, RotateOrtho, Transform, Mirror, Translate, Scale};

use std::ops::Mul;
use crate::types::FloatType;
use crate::matrix3d::Matrix3d;
use num_traits::{Zero, Float};

#[derive(Clone, Hash, PartialEq, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Matrix2dTransform<T: CoordinateType> {
    matrix: Matrix2d<T>
}

impl<T: CoordinateType> Matrix2dTransform<T> {
    pub fn new(matrix: Matrix2d<T>) -> Self {
        Matrix2dTransform { matrix }
    }

    /// Create a rotation by an integer multiple of 90 degrees.
    pub fn new_rotation(angle: Angle) -> Self {
        let zero = T::zero();
        let one = T::one();
        let minus_one = zero - one;

        let matrix = match angle {
            Angle::R0 => Matrix2d::new(one, zero, zero, one),
            Angle::R90 => Matrix2d::new(zero, minus_one, one, zero),
            Angle::R180 => Matrix2d::new(minus_one, zero, zero, minus_one),
            Angle::R270 => Matrix2d::new(zero, one, minus_one, zero),
        };

        Matrix2dTransform::new(matrix)
    }

    /// Create a scaling by a factor.
    pub fn new_scaling(factor: T) -> Self {
        let zero = T::zero();
        Matrix2dTransform::new(
            Matrix2d::new(factor, zero, zero, factor)
        )
    }

    /// Mirror at the x-axis.
    pub fn new_mirror_x() -> Self {
        let zero = T::zero();
        let one = T::one();
        let minus_one = zero - one;
        Matrix2dTransform::new(
            Matrix2d::new(minus_one, zero, zero, one)
        )
    }

    /// Mirror at the y-axis.
    pub fn new_mirror_y() -> Self {
        let zero = T::zero();
        let one = T::one();
        let minus_one = zero - one;
        Matrix2dTransform::new(
            Matrix2d::new(one, zero, zero, minus_one)
        )
    }

    // pub fn is_unitary(&self) -> bool {
    //     // TODO
    //     true
    // }

    /// Apply the transformation to a single point.
    pub fn transform_point(&self, p: Point<T>) -> Point<T> {
        self.matrix.mul_vector(p)
    }

    /// Return the matrix describing this transformation.
    pub fn to_matrix2d(&self) -> Matrix2d<T> {
        self.matrix.clone()
    }

    /// Get the inverse transformation.
    pub fn inverted(&self) -> Self {
        unimplemented!()
    }
}

#[test]
fn test_matrix_transform_rotations() {
    let p = Point::new(1, 0);

    assert_eq!(Matrix2dTransform::new_rotation(Angle::R0).transform_point(p), p);
    assert_eq!(Matrix2dTransform::new_rotation(Angle::R90).transform_point(p), Point::new(0, 1));
    assert_eq!(Matrix2dTransform::new_rotation(Angle::R180).transform_point(p), Point::new(-1, 0));
    assert_eq!(Matrix2dTransform::new_rotation(Angle::R270).transform_point(p), Point::new(0, -1));
}


#[derive(Clone, Hash, PartialEq, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Rot90Transform {
    angle: Angle
}

impl Rot90Transform {
    pub fn new(angle: Angle) -> Self {
        Rot90Transform { angle }
    }
    pub fn is_unitary(&self) -> bool {
        true
    }
    /// Apply the transformation to a single point.
    pub fn transform_point<T: CoordinateType>(&self, p: Point<T>) -> Point<T> {
        p.rotate_ortho(self.angle)
    }

    pub fn magnification<T: CoordinateType>(&self) -> T {
        T::one()
    }
    pub fn try_magnification<T: CoordinateType>(&self) -> Option<T> {
        Some(self.magnification())
    }
}

// pub trait Translate<T>
//     where T: CoordinateType {
//     fn translate(&self, v: Vector<T>) -> Self;
// }
//
// pub trait Scale<T>
//     where T: CoordinateType {
//     fn scale(&self, factor: T) -> Self;
// }
//
// pub trait Transform<T>
//     where T: CoordinateType {
//     /// Point wise transformation.
//     fn transform<F>(&self, transformation: F) -> Self
//         where F: Fn(Point<T>) -> Point<T>;
// }
//
// impl<S, T> Scale<T> for S
//     where T: CoordinateType, S: Transform<T> {
//     fn scale(&self, factor: T) -> S {
//         self.transform(|p: Point<T>| p * factor)
//     }
// }
//
// impl<S, T> Translate<T> for S
//     where T: CoordinateType, S: Transform<T> {
//     fn translate(&self, v: Vector<T>) -> S {
//         self.transform(|p: Point<T>| p + v)
//     }
// }

/// Describes a geometric transformation that consists of a optional mirroring along the x-axis
/// followed by a rotation by a multiple of 90 degrees
/// followed by a displacement.
#[derive(Clone, Default, PartialEq, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SimpleTransform<T: CoordinateType> {
    pub mirror: bool,
    pub rotation: Angle,
    pub magnification: T,
    pub displacement: Vector<T>,
}

impl<T: CoordinateType> SimpleTransform<T> {
    pub fn new(mirror: bool, rotation: Angle, magnification: T, displacement: Vector<T>) -> Self {
        SimpleTransform {
            mirror,
            rotation,
            magnification,
            displacement,
        }
    }

    pub fn transform_point(&self, p: Point<T>) -> Point<T> {
        if self.mirror {
            p.mirror_x()
        } else {
            p
        }
            .rotate_ortho(self.rotation)
            .scale(self.magnification)
            .translate(self.displacement)
    }
}


#[derive(Clone, PartialEq, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ComplexTransform<T: CoordinateType> {
    mirror: bool,
    rotation: FloatType,
    displacement: Vector<T>,
}

/// Affine transformation represented as a 3x3 matrix like:
/// ```txt
/// m11 m12 0
/// m21 m22 0
/// m31 m32 1
/// ```
#[derive(Clone, Hash, PartialEq, Eq, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Matrix3dTransform<T: CoordinateType> {
    pub m11: T,
    pub m12: T,
    pub m21: T,
    pub m22: T,
    pub m31: T,
    pub m32: T,
}

impl<T: CoordinateType> Matrix3dTransform<T> {
    pub fn new(m11: T, m12: T, m21: T, m22: T, m31: T, m32: T) -> Self {
        Matrix3dTransform {
            m11,
            m12,
            m21,
            m22,
            m31,
            m32,
        }
    }

    /// Get the identity transform.
    pub fn identity() -> Self {
        Self::translate(Vector::zero())
    }

    /// Create a translation by a vector.
    pub fn translate(t: Vector<T>) -> Self {
        Self::new(
            T::one(), T::zero(),
            T::zero(), T::one(),
            t.x, t.y,
        )
    }

    /// Create a rotation by an integer multiple of 90 degrees.
    pub fn rotate90(angle: Angle) -> Self {
        let zero = T::zero();
        let one = T::one();
        let minus_one = zero - one;

        let (a, b, c, d) = match angle {
            Angle::R0 => (one, zero, zero, one),
            Angle::R90 => (zero, one, minus_one, zero),
            Angle::R180 => (minus_one, zero, zero, minus_one),
            Angle::R270 => (zero, minus_one, one, zero),
        };

        Self::new(
            a, b,
            c, d,
            T::zero(), T::zero(),
        )
    }

    /// Create a scaling by a factor.
    pub fn scale(factor: T) -> Self {
        let zero = T::zero();
        Self::new(
            factor, zero,
            zero, factor,
            zero, zero,
        )
    }

    /// Mirror at the x-axis.
    pub fn mirror_x() -> Self {
        let zero = T::zero();
        let one = T::one();
        let minus_one = zero - one;
        Self::new(
            minus_one, zero,
            zero, zero,
            zero, zero,
        )
    }

    /// Mirror at the y-axis.
    pub fn mirror_y() -> Self {
        let zero = T::zero();
        let one = T::one();
        let minus_one = zero - one;
        Self::new(
            one, zero,
            zero, minus_one,
            zero, zero,
        )
    }

    /// Apply the transformation to a single point.
    pub fn transform_point(&self, p: Point<T>) -> Point<T> {
        Point::new(
            p.x * self.m11 + p.y * self.m21 + self.m31,
            p.x * self.m12 + p.y * self.m22 + self.m32,
        )
    }


    /// Return the matrix describing this transformation.
    pub fn to_matrix3d(&self) -> Matrix3d<T> {
        Matrix3d::new(
            self.m11, self.m12, T::zero(),
            self.m21, self.m22, T::zero(),
            self.m31, self.m32, T::zero(),
        )
    }

    /// Get the inverse transformation if it exists.
    pub fn try_invert(&self) -> Option<Self> {
        unimplemented!()
    }
}

#[test]
fn test_identity() {
    let p = Point::new(1, 2);
    let tf = Matrix3dTransform::identity();
    assert_eq!(tf.transform_point(p), p);
}

#[test]
fn test_translate() {
    let p = Point::new(1, 2);
    let tf = Matrix3dTransform::translate(Vector::new(10, 100));
    assert_eq!(tf.transform_point(p), Point::new(11, 102));
}

#[test]
fn test_rotate90() {
    let p = Point::new(1, 2);
    let tf = Matrix3dTransform::rotate90(Angle::R0);
    assert_eq!(tf.transform_point(p), Point::new(1, 2));
    let tf = Matrix3dTransform::rotate90(Angle::R90);
    assert_eq!(tf.transform_point(p), Point::new(-2, 1));
    let tf = Matrix3dTransform::rotate90(Angle::R180);
    assert_eq!(tf.transform_point(p), Point::new(-1, -2));
    let tf = Matrix3dTransform::rotate90(Angle::R270);
    assert_eq!(tf.transform_point(p), Point::new(2, -1));
}

#[test]
fn test_scale() {
    let p = Point::new(1, 2);
    let tf = Matrix3dTransform::scale(2);
    assert_eq!(tf.transform_point(p), Point::new(2, 4));
}

impl<T: CoordinateType + Float> Matrix3dTransform<T> {
    /// Create a rotation by an arbitrary angle (in radians).
    pub fn rotation(phi: T) -> Self {
        let zero = T::zero();
        let cos = phi.cos();
        let sin = phi.sin();
        Self::new(
            cos, sin,
            zero - sin, cos,
            T::zero(), T::zero(),
        )
    }
}

#[test]
fn test_rotate() {
    let p = Point::new(1.0, 0.0);
    let pi = std::f64::consts::PI;
    let tf = Matrix3dTransform::rotation(pi);
    assert!(
        (tf.transform_point(p) - Point::new(-1.0, 0.0)).norm2_squared() < 1e-6
    );
    let tf = Matrix3dTransform::rotation(pi*0.5);
    assert!(
        (tf.transform_point(p) - Point::new(0.0, 1.0)).norm2_squared() < 1e-6
    );
}