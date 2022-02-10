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

//! Transforms are used to describe the location, rotation, scaling and mirroring
//! of geometric shapes.

use crate::CoordinateType;

use crate::vector::Vector;
use crate::point::Point;
use crate::matrix2d::*;
use crate::traits::{RotateOrtho, Mirror, Translate, Scale, TryCastCoord};
use crate::types::Angle;

use std::ops::Mul;
use crate::types::FloatType;
use crate::matrix3d::Matrix3d;
use num_traits::{Zero, Float, NumCast};

/// General geometric transformation.
pub trait Transformation {
    /// Coordinate type of input points.
    type SourceCoord;
    /// Coordinate type of output points.
    type DestinationCoord;

    /// Apply the transform to a point.
    fn transform_point(&self, p: Point<Self::SourceCoord>) -> Point<Self::DestinationCoord>;
}

/// Geometric transformation which preserves parallelism.
/// Adds 'shear' to the [`SimilarityTransform`].
pub trait AffineTransform: Transformation {
    /// Get the both eigen-vectors which characterize the 'shear'.
    /// Eigen-vectors are the vectors whose direction is not changed when applying the transformation.
    fn eigen_vectors(&self) -> (Vector<Self::SourceCoord>, Vector<Self::SourceCoord>);

    // fn eigen_values(&self) -> (); // TODO
}

/// Geometric transformation which preserves angles and ratios of distances.
/// Adds resizing to the [`IsometricTransform`].
pub trait SimilarityTransform: AffineTransform {
    /// Transform a distance.
    fn transform_distance(&self, distance: Self::SourceCoord) -> Self::DestinationCoord;
}

/// Geometric transformation which preserves angles and ratios of distances.
/// Adds resizing by integer numbers to the [`IsometricRTransform`].
pub trait SimilarityRTransform: SimilarityTransform<SourceCoord=Self::Coord, DestinationCoord=Self::Coord> {
    /// Type or source and destination coordinates.
    type Coord;
    // TODO
}

/// Geometric transformation which preserves angles and distances (e.g. euclidean transform).
pub trait IsometricTransform: SimilarityTransform {
    // TODO
}

/// Geometric transformation which preserves angles and distances (e.g. euclidean transform) but
/// allows only rotations by a multiple of 90 degrees.
pub trait IsometricRTransform: IsometricTransform + SimilarityRTransform {
    // TODO

    /// Return the rotation by a multiple of 90 degrees.
    fn rotation(&self) -> Angle;
}

/// Geometric transformation which preserves oriented angles and distances (i.e. translation).
pub trait DisplacementTransform: IsometricRTransform {
    /// Get the displacement vector.
    fn displacement(&self) -> Vector<Self::SourceCoord>;
}

// // Blanket implementation.
// impl<Tf: DisplacementTransform> IsometricRTransform for Tf {
//     fn rotation(&self) -> Angle {
//         Angle::R0
//     }
// }
//
// /// A geometric displacement transformation.
// pub struct Displacement<T> {
//     disp: Vector<T>
// }
//
// impl<T: Copy> DisplacementTransform<Coord=T> for Displacement<T> {
//     fn displacement(&self) -> Vector<Self::SourceCoord> {
//         *self.disp
//     }
// }


/// Description of a transformation in the euclidean plane by a 2x2 matrix `A`.
/// Transforming a point `p` is computed by the matrix product `A*p`.
#[derive(Clone, Hash, PartialEq, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Matrix2dTransform<T: CoordinateType> {
    matrix: Matrix2d<T>
}

impl<T: CoordinateType> Matrix2dTransform<T> {
    /// Create a new transform based on a matrix.
    pub fn new(matrix: Matrix2d<T>) -> Self {
        Matrix2dTransform { matrix }
    }

    /// Create a rotation by an integer multiple of 90 degrees.
    pub fn new_rotation90(angle: Angle) -> Self {
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
        self.matrix.mul_column_vector(p.into()).into()
    }

    /// Return the matrix describing this transformation.
    pub fn to_matrix2d(&self) -> Matrix2d<T> {
        self.matrix.clone()
    }

    /// Get the inverse transformation.
    pub fn try_invert(&self) -> Option<Self> {
        self.matrix.try_inverse()
            .map(|inv| Self { matrix: inv })
    }
}

#[test]
fn test_matrix_transform_rotations() {
    let p = Point::new(1, 0);

    assert_eq!(Matrix2dTransform::new_rotation90(Angle::R0).transform_point(p), p);
    assert_eq!(Matrix2dTransform::new_rotation90(Angle::R90).transform_point(p), Point::new(0, 1));
    assert_eq!(Matrix2dTransform::new_rotation90(Angle::R180).transform_point(p), Point::new(-1, 0));
    assert_eq!(Matrix2dTransform::new_rotation90(Angle::R270).transform_point(p), Point::new(0, -1));
}

/// Transformation that consists only of a rotation by a multiple of 90 degrees
/// around the origin `(0, 0)`.
#[derive(Clone, Hash, PartialEq, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Rot90Transform {
    angle: Angle
}

impl Rot90Transform {
    /// Create a new rotation transformation.
    pub fn new(angle: Angle) -> Self {
        Rot90Transform { angle }
    }
    /// This transformation is always unitary. Returns always `true`.
    pub fn is_unitary(&self) -> bool {
        true
    }
    /// Apply the transformation to a single point.
    pub fn transform_point<T: CoordinateType>(&self, p: Point<T>) -> Point<T> {
        p.rotate_ortho(self.angle)
    }

    /// Get the magnification of this transformation. Always `1`.
    pub fn magnification<T: CoordinateType>(&self) -> T {
        T::one()
    }

    /// Get the magnification of this transformation. Always `Some(1)`.
    pub fn try_magnification<T: CoordinateType>(&self) -> Option<T> {
        Some(self.magnification())
    }
}

/// Describes a geometric transformation that consists of a optional mirroring along the x-axis
/// followed by a rotation by a multiple of 90 degrees
/// followed by a displacement.
#[derive(Copy, Clone, PartialEq, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SimpleTransform<T> {
    /// Mirror on the x-axis?
    pub mirror: bool,
    /// Rotation by a multiple of 90 degrees.
    pub rotation: Angle,
    /// Enlargement.
    pub magnification: T,
    /// Translation.
    pub displacement: Vector<T>,
}

impl<T: CoordinateType> Default for SimpleTransform<T> {
    fn default() -> Self {
        SimpleTransform::identity()
    }
}

impl<T: CoordinateType> SimpleTransform<T> {
    /// Create a new transformation.
    pub fn new(mirror: bool, rotation: Angle, magnification: T, displacement: Vector<T>) -> Self {
        SimpleTransform {
            mirror,
            rotation,
            magnification,
            displacement,
        }
    }

    /// Get the identity transform.
    pub fn identity() -> Self {
        Self::translate(Vector::zero())
    }

    /// Create a translation by a vector.
    pub fn translate<V: Into<Vector<T>>>(v: V) -> Self {
        let t = v.into();
        Self::new(
            false, Angle::R0,
            T::one(), t,
        )
    }

    /// Create a rotation by an integer multiple of 90 degrees.
    /// Rotation center is (0, 0).
    pub fn rotate90(angle: Angle) -> Self {
        Self::new(
            false, angle,
            T::one(), Vector::zero(),
        )
    }

    /// Create a rotation arount `rotation_center` by an integer multiple of 90 degrees.
    pub fn rotate90_around(angle: Angle, rotation_center: Point<T>) -> Self {
        Self::translate(Point::zero() - rotation_center)
            .then(&Self::rotate90(angle))
            .then(&Self::translate(rotation_center))
    }

    /// Rotate 90 degrees counter-clock wise.
    pub fn rotate_ccw90() -> Self {
        Self::rotate90(Angle::R90)
    }

    /// Rotate 90 degrees counter-clock wise.
    pub fn rotate_cw90() -> Self {
        Self::rotate90(Angle::R270)
    }

    /// Create a scaling by a factor.
    pub fn scale(factor: T) -> Self {
        Self::new(
            false, Angle::R0,
            factor, Vector::zero(),
        )
    }

    /// Create a transformation that mirrors at the x-axis.
    pub fn mirror_x() -> Self {
        Self::new(
            true, Angle::R0,
            T::one(), Vector::zero(),
        )
    }

    /// Create a transformation that mirrors at the y-axis.
    pub fn mirror_y() -> Self {
        Self::new(
            true, Angle::R180,
            T::one(), Vector::zero(),
        )
    }

    /// Transform a distance.
    pub fn transform_distance(&self, d: T) -> T {
        d * self.magnification
    }

    /// Transform a single point.
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

    /// Convert to a matrix transformation.
    pub fn to_matrix_transform(&self) -> Matrix3dTransform<T> {
        if self.mirror {
            Matrix3dTransform::mirror_x()
        } else {
            Matrix3dTransform::identity()
        }
            .then_rotate90(self.rotation)
            .then_scale(self.magnification)
            .then_translate(self.displacement)
    }

    /// Return a new transformation that is equal to applying
    /// first `self` then `t`.
    pub fn then(&self, t: &Self) -> Self {
        let d = t.transform_point(self.displacement.into());
        let r = if t.mirror {
            -self.rotation
        } else {
            self.rotation
        };
        Self {
            mirror: self.mirror ^ t.mirror,
            rotation: r + t.rotation,
            magnification: self.magnification * t.magnification,
            displacement: d.v(),
        }
    }
}

#[test]
fn test_simple_transform_combine() {
    let t1 = SimpleTransform::new(false, Angle::R90,
                                  1, (1, 2).into());
    let t2 = SimpleTransform::new(true, Angle::R90,
                                  1, (3, 4).into());

    let p = Point::new(10, 11);
    assert_eq!(t2.transform_point(t1.transform_point(p)), t1.then(&t2).transform_point(p));
    assert_eq!(t1.transform_point(t2.transform_point(p)), t2.then(&t1).transform_point(p));
}

/// Transformation described by a mirroring at the `x` axis,
/// then a rotation around the origin, then a scaling, then a translation.
/// This transformation allows rotations by arbitrary angles.
#[derive(Clone, PartialEq, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ComplexTransform<T: CoordinateType> {
    /// Mirror at `x`-axis?
    mirror: bool,
    /// Rotation by an arbitrary angle (radians).
    rotation: FloatType,
    /// Scaling.
    magnification: FloatType,
    /// Translation.
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
    /// m11
    pub m11: T,
    /// m12
    pub m12: T,
    /// m21
    pub m21: T,
    /// m22
    pub m22: T,
    /// m31. Used to express the `x` component of the translation.
    pub m31: T,
    /// m32. Used to express the `y` component of the translation.
    pub m32: T,
}

impl<T: CoordinateType> Matrix3dTransform<T> {
    /// Create a new transform based on the matrix elements.
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
    pub fn translate<V: Into<Vector<T>>>(v: V) -> Self {
        let t = v.into();
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
            zero, one,
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


    /// Return the 3x3 matrix describing this transformation.
    pub fn to_matrix3d(&self) -> Matrix3d<T> {
        Matrix3d::new(
            self.m11, self.m12, T::zero(),
            self.m21, self.m22, T::zero(),
            self.m31, self.m32, T::zero(),
        )
    }

    /// Return the 2x2 matrix that describes this transformation
    /// without any translation.
    pub fn to_matrix2d(&self) -> Matrix2d<T> {
        Matrix2d::new(
            self.m11, self.m12,
            self.m21, self.m22,
        )
    }

    /// Compute the determinant of the 3x3 matrix that describes this transformation.
    pub fn determinant(&self) -> T {
        /*
         ```python
         import sympy as sp
         m = sp.Matrix([[sp.Symbol(f"self.m{i}{j}") for j in range(1,4)] for i in range(1,4)])
         m[0,2] = 0
         m[1,2] = 0
         m[2,2] = 1
         print(m.det())
         ```
          */
        self.m11 * self.m22 - self.m12 * self.m21
    }

    /// Get the inverse transformation if it exists.
    pub fn try_invert(&self) -> Option<Self> {
        /*
         ```python
         import sympy as sp
         m = sp.Matrix([[sp.Symbol(f"a.m{i}{j}") for j in range(1,4)] for i in range(1,4)])
         m[0,2] = 0
         m[1,2] = 0
         m[2,2] = 1
         det = sp.Symbol('det')

         # Compute inverse multiplied with determinant.
         inv_det = m.inv() * m.det()
         inv = inv_det / det
         # Print inverse with factored-out determinant.
         print(repr(inv).replace('[', '').replace(']', ''))
         ```
          */
        let a = self;

        // Compute determinant.
        let det = a.determinant();
        if !det.is_zero() {
            Some(Self::new(
                a.m22 / det, T::zero() - a.m12 / det,
                T::zero() - a.m21 / det, a.m11 / det,
                (a.m21 * a.m32 - a.m22 * a.m31) / det, (a.m12 * a.m31 - a.m11 * a.m32) / det,
            ))
        } else {
            None
        }
    }

    /// Return a new transformation that is equal to applying
    /// first `self` then `t`.
    pub fn then(&self, t: &Self) -> Self {
        Self::new(
            self.m11 * t.m11 + self.m12 * t.m21,
            self.m11 * t.m12 + self.m12 * t.m22,
            self.m21 * t.m11 + self.m22 * t.m21,
            self.m21 * t.m12 + self.m22 * t.m22,
            self.m31 * t.m11 + self.m32 * t.m21 + t.m31,
            self.m31 * t.m12 + self.m32 * t.m22 + t.m32,
        )
    }

    /// Create a new transformation with an additional scaling.
    pub fn then_scale(&self, factor: T) -> Self {
        self.then(&Self::scale(factor))
    }

    /// Create a new transformation with an additional translation.
    pub fn then_translate<V: Into<Vector<T>>>(&self, v: V) -> Self {
        self.then(&Self::translate(v))
    }

    /// Create a new transformation with an additional rotation by a multiple of 90 degrees.
    pub fn then_rotate90(&self, angle: Angle) -> Self {
        self.then(&Self::rotate90(angle))
    }

    /// Create a new transformation with an additional mirroring at the x-axis.
    pub fn then_mirror_x(&self) -> Self {
        self.then(&Self::mirror_x())
    }

    /// Create a new transformation with an additional mirroring at the y-axis.
    pub fn then_mirror_y(&self) -> Self {
        self.then(&Self::mirror_y())
    }

    /// Get the translation part of this affine transformation.
    pub fn get_translation(&self) -> Vector<T> {
        Vector::new(self.m31, self.m32)
    }
}

impl<T: CoordinateType> Mul for Matrix3dTransform<T> {
    type Output = Matrix3dTransform<T>;

    /// Shortcut for `self.then(&rhs)`.
    fn mul(self, rhs: Self) -> Self::Output {
        self.then(&rhs)
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
    assert_eq!(tf.get_translation(), Vector::new(10, 100));
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

    /// Create a new transformation with an additional rotation.
    pub fn then_rotate(&self, phi: T) -> Self {
        self.then(&Self::rotation(phi))
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
    let tf = Matrix3dTransform::rotation(pi * 0.5);
    assert!(
        (tf.transform_point(p) - Point::new(0.0, 1.0)).norm2_squared() < 1e-6
    );
}

#[test]
fn test_then() {
    let tf1 = Matrix3dTransform::translate((1, 2));
    let id: Matrix3dTransform<i32> = Matrix3dTransform::identity();
    assert_eq!(tf1.then(&id), tf1);

    let tf1_tf1 = Matrix3dTransform::translate((2, 4));
    assert_eq!(tf1.then(&tf1), tf1_tf1);

    let tf3 = Matrix3dTransform::rotate90(Angle::R90);
    assert_eq!(tf3.then(&tf3).then(&tf3).then(&tf3), id);
}

#[test]
fn test_invert() {
    let id: Matrix3dTransform<i32> = Matrix3dTransform::identity();
    assert_eq!(id.try_invert(), Some(id));

    let tf1 = Matrix3dTransform::translate((1, 2));
    let tf1_inv = tf1.try_invert().unwrap();
    assert_eq!(tf1.then(&tf1_inv), Matrix3dTransform::identity());

    let tf2 = Matrix3dTransform::translate((1, 2))
        .then_rotate90(Angle::R90)
        .then_mirror_x();
    assert!(tf2.to_matrix2d().is_unitary());

    let tf2_inv = tf2.try_invert().unwrap();
    assert_eq!(tf2.then(&tf2_inv), Matrix3dTransform::identity());

    assert_eq!(tf2.try_invert().unwrap().try_invert(), Some(tf2));
}

#[test]
fn test_invert_float() {
    let tf = Matrix3dTransform::rotation(1.234)
        .then_scale(12345.6)
        .then_translate((1.2, 34.5));
    let tf_inv = tf.try_invert().unwrap();
    let p = Point::new(42.42, -1.0);
    let p2 = tf_inv.transform_point(tf.transform_point(p));
    assert!((p - p2).norm2_squared() < 1e-6); // Test for approximate equality.
}

impl<T: CoordinateType + NumCast, Dst: CoordinateType + NumCast> TryCastCoord<T, Dst> for SimpleTransform<T> {
    type Output = SimpleTransform<Dst>;

    fn try_cast(&self) -> Option<Self::Output> {
        match (self.displacement.try_cast(), Dst::from(self.magnification)) {
            (Some(displacement), Some(magnification)) => Some(
                Self::Output {
                    mirror: self.mirror,
                    displacement,
                    magnification,
                    rotation: self.rotation,
                }),
            _ => None
        }
    }
}
