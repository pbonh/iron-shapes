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

//! Data structures and functions for 3x3 matrices.

use crate::CoordinateType;

/// 3x3 matrix of the form.
/// ```txt
/// [[ m11, m12, m13 ],
///  [ m21, m22, m23 ],
///  [ m31, m32, m33 ]]
/// ```
#[derive(Clone, Hash, PartialEq, Eq, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Matrix3d<T: CoordinateType> {
    /// m11
    pub m11: T,
    /// m12
    pub m12: T,
    /// m13
    pub m13: T,
    /// m21
    pub m21: T,
    /// m22
    pub m22: T,
    /// m23
    pub m23: T,
    /// m31
    pub m31: T,
    /// m32
    pub m32: T,
    /// m33
    pub m33: T,
}

impl<T> Matrix3d<T>
    where T: CoordinateType
{
    /// Create a new 3x3 matrix with entries of the form:
    /// ```txt
    /// [[ m11, m12, m13 ],
    ///  [ m21, m22, m23 ],
    ///  [ m31, m32, m33 ]]
    /// ```
    pub fn new(m11: T, m12: T, m13: T, m21: T, m22: T, m23: T, m31: T, m32: T, m33: T) -> Self {
        Matrix3d {
            m11,
            m12,
            m13,
            m21,
            m22,
            m23,
            m31,
            m32,
            m33,
        }
    }

    /// Return the identity matrix.
    pub fn identity() -> Self {
        let _0 = T::zero();
        let _1 = T::one();
        Self::new(_1, _0, _0,
                  _0, _1, _0,
                  _0, _0, _1)
    }

    /// Compute product of the matrix with a scalar.
    pub fn mul_scalar(&self, rhs: T) -> Self {
        Matrix3d::new(
            self.m11 * rhs, self.m12 * rhs, self.m13 * rhs,
            self.m21 * rhs, self.m22 * rhs, self.m23 * rhs,
            self.m31 * rhs, self.m32 * rhs, self.m33 * rhs,
        )
    }

    /// Element-wise addition of two matrices.
    pub fn add(&self, rhs: &Self) -> Self {
        Matrix3d::new(
            self.m11 + rhs.m11, self.m12 + rhs.m12, self.m13 + rhs.m13,
            self.m21 + rhs.m21, self.m22 + rhs.m22, self.m23 + rhs.m23,
            self.m31 + rhs.m31, self.m32 + rhs.m32, self.m33 + rhs.m33,
        )
    }

    /// Compute multiplication with a column vector `A*rhs`.
    pub fn mul_column_vector(&self, rhs: &(T, T, T)) -> (T, T, T) {
        (
            self.m11 * rhs.0 + self.m12 * rhs.1 + self.m13 * rhs.2,
            self.m21 * rhs.0 + self.m22 * rhs.1 + self.m23 * rhs.2,
            self.m31 * rhs.0 + self.m32 * rhs.1 + self.m33 * rhs.2,
        )
    }

    /// Matrix-matrix multiplication.
    pub fn mul_matrix(&self, rhs: &Self) -> Self {
        let a = self;
        let b = rhs;
        let c11 = a.m11 * b.m11 + a.m12 * b.m21 + a.m13 * b.m31;
        let c12 = a.m11 * b.m12 + a.m12 * b.m22 + a.m13 * b.m32;
        let c13 = a.m11 * b.m13 + a.m12 * b.m23 + a.m13 * b.m33;
        let c21 = a.m21 * b.m11 + a.m22 * b.m21 + a.m23 * b.m31;
        let c22 = a.m21 * b.m12 + a.m22 * b.m22 + a.m23 * b.m32;
        let c23 = a.m21 * b.m13 + a.m22 * b.m23 + a.m23 * b.m33;
        let c31 = a.m31 * b.m11 + a.m32 * b.m21 + a.m33 * b.m31;
        let c32 = a.m31 * b.m12 + a.m32 * b.m22 + a.m33 * b.m32;
        let c33 = a.m31 * b.m13 + a.m32 * b.m23 + a.m33 * b.m33;
        Self::new(
            c11, c12, c13,
            c21, c22, c23,
            c31, c32, c33,
        )
    }

    /// Compute the transpose of the matrix.
    pub fn transpose(&self) -> Self {
        Self::new(
            self.m11, self.m21, self.m31,
            self.m12, self.m22, self.m32,
            self.m13, self.m23, self.m33,
        )
    }

    /// Test if this matrix is the identity matrix.
    pub fn is_identity(&self) -> bool {
        self == &Self::identity()
    }

    /// Test if this matrix is unitary.
    pub fn is_unitary(&self) -> bool {
        self.mul_matrix(&self.transpose()).is_identity()
    }

    /// Compute the determinant of this matrix.
    pub fn determinant(&self) -> T {
        /*
        ```python
        import sympy as sp
        m = sp.Matrix([[sp.Symbol(f"a.m{i}{j}") for j in range(1,4)] for i in range(1,4)])
        m.det()
        ```
         */
        let a = self;
        a.m11 * a.m22 * a.m33 - a.m11 * a.m23 * a.m32 - a.m12 * a.m21 * a.m33
            + a.m12 * a.m23 * a.m31 + a.m13 * a.m21 * a.m32 - a.m13 * a.m22 * a.m31
    }

    /// Compute the inverse matrix.
    /// `None` will be returned if the determinant is zero and the matrix
    /// is not invertible.
    pub fn try_inverse(&self) -> Option<Self> {
        /*
        ```python
        import sympy as sp
        m = sp.Matrix([[sp.Symbol(f"a.m{i}{j}") for j in range(1,4)] for i in range(1,4)])
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
                (a.m22 * a.m33 - a.m23 * a.m32) / det, (a.m13 * a.m32 - a.m12 * a.m33) / det, (a.m12 * a.m23 - a.m13 * a.m22) / det,
                (a.m23 * a.m31 - a.m21 * a.m33) / det, (a.m11 * a.m33 - a.m13 * a.m31) / det, (a.m13 * a.m21 - a.m11 * a.m23) / det,
                (a.m21 * a.m32 - a.m22 * a.m31) / det, (a.m12 * a.m31 - a.m11 * a.m32) / det, (a.m11 * a.m22 - a.m12 * a.m21) / det,
            ))
        } else {
            None
        }
    }
}

impl<T: CoordinateType> Default for Matrix3d<T> {
    fn default() -> Self {
        Self::identity()
    }
}

#[test]
fn test_matrix_multiplication() {
    let id: Matrix3d<i32> = Matrix3d::identity();
    assert_eq!(id.mul_matrix(&id), id);
}

#[test]
fn test_mul_column_vector() {
    let id: Matrix3d<i32> = Matrix3d::identity();
    let v = (1, 2, 3);
    assert_eq!(id.mul_column_vector(&v), v);
}

#[test]
fn test_inverse() {
    let m = Matrix3d::new(1.0, 4.0, 2.0,
                          1.0, 2.0, 4.0,
                          2.0, 0.0, 4.0);
    let i = m.try_inverse().unwrap();
    assert_eq!(m.mul_matrix(&i), Matrix3d::identity());
    assert_eq!(i.mul_matrix(&m), Matrix3d::identity());
}