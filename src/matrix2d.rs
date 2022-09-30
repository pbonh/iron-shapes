// Copyright (c) 2018-2020 Thomas Kramer.
// SPDX-FileCopyrightText: 2018-2022 Thomas Kramer
//
// SPDX-License-Identifier: AGPL-3.0-or-later

//! Data structures and functions for 2x2 matrices.

use crate::CoordinateType;

use crate::vector::Vector;

/// A 2x2 matrix of the form:
/// ```txt
/// [[ m11, m12 ],
///  [ m21, m22 ]]
/// ```
#[derive(Clone, Hash, PartialEq, Eq, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Matrix2d<T: CoordinateType> {
    /// m11
    pub(crate) m11: T,
    /// m12
    pub(crate) m12: T,
    /// m21
    pub(crate) m21: T,
    /// m22
    pub(crate) m22: T,
}

impl<T> Matrix2d<T>
where
    T: CoordinateType,
{
    /// Create a new 2x2 matrix with entries of the form:
    /// ```txt
    /// [[ m11, m12 ],
    ///  [ m21, m22 ]]
    /// ```
    pub fn new(m11: T, m12: T, m21: T, m22: T) -> Self {
        Matrix2d { m11, m12, m21, m22 }
    }

    /// Return the identity matrix.
    pub fn identity() -> Self {
        Self::new(T::one(), T::zero(), T::zero(), T::one())
    }

    /// Compute product of the matrix with a scalar.
    pub fn mul_scalar(&self, rhs: T) -> Self {
        Matrix2d::new(
            self.m11 * rhs,
            self.m12 * rhs,
            self.m21 * rhs,
            self.m22 * rhs,
        )
    }

    /// Compute matrix-column-vector multiplication.
    /// The vector is interpreted as column vector.
    pub fn mul_column_vector(&self, rhs: Vector<T>) -> Vector<T> {
        Vector::new(
            rhs.x * self.m11 + rhs.y * self.m12,
            rhs.x * self.m21 + rhs.y * self.m22,
        )
    }

    /// Matrix-matrix multiplication.
    pub fn mul_matrix(&self, rhs: &Self) -> Self {
        let a = self;
        let b = rhs;
        let c11 = a.m11 * b.m11 + a.m12 * b.m21;
        let c12 = a.m11 * b.m12 + a.m12 * b.m22;
        let c21 = a.m21 * b.m11 + a.m22 * b.m21;
        let c22 = a.m21 * b.m12 + a.m22 * b.m22;
        Self::new(c11, c12, c21, c22)
    }

    /// Compute the transpose of the matrix.
    pub fn transpose(&self) -> Self {
        Self::new(self.m11, self.m21, self.m12, self.m22)
    }

    /// Compute the determinant of this matrix.
    pub fn determinant(&self) -> T {
        self.m11 * self.m22 - self.m12 * self.m21
    }

    /// Test if this matrix is the identity matrix.
    pub fn is_identity(&self) -> bool {
        self == &Self::identity()
    }

    /// Test if this matrix is unitary.
    pub fn is_unitary(&self) -> bool {
        self.mul_matrix(&self.transpose()).is_identity()
    }

    /// Compute the inverse matrix.
    /// `None` will be returned if the determinant is zero and the matrix
    /// is not invertible.
    pub fn try_inverse(&self) -> Option<Self> {
        // Compute determinant.
        let det = self.determinant();
        if !det.is_zero() {
            let z = T::zero();
            Some(Self::new(
                self.m22 / det,
                z - self.m12 / det,
                z - self.m21 / det,
                self.m11 / det,
            ))
        } else {
            None
        }
    }
}

impl<T: CoordinateType> Default for Matrix2d<T> {
    fn default() -> Self {
        Self::identity()
    }
}

#[test]
fn test_matrix_multiplication() {
    let a = Matrix2d::new(1.0, 2.0, 3.0, 4.0);
    let b = Matrix2d::new(5.0, 6.0, 7.0, 8.0);
    let id = Matrix2d::identity();
    assert_eq!(id.mul_matrix(&id), id);
    assert_eq!(b.mul_matrix(&id), b);
    assert_eq!(id.mul_matrix(&b), b);
    assert_eq!(
        a.mul_matrix(&b),
        Matrix2d::new(19.0, 22.0, 15.0 + 28.0, 18.0 + 32.0)
    );
}

#[test]
fn test_inverse() {
    let m = Matrix2d::new(2.0, 1.0, 4.0, 8.0);
    let i = m.try_inverse().unwrap();
    assert_eq!(m.mul_matrix(&i), Matrix2d::identity());
    assert_eq!(i.mul_matrix(&m), Matrix2d::identity());
}
