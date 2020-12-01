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

use std::ops::Mul;
use std::convert::identity;

#[derive(Clone, Hash, PartialEq, Eq, Debug)]
pub struct Matrix2d<T: CoordinateType> {
    rows: [Vector<T>; 2]
}

impl<T> Matrix2d<T>
    where T: CoordinateType
{
    /// Create a new 2x2 matrix with entries of the form:
    /// ```txt
    /// [[ a, b ],
    ///  [ c, c ]]
    /// ```
    pub fn new(a: T, b: T, c: T, d: T) -> Self {
        Matrix2d {
            rows: [Vector::new(a, b), Vector::new(c, d)]
        }
    }

    /// Return the identity matrix.
    pub fn identity() -> Self {
        Self::new(T::one(), T::zero(), T::zero(), T::one())
    }

    /// Compute product of the matrix with a scalar.
    pub fn mul_scalar(&self, rhs: T) -> Self {
        Matrix2d {
            rows: [self.rows[0] * rhs, self.rows[1] * rhs]
        }
    }

    /// Compute matrix-vector multiplication.
    /// The vector is interpreted as column vector.
    pub fn mul_vector(&self, rhs: Vector<T>) -> Vector<T> {
        Vector {
            x: self.rows[0].dot(rhs),
            y: self.rows[1].dot(rhs),
        }
    }

    /// Matrix-matrix multiplication.
    pub fn mul(&self, rhs: &Self) -> Self {
        // TODO: Make more efficient.
        let t = rhs.transpose();
        Matrix2d {
            rows: [self.mul_vector(t.rows[0]),
                self.mul_vector(t.rows[1])]
        }.transpose()
    }

    /// Compute the transpose of the matrix.
    pub fn transpose(&self) -> Self {
        Matrix2d::new(self.rows[0].x, self.rows[1].x,
                      self.rows[0].y, self.rows[1].y)
    }

    /// Compute the determinant of this matrix.
    pub fn determinant(&self) -> T {
        let a = self.rows[0].x;
        let b = self.rows[0].y;
        let c = self.rows[1].x;
        let d = self.rows[1].y;
        a * d - b * c
    }

    /// Compute the inverse matrix.
    /// `None` will be returned if the determinant is zero and the matrix
    /// is not invertible.
    pub fn try_inverse(&self) -> Option<Self> {
        let a = self.rows[0].x;
        let b = self.rows[0].y;
        let c = self.rows[1].x;
        let d = self.rows[1].y;

        // Compute determinant.
        let det = a * d - b * c;
        if !det.is_zero() {
            let z = T::zero();
            Some(Self::new(d / det, z - b / det, z - c / det, a / det))
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
    assert_eq!(id.mul(&id), id);
    assert_eq!(b.mul(&id), b);
    assert_eq!(id.mul(&b), b);
    assert_eq!(a.mul(&b), Matrix2d::new(19.0, 22.0, 15.0+28.0, 18.0+32.0));
}

#[test]
fn test_inverse() {
    let m = Matrix2d::new(2.0, 1.0, 4.0, 8.0);
    let i = m.try_inverse().unwrap();
    assert_eq!(m.mul(&i), Matrix2d::identity());
    assert_eq!(i.mul(&m), Matrix2d::identity());
}