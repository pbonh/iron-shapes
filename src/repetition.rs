// Copyright (c) 2018-2020 Thomas Kramer.
// SPDX-FileCopyrightText: 2018-2022 Thomas Kramer
//
// SPDX-License-Identifier: AGPL-3.0-or-later

//! Describe repetitions of geometrical objects.
//!
//! Regular repetitions (arrays) are based on two lattice vectors, irregular repetitions are based
//! on a list of offsets.

use crate::prelude::{CoordinateType, Vector};
use num_traits::Zero;

/// Describe a equi-spaced n*m two-dimensional repetition as a lattice.
/// The offsets are computed as `(i*a, j*b)` for `i` in `0..n` and `j` in `0..m`.
/// `a` and `b` the distance vectors between two neighbouring points.
#[derive(PartialEq, Eq, Copy, Clone, Debug, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct RegularRepetition<T>
    where T: CoordinateType {
    /// First lattice vector.
    a: Vector<T>,
    /// Second lattice vector.
    b: Vector<T>,
    /// First dimension.
    n: u32,
    /// Second dimension.
    m: u32,
}

impl<T> RegularRepetition<T>
    where T: CoordinateType {
    /// Create a new lattice based repetition.
    ///
    /// # Parameters
    /// * `a, b`: Lattice vectors.
    /// * `n, m`: Number of repetitions in directions `a` and `b`.
    pub fn new(a: Vector<T>, b: Vector<T>, n: u32, m: u32) -> Self {
        return RegularRepetition { a, b, n, m };
    }

    /// Create a repetition along the x and y axis.
    ///
    /// # Example
    ///
    /// ```
    /// use iron_shapes::prelude::RegularRepetition;
    ///
    /// let rep = RegularRepetition::new_rectilinear(1, 1, 1, 2);
    /// assert_eq!(rep.len(), 2);
    /// let offsets: Vec<_> = rep.iter().collect();
    ///
    /// assert_eq!(offsets, [(0, 0).into(), (0, 1).into()]);
    /// ```
    pub fn new_rectilinear(spacing_x: T, spacing_y: T, num_x: u32, num_y: u32) -> Self {
        Self::new(Vector::new(spacing_x, T::zero()), Vector::new(T::zero(), spacing_y),
                  num_x, num_y)
    }

    /// Iterate over each offsets of this repetition.
    pub fn iter(self) -> impl Iterator<Item=Vector<T>> {
        let mut current = Vector::zero();
        (0..self.m).flat_map(move |_| {
            let mut row = current;
            current = current + self.b;

            (0..self.n).map(move |_| {
                let pos = row;
                row = row + self.a;
                pos
            })
        })
    }

    /// Return the number of offsets in this repetition.
    pub fn len(&self) -> usize {
        (self.n as usize) * (self.m as usize)
    }
}

/// Describe a non-equispaced repetition by storing a list of offsets.
#[derive(PartialEq, Eq, Clone, Debug, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct IrregularRepetition<T>
    where T: CoordinateType {
    /// Offset vectors of the repetition.
    offsets: Vec<Vector<T>>
}

impl<T> IrregularRepetition<T>
    where T: CoordinateType {

    /// Create a new irregular repetition from a list of offsets.
    pub fn new(offsets: Vec<Vector<T>>) -> Self {
        return IrregularRepetition { offsets };
    }

    /// Iterate over each offsets of this repetition.
    pub fn iter(&self) -> impl Iterator<Item=Vector<T>> + '_ {
        self.offsets.iter().copied()
    }

    /// Return the number of offsets in this repetition.
    pub fn len(&self) -> usize {
        self.offsets.len()
    }
}

/// Describe the regular or irregular repetition of a geometrical object.
#[derive(PartialEq, Eq, Clone, Debug, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum Repetition<T>
    where T: CoordinateType {
    /// Lattice based repetition.
    Regular(RegularRepetition<T>),
    /// Repetition with random offsets.
    Irregular(IrregularRepetition<T>),
}

impl<T> Repetition<T>
    where T: CoordinateType {
    /// Return the number of offsets in this repetition.
    pub fn len(&self) -> usize {
        match self {
            Repetition::Regular(r) => r.len(),
            Repetition::Irregular(r) => r.len()
        }
    }


    /// Iterate over each offsets of this repetition.
    pub fn iter(&self) -> Box<dyn Iterator<Item=Vector<T>> + '_> {
        match self {
            Repetition::Regular(r) => Box::new(r.iter()),
            Repetition::Irregular(r) => Box::new(r.iter())
        }
    }
}