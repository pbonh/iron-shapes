// Copyright (c) 2018-2020 Thomas Kramer.
// SPDX-FileCopyrightText: 2018-2022 Thomas Kramer
//
// SPDX-License-Identifier: AGPL-3.0-or-later

//! This crate provides basic data structures for Euclidean geometry in the plane.

// Enforce documentation.
#![deny(missing_docs)]

#[cfg(feature = "serde")]
#[macro_use]
extern crate serde;
extern crate itertools;

use num_traits::Num;

/// Trait for types that can be used as coordinates in the euclidean plane.
/// In practice this are integers, floats and possible rational numbers.
pub trait CoordinateType: Num + Copy + PartialOrd + Sized {}
impl<T: Num + Copy + PartialOrd> CoordinateType for T {}

/// Trait for types that can be used as areas.
/// In practice this are integers, floats and possible rational numbers.
pub trait AreaType: Num + Copy + PartialOrd + Sized {}
impl<T: Num + Copy + PartialOrd> AreaType for T {}

pub mod prelude;

pub mod algorithms;
mod cmp;
pub mod concepts;
pub mod concepts_impl;
pub mod edge;
pub mod edge_integer;
pub mod edge_rational;
pub mod interval;
pub mod isotropy;
pub mod math;
pub mod matrix2d;
pub mod matrix3d;
pub mod multi_polygon;
pub mod path;
pub mod point;
pub mod point_string;
pub mod polygon;
pub mod rect;
pub mod redge;
pub mod repetition;
pub mod shape;
pub mod simple_polygon;
pub mod simple_rpolygon;
pub mod text;
pub mod traits;
pub mod transform;
pub mod types;
pub mod vector;
