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
 
mod cmp;
pub mod math;
pub mod types;
pub mod traits;
pub mod vector;
pub mod point;
pub mod point_string;
pub mod edge;
pub mod edge_rational;
pub mod edge_integer;
pub mod rect;
pub mod simple_polygon;
pub mod simple_rpolygon;
pub mod polygon;
pub mod matrix2d;
pub mod matrix3d;
pub mod multi_polygon;
pub mod path;
pub mod shape;
pub mod text;
pub mod transform;
pub mod redge;
pub mod repetition;
pub mod algorithms;