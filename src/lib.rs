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
pub mod polygon;
pub mod matrix2d;
pub mod matrix3d;
pub mod multi_polygon;
pub mod path;
pub mod shape;
pub mod transform;
