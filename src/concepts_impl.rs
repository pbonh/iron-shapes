// SPDX-FileCopyrightText: 2022 Thomas Kramer
//
// SPDX-License-Identifier: AGPL-3.0-or-later

//! Implement the geometry concepts for the default data types.

use crate::concepts::*;
use crate::isotropy::*;
use crate::prelude::Point;

impl PointConcept<I32Coordinates> for Point<i32> {

    fn new(x: i32, y: i32) -> Self {
        Point::new(x, y)
    }

    fn get(&self, orient: Orientation2D) -> i32 {
        match orient {
            Orientation2D::Horizontal => self.x,
            Orientation2D::Vertical => self.y
        }
    }

    fn set(&mut self, orient: Orientation2D, value: i32) {
        let c = match orient {
            Orientation2D::Horizontal => &mut self.x,
            Orientation2D::Vertical => &mut self.y,
        };
        *c = value
    }
}

/// Coordinate types based on `i32` as main coordinate type.
pub struct I32Coordinates;

impl CoordinateConcept for I32Coordinates {
    type Coord = i32;
    type Area = f64;
    type ManhattanArea = i64;
    type UnsignedArea = f64;
    type CoordinateDifference = i32;
    type CoordinateDistance = f64;
}