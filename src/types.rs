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

//! Commonly used type definitions and constants.

use std::ops::{Add, Neg};

/// Precision for distance related decisions.
pub const PREC_DISTANCE: DistanceType = 1e-5;

/// Default floating point type.
pub type FloatType = f64;

/// Default type for euclidean distances.
pub type DistanceType = FloatType;

/// Angle expressed as a multiple of 90 degrees.
#[derive(Copy, Clone, Hash, PartialEq, Eq, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum Angle {
    /// 0 Degrees.
    R0,
    /// 90 Degrees.
    R90,
    /// 180 Degrees.
    R180,
    /// 270 Degrees.
    R270,
}

impl Angle {
    /// Describe the angle as a integer multiple of 90 degrees.
    pub fn as_int(&self) -> u32 {
        match self {
            Angle::R0 => 0,
            Angle::R90 => 1,
            Angle::R180 => 2,
            Angle::R270 => 3,
        }
    }

    /// Convert an integer to an angle.
    /// The integer specifies the number of 90 degree rotations.
    pub fn from_u32(a: u32) -> Self {
        match a % 4 {
            0 => Angle::R0,
            1 => Angle::R90,
            2 => Angle::R180,
            3 => Angle::R270,
            _ => panic!()
        }
    }
}

impl Add for Angle {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::from_u32(self.as_int() + rhs.as_int())
    }
}

impl Neg for Angle {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::from_u32(4 - self.as_int())
    }
}

impl Default for Angle {
    fn default() -> Self {
        Angle::R0
    }
}

/// Location relative to a directed line or edge.
/// Something can be on the left of the line, right of the line or on top of the line (center).
#[derive(Clone, Copy, Hash, PartialEq, Debug)]
pub enum Side {
    /// Location on the left side.
    Left,
    /// Neither on the left nor right, but on top.
    Center,
    /// Location on the right side.
    Right,
}

impl Side {
    /// Test if this is the left side.
    pub fn is_left(&self) -> bool {
        *self == Side::Left
    }
    /// Test if this is the right side.
    pub fn is_right(&self) -> bool {
        *self == Side::Right
    }
    /// Test if this is `Center`.
    pub fn is_center(&self) -> bool {
        *self == Side::Center
    }

    /// Get the other side.
    pub fn other(&self) -> Self {
        match self {
            Side::Left => Side::Right,
            Side::Center => Side::Center,
            Side::Right => Side::Left
        }
    }
}

/// Relative orientation of two geometrical objects such as vectors.
#[derive(Clone, Copy, Hash, PartialEq, Debug)]
pub enum Orientation {
    /// Clock-wise orientation.
    ClockWise,
    /// Counter-clock-wise orientation.
    CounterClockWise,
    /// Neither clock-wise nor counter-clock-wise.
    Straight,
}

impl Orientation {
    /// Test if the orientation is equal to `ClockWise`.
    pub fn is_clock_wise(self) -> bool {
        self == Self::ClockWise
    }
    /// Test if the orientation is equal to `CounterClockWise`.
    pub fn is_counter_clock_wise(self) -> bool {
        self == Self::CounterClockWise
    }
    /// Test if the orientation is equal to `Straight`.
    pub fn is_straight(self) -> bool {
        self == Self::Straight
    }
}

/// This is a result type for containment checks.
/// * `No` Not inside.
/// * `OnBounds` Lies on the boundaries.
/// * `WithinBounds` Fully inside.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum ContainsResult {
    /// Does not contain the point.
    No,
    /// Contains the point but on the borders/end-points.
    OnBounds,
    /// Fully contains the point.
    WithinBounds,
}

impl ContainsResult {
    /// Tells if the point is contained but does not lie on the bounds.
    pub fn is_within_bounds(&self) -> bool {
        match self {
            ContainsResult::WithinBounds => true,
            _ => false
        }
    }

    /// Tells if the point is contained or lies on the bounds.
    pub fn inclusive_bounds(&self) -> bool {
        match self {
            ContainsResult::WithinBounds | ContainsResult::OnBounds => true,
            _ => false
        }
    }

    /// Check if the point neither is on the bounds nor within the bounds.
    pub fn on_bounds(&self) -> bool {
        match self {
            ContainsResult::OnBounds => true,
            _ => false
        }
    }

    /// Check if the point lies on the bounds.
    pub fn is_no(&self) -> bool {
        match self {
            ContainsResult::No => true,
            _ => false
        }
    }

    /// Returns the stronger result of the both.
    /// Ordering from weak to strong is `No`, `OnBounds`, `WithinBounds`
    pub fn max(self, other: Self) -> Self {
        match (self, other) {
            (ContainsResult::WithinBounds, _) => ContainsResult::WithinBounds,
            (_, ContainsResult::WithinBounds) => ContainsResult::WithinBounds,
            (ContainsResult::OnBounds, _) => ContainsResult::OnBounds,
            (_, ContainsResult::OnBounds) => ContainsResult::OnBounds,
            (ContainsResult::No, _) => ContainsResult::No,
        }
    }

    /// Returns the weaker result of the both.
    /// Ordering from weak to strong is `No`, `OnBounds`, `WithinBounds`
    pub fn min(self, other: Self) -> Self {
        match (self, other) {
            (ContainsResult::No, _) => ContainsResult::No,
            (_, ContainsResult::No) => ContainsResult::No,
            (ContainsResult::OnBounds, _) => ContainsResult::OnBounds,
            (_, ContainsResult::OnBounds) => ContainsResult::OnBounds,
            (ContainsResult::WithinBounds, _) => ContainsResult::WithinBounds,
        }
    }
}