// SPDX-FileCopyrightText: 2022 Thomas Kramer
//
// SPDX-License-Identifier: AGPL-3.0-or-later

//! Describe alignment and orientation along coordinate axes.

use num_traits::Signed;

/// Direction in one dimension.
#[derive(Copy, Clone, Debug, Hash, Ord, PartialOrd, Eq, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum Direction1D {
    /// Towards negative coordinates.
    Low,
    /// Towards positive coordinates.
    High,
}


impl Direction1D {
    /// Return the opposite direction.
    pub fn reversed(self) -> Self {
        match self {
            Direction1D::Low => Direction1D::High,
            Direction1D::High => Direction1D::Low
        }
    }

    /// Map `Low` to `-1` and `High` to `+1`.
    pub fn sign<N>(&self) -> N
        where N: Signed {
        match self {
            Direction1D::Low => N::one().neg(),
            Direction1D::High => N::one()
        }
    }
}

impl From<Direction2D> for Direction1D {
    /// Downcast 2D direction to 1D direction by preserving the sign.
    fn from(d2: Direction2D) -> Self {
        match d2 {
            Direction2D::Left | Direction2D::Up => Direction1D::High,
            Direction2D::Right | Direction2D::Down => Direction1D::Low,
        }
    }
}

/// Orientation in two dimensions along one of the axes.
#[derive(Copy, Clone, Debug, Hash, Ord, PartialOrd, Eq, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum Orientation2D {
    /// Oriented parallel to the x axis.
    Horizontal,
    /// Oriented parallel to the y axis.
    Vertical,
}

impl Orientation2D {

    /// Get the orthogonal orientation.
    pub fn other(self) -> Self {
        match self {
            Orientation2D::Horizontal => Self::Vertical,
            Orientation2D::Vertical => Self::Horizontal
        }
    }

    /// Get the positive or negative direction aligned with this orientation.
    pub fn direction(self, d1: Direction1D) -> Direction2D {
        use Orientation2D::*;
        use Direction1D::*;
        use Direction2D::*;
        match (self, d1) {
            (Horizontal, Low) => Right,
            (Horizontal, High) => Left,
            (Vertical, Low) => Down,
            (Vertical, High) => Up
        }
    }

    /// Check if orientation is vertical.
    pub fn is_vertical(&self) -> bool {
        self == &Self::Vertical
    }

    /// Check if orientation is horizontal.
    pub fn is_horizontal(&self) -> bool {
        self == &Self::Horizontal
    }
}

/// Directions along the coordinate axes in two dimensions.
#[derive(Copy, Clone, Debug, Hash, Ord, PartialOrd, Eq, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum Direction2D {
    /// Positive x direction.
    Left,
    /// Positive y direction.
    Up,
    /// Negative x direction.
    Right,
    /// Negative y direction.
    Down,
}

impl Direction2D {

    /// Get the opposite direction.
    pub fn reversed(self) -> Self {
        match self {
            Direction2D::Left => Self::Right,
            Direction2D::Up => Self::Down,
            Direction2D::Right => Self::Left,
            Direction2D::Down => Self::Up
        }
    }

    /// Get the orientation of this direction.
    pub fn orientation(&self) -> Orientation2D {
        match self {
            Direction2D::Left | Direction2D::Right => Orientation2D::Horizontal,
            Direction2D::Up | Direction2D::Down => Orientation2D::Vertical,
        }
    }

    /// Get the sign of the direction.
    pub fn sign<N>(self) -> N
        where N: Signed {
        let d1: Direction1D = self.into();
        d1.sign()
    }
}