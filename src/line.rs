// Copyright (c) 2018-2020 Thomas Kramer.
// SPDX-FileCopyrightText: 2018-2022 Thomas Kramer
//
// SPDX-License-Identifier: AGPL-3.0-or-later

use crate::vector::Vector;
use crate::point::Point;

use crate::CoordinateType;

use std::f64;

/// Represents a line of infinite length in Hesse normal form.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Line<T>
    where T: CoordinateType {
    pub n: Vector<T>
}

impl<T: CoordinateType> Line<T>
{
    //    pub fn new()

    pub fn distance(self, &p: Point<T>) -> f64 {
        let l = self.n.length();
    }
}