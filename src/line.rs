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
use crate::vector::Vector;
use crate::point::Point;

use crate::CoordinateType;

use std::f64;

/// Represents a line of infinite length in Hesse normal form.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
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