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
use crate::cmp::{max, min};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Interval<T>(Option<(T, T)>)
    where T: CoordinateType;

impl<T: CoordinateType> Interval<T> {
    pub fn new(lower: T, upper: T) -> Self {
        if lower <= upper {
            Interval(Some((lower, upper)))
        } else {
            Interval(None)
        }
    }

    pub fn intersection(&self, other: &Self) -> Self {
        match self {
            Interval(None) => Interval(None),
            Interval(Some((l1, u1))) => match other {
                Interval(None) => Interval(None),
                Interval(Some((l2, u2))) =>
                    Interval::new(max(*l1, *l2), min(*u1, *u2))
            }
        }
    }
}

#[test]
fn test_interval_intersection() {
    let a = Interval::new(0, 10);
    let b = Interval::new(5, 15);
    assert_eq!(a.intersection(&b), Interval::new(5, 10));

    let a = Interval::new(0, 10);
    let b = Interval::new(11, 15);
    assert_eq!(a.intersection(&b), Interval(None));
}