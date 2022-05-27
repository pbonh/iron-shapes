// SPDX-FileCopyrightText: 2018-2022 Thomas Kramer
//
// SPDX-License-Identifier: AGPL-3.0-or-later


//! A one dimensional interval.

use crate::CoordinateType;
use crate::cmp::{max, min};

/// A one dimensional interval which is represented by a start and end coordinate.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Interval<T>(Option<(T, T)>)
    where T: CoordinateType;

impl<T: CoordinateType> Interval<T> {

    /// Create a new interval.
    pub fn new(lower: T, upper: T) -> Self {
        if lower <= upper {
            Interval(Some((lower, upper)))
        } else {
            Interval(None)
        }
    }

    /// Compute the intersection of two intervals.
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

    /// Test if the interval is empty.
    pub fn is_empty(&self) -> bool {
        self.0.is_none()
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