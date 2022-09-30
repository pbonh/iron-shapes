// Copyright (c) 2018-2020 Thomas Kramer.
// SPDX-FileCopyrightText: 2018-2022 Thomas Kramer
//
// SPDX-License-Identifier: AGPL-3.0-or-later

use std::cmp::PartialOrd;

/// Return the smallest argument. If both are equal return `a`.
#[inline]
pub fn min<T: PartialOrd>(a: T, b: T) -> T {
    if a <= b {
        a
    } else {
        b
    }
}

/// Return the largest argument. If both are equal return `a`.
#[inline]
pub fn max<T: PartialOrd>(a: T, b: T) -> T {
    if a >= b {
        a
    } else {
        b
    }
}
