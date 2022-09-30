// Copyright (c) 2018-2020 Thomas Kramer.
// SPDX-FileCopyrightText: 2018-2022 Thomas Kramer
//
// SPDX-License-Identifier: AGPL-3.0-or-later

//! Compute the convex hull of a list of points.

use crate::prelude::{Edge, Point, Side, SimplePolygon};
use num_traits::Num;

/// Get the convex hull of the polygon.
///
/// Implements Andrew's Monotone Chain algorithm.
/// See: <http://geomalgorithms.com/a10-_hull-1.html>
pub fn convex_hull<T>(points: Vec<Point<T>>) -> SimplePolygon<T>
where
    T: Num + Ord + Copy,
{
    // Sort points by coordinates.
    let p = {
        let mut p = points;
        p.sort();
        p
    };

    assert!(p.len() >= 2);

    // minmin = index of P with min x first and min y second
    let minmin = 0;
    let p_minmin = p[minmin];
    // maxmax = index of P with max x first and max y second
    let maxmax = p.len() - 1;
    let p_maxmax = p[maxmax];
    // minmax = index of P with min x first and max y second
    let minmax = p
        .iter()
        .enumerate()
        .take_while(|(_, p)| p.x == p_minmin.x)
        .last()
        .unwrap()
        .0;
    let p_minmax = p[minmax];
    // maxmin = index of P with max x first and min y second
    let maxmin = p
        .iter()
        .enumerate()
        .rev()
        .take_while(|(_, p)| p.x == p_maxmax.x)
        .last()
        .unwrap()
        .0;
    let p_maxmin = p[maxmin];

    debug_assert!(minmin <= minmax);
    debug_assert!(minmax <= maxmin || p_minmax.x == p_maxmin.x);
    debug_assert!(maxmin <= maxmax);

    debug_assert!(p.iter().all(|&p| p_minmin <= p));
    debug_assert!(p.iter().all(|&p| p_maxmax >= p));
    debug_assert!(p
        .iter()
        .all(|&p| p_minmax.x < p.x || (p_minmax.x == p.x && p_minmax.y >= p.y)));
    debug_assert!(p
        .iter()
        .all(|&p| p_maxmin.x > p.x || (p_maxmin.x == p.x && p_maxmin.y <= p.y)));

    // Handle degenerate case where all x coordinates are equal.
    if p_minmin.x == p_maxmax.x {
        if p_minmin.y == p_maxmax.y {
            SimplePolygon::new(vec![p_minmin])
        } else {
            SimplePolygon::new(vec![p_minmin, p_maxmax])
        }
    } else {
        let build_half_hull = |l: Edge<T>, points: &[Point<T>]| {
            // Push starting point on stack.
            let mut stack = vec![l.start];

            for &pi in points {
                // Skip all points that are not strictly right of `l`.
                if l.side_of(pi) == Side::Right {
                    while stack.len() >= 2 {
                        let pt1 = stack[stack.len() - 1];
                        let pt2 = stack[stack.len() - 2];

                        if Edge::new(pt2, pt1).side_of(pi) == Side::Left {
                            // `pi` is strictly left of the line defined by the top two elements
                            // on the stack.
                            break;
                        }
                        stack.pop();
                    }
                    stack.push(pi);
                }
            }

            stack
        };

        // Compute the lower hull.
        let l_min = Edge::new(p_minmin, p_maxmin);
        let mut lower_half_hull = build_half_hull(l_min, &p[minmax + 1..maxmin]);

        // Push p_maxmin on stack if necessary.
        if p_maxmin != p_maxmax {
            lower_half_hull.push(p_maxmin);
        }

        // Compute the upper hull.
        let l_max = Edge::new(p_maxmax, p_minmax);
        let mut upper_half_hull = build_half_hull(l_max, &p[minmax + 1..maxmin]);

        // Push p_minmax on stack if necessary.
        if p_minmax != p_minmin {
            upper_half_hull.push(p_minmax);
        }

        // Join both hulls.
        lower_half_hull.append(&mut upper_half_hull);

        SimplePolygon::new(lower_half_hull)
    }
}
