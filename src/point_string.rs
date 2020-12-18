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

//! A point string is a finite sequence of points.

use crate::point::Point;
use crate::edge::Edge;
use crate::rect::Rect;

use crate::CoordinateType;

use crate::traits::{MapPointwise, TryCastCoord};
pub use crate::traits::TryBoundingBox;

use std::iter::FromIterator;
use std::slice::Iter;

use num_traits::{Float, NumCast};
use itertools::Itertools;

/// A point string is a finite sequence of points.
/// TODO: Implement `Deref` for accessing the list of points.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct PointString<T: CoordinateType> {
    /// The points defining this point string.
    pub points: Vec<Point<T>>
}

impl<T: CoordinateType> PointString<T> {
    /// Create new point string by taking vertices from a type that implements `Into<PointString<T>>`.
    pub fn new<I>(i: I) -> Self
        where I: Into<Self> {
        i.into()
    }

    /// Get the number of vertices.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Shortcut for `self.points.iter()`.
    pub fn iter(&self) -> Iter<Point<T>> {
        self.points.iter()
    }

    /// Get the sequence of edges of the point string starting from the first point to the last.
    /// # Examples
    ///
    /// ```
    /// use iron_shapes::point_string::PointString;
    /// use iron_shapes::edge::Edge;
    /// let coords = vec![(0, 0), (1, 0), (2, 0)];
    ///
    /// let point_string = PointString::new(coords);
    ///
    /// let edges: Vec<_> = point_string.edges().collect();
    ///
    /// assert_eq!(edges, vec![Edge::new((0, 0), (1, 0)), Edge::new((1, 0), (2, 0))]);
    /// ```
    pub fn edges(&self) -> impl Iterator<Item=Edge<T>> + '_ {
        self.iter()
            .zip(self.iter().skip(1))
            .map(|(a, b)| Edge::new(a, b))
    }

    /// Same as `edges` but in reverse order.
    /// Get the sequence of edges of the point string starting from the last point to the first.
    /// # Examples
    ///
    /// ```
    /// use iron_shapes::point_string::PointString;
    /// use iron_shapes::edge::Edge;
    /// let coords = vec![(0, 0), (1, 0), (2, 0)];
    ///
    /// let point_string = PointString::new(coords);
    ///
    /// let edges: Vec<_> = point_string.edges_reversed().collect();
    ///
    /// assert_eq!(edges, vec![Edge::new((2, 0), (1, 0)), Edge::new((1, 0), (0, 0))]);
    /// ```
    pub fn edges_reversed(&self) -> impl Iterator<Item=Edge<T>> + '_ {
        self.iter().rev()
            .zip(self.iter().rev().skip(1))
            .map(|(a, b)| Edge::new(a, b))
    }
}

impl<T: CoordinateType + NumCast> PointString<T> {
    /// Compute geometrical length of the path defined by the point string.
    /// # Examples
    ///
    /// ```
    /// use iron_shapes::point_string::PointString;
    /// let coords = vec![(0, 0), (1, 0), (2, 0)];
    ///
    /// let point_string = PointString::new(coords);
    ///
    /// assert_eq!(point_string.path_length::<f64>(), 2.0);
    /// ```
    pub fn path_length<F: Float>(&self) -> F {
        self.iter()
            .zip(self.iter().skip(1))
            .map(|(a, b)| a.distance(b))
            .fold(F::zero(), |a, b| a + b)
    }
}

impl<T: CoordinateType + NumCast, Dst: CoordinateType + NumCast> TryCastCoord<T, Dst> for PointString<T> {
    type Output = PointString<Dst>;

    fn try_cast(&self) -> Option<Self::Output> {
        let new_points: Vec<_> = self.points.iter()
            .map(|p| p.try_cast())
            .while_some()
            .collect();
        if new_points.len() == self.points.len() {
            Some(PointString::new(new_points))
        } else {
            // Some points could not be casted.
            None
        }
    }
}

/// Create a point string from something that can be turned into an iterator of values convertible to [`Point`]s.
impl<I, T, P> From<I> for PointString<T>
    where T: CoordinateType,
          I: IntoIterator<Item=P>,
          Point<T>: From<P>
{
    fn from(iter: I) -> Self {
        let points: Vec<Point<T>> = iter.into_iter().map(
            |x| x.into()
        ).collect();

        PointString { points }
    }
}
//
// /// Create a point string from a [`Vec`] of values convertible to [`Point`]s.
// impl<'a, T, P> From<&'a Vec<P>> for PointString<T>
//     where T: CoordinateType,
//           Point<T>: From<&'a P>
// {
//     fn from(vec: &'a Vec<P>) -> Self {
//         let points: Vec<Point<T>> = vec.into_iter().map(
//             |x| x.into()
//         ).collect();
//
//         PointString { points }
//     }
// }
//
// /// Create a point string from a [`Vec`] of values convertible to [`Point`]s.
// impl<T, P> From<Vec<P>> for PointString<T>
//     where T: CoordinateType,
//           Point<T>: From<P>
// {
//     fn from(vec: Vec<P>) -> Self {
//         let points: Vec<Point<T>> = vec.into_iter().map(
//             |x| x.into()
//         ).collect();
//
//         PointString { points }
//     }
// }
//

/// Create a point string from a iterator of values convertible to [`Point`]s.
impl<T, P> FromIterator<P> for PointString<T>
    where T: CoordinateType,
          P: Into<Point<T>>
{
    fn from_iter<I>(iter: I) -> Self
        where I: IntoIterator<Item=P>
    {
        let points: Vec<Point<T>> = iter.into_iter().map(
            |x| x.into()
        ).collect();

        PointString { points }
    }
}

impl<T> MapPointwise<T> for PointString<T>
    where T: CoordinateType {
    fn transform<F: Fn(Point<T>) -> Point<T>>(&self, tf: F) -> Self {
        let points = self.points.iter()
            .map(|&p| tf(p))
            .collect();

        PointString {
            points
        }
    }
}

impl<T> TryBoundingBox<T> for PointString<T>
    where T: CoordinateType {
    /// Compute the bounding box of all the points in this string.
    /// Returns `None` if the string is empty.
    /// # Examples
    ///
    /// ```
    /// use iron_shapes::point_string::PointString;
    /// use iron_shapes::traits::TryBoundingBox;
    /// use iron_shapes::rect::Rect;
    /// let coords = vec![(0, 0), (1, 0), (2, 1), (-1, -3)];
    ///
    /// let point_string = PointString::new(coords);
    ///
    /// assert_eq!(point_string.try_bounding_box(), Some(Rect::new((2, 1), (-1, -3))));
    /// ```
    fn try_bounding_box(&self) -> Option<Rect<T>> {
        if self.points.is_empty() {
            None
        } else {
            let mut x_min = self.points[0].x;
            let mut x_max = x_min;
            let mut y_min = self.points[0].y;
            let mut y_max = y_min;

            for p in self.iter().skip(1) {
                let (x, y) = p.into();
                if x < x_min {
                    x_min = x;
                }
                if x > x_max {
                    x_max = x;
                }
                if y < y_min {
                    y_min = y;
                }
                if y > y_max {
                    y_max = y;
                }
            }

            Some(Rect::new((x_min, y_min), (x_max, y_max)))
        }
    }
}
