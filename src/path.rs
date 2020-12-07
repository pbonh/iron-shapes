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
use crate::vector::Vector;
use crate::point::Point;
use crate::point_string::PointString;

use crate::CoordinateType;

use crate::traits::{Scale, Translate, TryBoundingBox};
pub use crate::traits::{BoundingBox, RotateOrtho};
pub use crate::types::Angle;

pub use crate::types::{Side, ContainsResult};

use num_traits::{Float, NumCast};
use crate::simple_polygon::SimplePolygon;
use std::iter::FromIterator;
use crate::edge::*;
use crate::rect::Rect;

/// Encoding for the type of the beginning and end of the path.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum PathEndType<T: CoordinateType> {
    /// Beginning and end of path are not extended.
    Flat,
    /// Define the extension length at the beginning and at the end of the path.
    Extended(T, T),
    /// Path ends are round (approximately semi-circles).
    Round,
}

///
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Path<T: CoordinateType> {
    pub points: PointString<T>,
    pub width: T,
    pub path_type: PathEndType<T>,
}

impl<T: CoordinateType> Path<T> {
    /// Create new path by taking vertices from a type that implements `Into<PointString<T>>`.
    pub fn new<I>(i: I, width: T) -> Self
        where I: Into<PointString<T>> {
        Path {
            points: i.into(),
            width,
            path_type: PathEndType::Flat,
        }
    }

    /// Create a path with extended beginning and end.
    pub fn new_extended<I>(i: I, width: T, ext_begin: T, ext_end: T) -> Self
        where I: Into<PointString<T>> {
        Path {
            points: i.into(),
            width,
            path_type: PathEndType::Extended(ext_begin, ext_end),
        }
    }

    /// Create a path with rounded beginning and end.
    pub fn new_rounded<I>(i: I, width: T) -> Self
        where I: Into<PointString<T>> {
        Path {
            points: i.into(),
            width,
            path_type: PathEndType::Round,
        }
    }

    /// Translate the path by an offset vector.
    pub fn translate(&self, v: Vector<T>) -> Self {
        Path {
            points: self.points.translate(v),
            width: self.width,
            path_type: self.path_type,
        }
    }

    /// Get number of vertices defining the path.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    pub fn scale(&self, factor: T) -> Self {
        Path {
            points: self.points.scale(factor),
            width: self.width * factor,
            path_type: self.path_type,
        }
    }

    pub fn rothate_ortho(&self, angle: Angle) -> Self {
        Path {
            points: self.points.rotate_ortho(angle),
            width: self.width,
            path_type: self.path_type,
        }
    }
}

impl<T: CoordinateType + NumCast> Path<T> {
    /// Compute approximate area occupied by the path.
    /// Simply computes length*width.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use iron_shapes::prelude::*;
    /// let path = Path::new(&[(0, 0), (0, 2)], 1);
    /// assert_eq!(path.area_approx::<f64>(), 2f64);
    /// ```
    pub fn area_approx<F: Float>(&self) -> F {
        let base_len: F = self.points.path_length();
        let w = F::from(self.width).unwrap();
        let l = match self.path_type {
            PathEndType::Extended(l1, l2) => base_len + F::from(l1 + l2).unwrap(),
            _ => base_len,
        };

        let base_area = l * w;

        // Add area of circle if path ends are round.
        match self.path_type {
            PathEndType::Round => base_area + F::from(std::f64::consts::PI).unwrap() * w * w,
            _ => base_area
        }
    }


    /// Convert the path into a polygon.
    /// The polygon can be self-intersecting.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use iron_shapes::prelude::*;
    /// let path = Path::new(&[(0, 0), (10, 0), (10, 20)], 2);
    /// let polygon = path.to_polygon_approx();
    /// assert_eq!(polygon, SimplePolygon::new(&[(0., 1.), (0., -1.), (11., -1.), (11., 20.), (9., 20.), (9., 1.)]));
    /// ```
    pub fn to_polygon_approx(&self) -> SimplePolygon<f64> {
        let mut points_forward: Vec<Point<f64>> = Vec::new();
        let mut points_backward: Vec<Point<f64>> = Vec::new();

        let edges: Vec<Edge<f64>> = self.points.edges()
            .map(|e| e.cast_to_float()).collect();


        // Construct rectangular start and end caps.
        let create_flat_cap = |center_edge: Edge<f64>, width: f64, extension: f64| -> Vec<Point<f64>> {
            let d = center_edge.vector().normalized();
            let n = d.rotate_ortho(Angle::R90);
            let p = center_edge.end;
            let w_half = width / 2.;
            if extension == 0. {
                let p1 = p - n * w_half;
                let p2 = p + n * w_half;
                vec![p1.into(), p2]
            } else {
                let p1 = p - n * w_half;
                let p2 = p1 + d * extension;
                let p4 = p + n * w_half;
                let p3 = p4 + d * extension;
                vec![p1.into(), p2.into(), p3, p4]
            }
        };

        // Calculate start/end extensions.
        let (start_ext, end_ext) = match self.path_type {
            PathEndType::Extended(start_ext, end_ext) => {
                let start_ext = NumCast::from(start_ext).unwrap();
                let end_ext = NumCast::from(end_ext).unwrap();
                (start_ext, end_ext)
            }
            PathEndType::Flat => (0., 0.),
            PathEndType::Round => unimplemented!("Not implemented for round path ends.")
        };

        // Path width.
        let width = NumCast::from(self.width).unwrap();

        // Create caps.
        let start_cap = edges.first()
            .map(|e| create_flat_cap(e.reversed(), width, start_ext))
            .unwrap_or_else(|| Vec::new());
        let end_cap = edges.last()
            .map(|e| create_flat_cap(*e, width, end_ext))
            .unwrap_or_else(|| Vec::new());

        // Pre-compute normals.
        let normals: Vec<Vector<f64>> = edges.iter()
            .map(|e| e.vector().normal())
            .collect();

        let edge_pairs = edges.iter().zip(edges.iter().skip(1));
        let normal_pairs = normals.iter().zip(normals.iter().skip(1));

        for ((&e1, &e2), (&n1, &n2)) in edge_pairs.zip(normal_pairs) {
            let border1f = e1.translate(-n1);
            let border1b = e1.translate(n1);
            let border2f = e2.translate(-n2);
            let border2b = e2.translate(n2);

            // Forward.
            let border_intersection_f = border1f.line_intersection_approx(&border2f, 1e-15);

            match border_intersection_f {
                LineIntersection::Collinear => {}
                LineIntersection::None => {}
                LineIntersection::Point(p, _) => {
                    points_forward.push(p)
                }
            }

            // Backward.
            let border_intersection_b = border1b.line_intersection_approx(&border2b, 1e-15);

            match border_intersection_b {
                LineIntersection::Collinear => {}
                LineIntersection::None => {}
                LineIntersection::Point(p, _) => {
                    points_backward.push(p)
                }
            }
        }

        // Concatenate forward and backward points including start and end cap.
        SimplePolygon::from_iter(
            start_cap.iter()
                .chain(points_forward.iter())
                .chain(end_cap.iter())
                .chain(points_backward.iter().rev())
        )
    }
}


impl<T: CoordinateType + NumCast, Dst: CoordinateType + NumCast> TryCastCoord<T, Dst> for Path<T> {
    type Output = Path<Dst>;

    fn try_cast(&self) -> Option<Self::Output> {
        if let Some(new_width) = Dst::from(self.width) {
            self.points.try_cast()
                .map(|ps| Path::new(ps, new_width))
        } else {
            // Failed to cast the width.
            None
        }
    }
}

impl<T: CoordinateType> TryBoundingBox<T> for Path<T> {
    // /// Compute the bounding box of this path.
    // fn bounding_box(&self) -> Rect<T> {
    //     // Compute the bounding box by first converting the path into a polygon
    //     // and then computing the bounding box of the polygon.
    //     // Since integer Paths do not support conversion to a polygon the path needs
    //     // to be converted to a float coordinate type.
    //     // TODO: Make this more efficient and preferably without type conversions.
    //     let float_path: Path<FloatType> = self.cast();
    //     let bbox = float_path.to_polygon_approx().bounding_box();
    //     bbox.cast()
    // }

    /// Compute the bounding box of this path.
    /// The returned bounding box is not necessarily the smallest bounding box.
    ///
    /// TODO: Find a better approximation.
    fn try_bounding_box(&self) -> Option<Rect<T>> {
        // Find the bounding box of a zero-width path.
        let bbox = self.points.try_bounding_box();
        let _1 = T::one();
        let _2 = _1 + _1;
        bbox.map(|bbox| {
            // Enlarge it by width/2 in all directions to make sure the path is fully contained
            // in the bounding box.
            let (p1, p2) = (bbox.lower_left(), bbox.upper_right());
            let w_half = (self.width+_1)/_2;
            let width = Vector::new(w_half, w_half);
            Rect::new(p1 - width, p2 + width)
        })
    }
}

// impl<T> Translate<T> for Path<T>
//     where T: CoordinateType {
//     fn translate(&self, v: Vector<T>) -> Self {
//         Path {
//             points: self.points.translate(v),
//             width: self.width,
//             path_type: self.path_type,
//         }
//     }
// }
//
// impl<T> Scale<T> for Path<T>
//     where T: CoordinateType {
//     fn scale(&self, factor: T) -> Self {
//         Path {
//             points: self.points.scale(factor),
//             width: self.width * factor,
//             path_type: self.path_type,
//         }
//     }
// }
