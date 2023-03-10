// Copyright (c) 2018-2020 Thomas Kramer.
// SPDX-FileCopyrightText: 2018-2022 Thomas Kramer
//
// SPDX-License-Identifier: AGPL-3.0-or-later

//! `Path` is essentially a chain of line segments but with a possibly non-zero width.
//! It can be thought of the shape resulting by a stroke of a thick pen along the line segments.

use crate::point::Point;
use crate::point_string::PointString;
use crate::vector::Vector;

use crate::CoordinateType;

pub use crate::traits::{BoundingBox, RotateOrtho};
use crate::traits::{MapPointwise, Scale, Translate, TryBoundingBox};
pub use crate::types::Angle;

pub use crate::types::{ContainsResult, Side};

use crate::edge::*;
use crate::rect::Rect;
use crate::simple_polygon::SimplePolygon;
use crate::transform::SimpleTransform;
use num_traits::{Float, Num, NumCast};
use std::iter::FromIterator;
use std::ops::{Add, Mul};

/// Encoding for the type of the beginning and end of the path.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum PathEndType<T> {
    /// Beginning and end of path are not extended.
    Flat,
    /// Define the extension length at the beginning and at the end of the path.
    Extended(T, T),
    /// Path ends are round (approximately semi-circles).
    Round,
}

/// `Path` is essentially a chain of line segments but with a possibly a non-zero width.
/// It can be thought of the shape resulting by a stroke of a thick pen along the line segments.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Path<T> {
    /// The vertices of the path which define the sequence of line segments.
    pub points: PointString<T>,
    /// Width of the path.
    pub width: T,
    /// Type of the path endings.
    pub path_type: PathEndType<T>,
}

impl<T> Path<T> {
    /// Get number of vertices defining the path.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Check if path has zero length.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }
}

impl<T: Copy> Path<T> {
    /// Create new path by taking vertices from a type that implements `Into<PointString<T>>`.
    pub fn new<I>(i: I, width: T) -> Self
    where
        I: Into<PointString<T>>,
    {
        Path {
            points: i.into(),
            width,
            path_type: PathEndType::Flat,
        }
    }

    /// Create a path with extended beginning and end.
    pub fn new_extended<I>(i: I, width: T, ext_begin: T, ext_end: T) -> Self
    where
        I: Into<PointString<T>>,
    {
        Path {
            points: i.into(),
            width,
            path_type: PathEndType::Extended(ext_begin, ext_end),
        }
    }

    /// Create a path with rounded beginning and end.
    pub fn new_rounded<I>(i: I, width: T) -> Self
    where
        I: Into<PointString<T>>,
    {
        Path {
            points: i.into(),
            width,
            path_type: PathEndType::Round,
        }
    }
}

impl<T: Copy + Add<Output = T>> Path<T> {
    /// Translate the path by an offset vector.
    pub fn translate(&self, v: Vector<T>) -> Self {
        Path {
            points: self.points.translate(v),
            width: self.width,
            path_type: self.path_type,
        }
    }
}

impl<T: Copy + Mul<Output = T>> Path<T> {
    /// Scale the path. Scaling center is the origin `(0, 0)`.
    pub fn scale(&self, factor: T) -> Self {
        Path {
            points: self.points.scale(factor),
            width: self.width * factor,
            path_type: self.path_type,
        }
    }
}

impl<T: CoordinateType> Path<T> {
    /// Rotate the path by a multiple of 90 degrees around the origin `(0, 0)`.
    pub fn rotate_ortho(&self, angle: Angle) -> Self {
        Path {
            points: self.points.rotate_ortho(angle),
            width: self.width,
            path_type: self.path_type,
        }
    }

    /// Get the transformed version of this path by applying `tf`.
    pub fn transform(&self, tf: &SimpleTransform<T>) -> Self {
        Self {
            points: self.points.transform(|p| tf.transform_point(p)),
            width: tf.transform_distance(self.width),
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
            _ => base_area,
        }
    }

    /// Convert the path into a polygon.
    /// The polygon can be self-intersecting.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use iron_shapes::prelude::*;
    /// let path = Path::new(&[(0, 0), (10, 0), (10, 20)], 4);
    /// let polygon = path.to_polygon_approx();
    /// assert_eq!(polygon, SimplePolygon::from(&[(0., 2.), (0., -2.), (12., -2.), (12., 20.), (8., 20.), (8., 2.)]));
    /// ```
    pub fn to_polygon_approx(&self) -> SimplePolygon<f64> {
        let mut points_forward: Vec<Point<f64>> = Vec::new();
        let mut points_backward: Vec<Point<f64>> = Vec::new();

        let edges: Vec<Edge<f64>> = self
            .points
            .edges()
            .filter(|e| e.start != e.end) // Skip zero-length edges.
            .map(|e| e.cast_to_float())
            .collect();

        // Construct rectangular start and end caps.
        let create_flat_cap =
            |center_edge: Edge<f64>, width: f64, extension: f64| -> Vec<Point<f64>> {
                let d = center_edge.vector().normalized();
                let n = d.rotate_ortho(Angle::R90);
                let p = center_edge.end;
                let w_half = width / 2.;
                if extension == 0. {
                    let p1 = p - n * w_half;
                    let p2 = p + n * w_half;
                    vec![p1, p2]
                } else {
                    let p1 = p - n * w_half;
                    let p4 = p + n * w_half;
                    let p2 = p1 + d * extension;
                    let p3 = p4 + d * extension;
                    vec![p1, p2, p3, p4]
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
            PathEndType::Round => unimplemented!("Not implemented for round path ends."),
        };

        // Path width.
        let width = NumCast::from(self.width).unwrap();
        let half_width = width * 0.5;

        // Create caps.
        let start_cap = edges
            .first()
            .map(|e| create_flat_cap(e.reversed(), width, start_ext))
            .unwrap_or_else(Vec::new);
        let end_cap = edges
            .last()
            .map(|e| create_flat_cap(*e, width, end_ext))
            .unwrap_or_else(Vec::new);

        // Pre-compute normals (scaled by half the width).
        let normals: Vec<Vector<f64>> = edges
            .iter()
            .map(|e| e.vector().normal() * half_width)
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
                LineIntersection::Point(p, _) => points_forward.push(p),
            }

            // Backward.
            let border_intersection_b = border1b.line_intersection_approx(&border2b, 1e-15);

            match border_intersection_b {
                LineIntersection::Collinear => {}
                LineIntersection::None => {}
                LineIntersection::Point(p, _) => points_backward.push(p),
            }
        }

        // Concatenate forward and backward points including start and end cap.
        SimplePolygon::from_iter(
            start_cap
                .iter()
                .chain(points_forward.iter())
                .chain(end_cap.iter())
                .chain(points_backward.iter().rev()),
        )
    }
}

impl<T: CoordinateType + NumCast, Dst: CoordinateType + NumCast> TryCastCoord<T, Dst> for Path<T> {
    type Output = Path<Dst>;

    fn try_cast(&self) -> Option<Self::Output> {
        let new_width = Dst::from(self.width);
        let new_points = self.points.try_cast();
        let new_path_type = match self.path_type {
            PathEndType::Extended(begin_ext, end_ext) => {
                let new_begin_ext = Dst::from(begin_ext);
                let new_end_ext = Dst::from(end_ext);
                match (new_begin_ext, new_end_ext) {
                    (Some(b), Some(e)) => Some(PathEndType::Extended(b, e)),
                    _ => None,
                }
            }
            PathEndType::Flat => Some(PathEndType::Flat),
            PathEndType::Round => Some(PathEndType::Round),
        };

        match (new_width, new_points, new_path_type) {
            (Some(width), Some(points), Some(path_type)) => Some(Path {
                points,
                width,
                path_type,
            }),
            _ => {
                // Failed to cast some values.
                None
            }
        }
    }
}

impl<T: Copy + PartialOrd + Num> TryBoundingBox<T> for Path<T> {
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
            let w_half = (self.width + _1) / _2;
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
