// SPDX-FileCopyrightText: 2022 Thomas Kramer
//
// SPDX-License-Identifier: AGPL-3.0-or-later

//! Abstraction concepts over geometric primitives.
//!
//! # References
//! * Inspired by [the Boost polygon library](https://www.boost.org/doc/libs/1_78_0/libs/polygon/doc/gtl_design_overview.htm),
//! [archived](https://web.archive.org/web/20220524201016/https://www.boost.org/doc/libs/1_78_0/libs/polygon/doc/gtl_design_overview.htm)
//! * and ["Geometry Template Library for STL-like 2D Operations"](https://www.boost.org/doc/libs/1_44_0/libs/polygon/doc/GTL_boostcon2009.pdf),
//! [archived](https://web.archive.org/web/20220524201159/https://www.boost.org/doc/libs/1_44_0/libs/polygon/doc/GTL_boostcon2009.pdf)

use num_traits::{Num, Signed, Float};
use crate::isotropy::*;

/// Define the coordinate concept.
/// This will be used to parametrize all other geometric concepts.
pub trait CoordinateConcept {
    type Coord: Num + Signed + Copy + PartialEq;
    type Area: Num + Signed + Copy + PartialEq;
    type ManhattanArea: Num + Copy + PartialEq;
    type UnsignedArea: Num + Copy + PartialEq;
    type CoordinateDifference: Num + Signed + Copy + PartialEq + From<Self::Coord>;
    type CoordinateDistance: Num + Copy + PartialEq + From<Self::CoordinateDifference> + Float;
}

pub struct IntCoordinates;

impl CoordinateConcept for IntCoordinates {
    type Coord = i32;
    type Area = f64;
    type ManhattanArea = i64;
    type UnsignedArea = f64;
    type CoordinateDifference = i32;
    type CoordinateDistance = f64;
}

// macro_rules! crd {
//  ($t:ty) => {<<Self as $t>::Coord as CoordinateConcept>}
// }

pub trait PointConcept<C: CoordinateConcept> {

    /// Construct a new point.
    fn new(x: C::Coord, y: C::Coord) -> Self;

    /// Get a coordinate value.
    fn get(&self, orient: Orientation2D) -> C::Coord;

    /// Set a coordinate value.
    fn set(&mut self, orient: Orientation2D, value: C::Coord);

    /// Get the x-coordinate value.
    fn x(&self) -> C::Coord {
        self.get(Orientation2D::Horizontal)
    }

    /// Get the y-coordinate value.
    fn y(&self) -> C::Coord {
        self.get(Orientation2D::Vertical)
    }

    /// Compute the x or y component of the vector from the point to the `other` point.
    fn projected_distance(&self, other: &Self, orient: Orientation2D) -> C::CoordinateDifference {
        let a: C::CoordinateDifference = self.get(orient).into();
        let b: C::CoordinateDifference = other.get(orient).into();
        b - a
    }

    fn manhattan_distance(&self, other: &Self)  -> C::CoordinateDifference {
        self.projected_distance(other, Orientation2D::Horizontal) + self.projected_distance(other, Orientation2D::Vertical)
    }

    fn distance_squared(&self, other: &Self) -> C::CoordinateDistance {
        let d_horiz: C::CoordinateDistance = self.projected_distance(other, Orientation2D::Horizontal).into();
        let d_vert: C::CoordinateDistance = self.projected_distance(other, Orientation2D::Vertical).into();

        d_horiz * d_horiz + d_vert * d_vert
    }

    fn euclidian_distance(&self, other: &Self) -> C::CoordinateDistance {
        self.distance_squared(other).sqrt()
    }
}

/// A polygon edge.
pub trait Segment<C: CoordinateConcept> {
    type Point: PointConcept<C>;

    fn get_point(&self, dir: Direction1D) -> Self::Point;

    fn start(&self) -> Self::Point {
        self.get_point(Direction1D::Low)
    }

    fn end(&self) -> Self::Point {
        self.get_point(Direction1D::High)
    }
}

pub trait IntoSegments<C: CoordinateConcept> {
    type Segment: Segment<C>;
    type SegmentIter: Iterator<Item=Self::Segment>;

    fn into_segments(self) -> Self::SegmentIter;
}

pub trait IntoPoints<C: CoordinateConcept> {
    type Point: PointConcept<C>;
    type PointIter: Iterator<Item=Self::Point>;

    fn into_points(self) -> Self::PointIter;
}

pub trait Rectangle<C: CoordinateConcept>: Polygon90<C> {
    type Interval;

    fn get(&self, orientation: Orientation2D) -> Self::Interval;

    // fn lower_left(&self) -> <Self as PolygonSet>::Point;
    // fn upper_right(&self) -> <Self as PolygonSet>::Point;
    //
    // fn lower_right(&self) -> <Self as PolygonSet>::Point {
    //     let ll = self.lower_left();
    //     let ur = self.upper_right();
    //
    //     <Self as PolygonSet>::Point::new(ur.x(), ll.y())
    // }

}

pub trait Polygon90<C: CoordinateConcept>: Polygon<C> + Polygon90WithHoles<C> {
    /// Iterator over alternating x/y coordinates of the points. Starts with an x coordinate.
    type CompactIterator: Iterator<Item=<Self as PolygonSet<C>>::Point>;

    fn compact_iter(&self) -> Self::CompactIterator;

    fn set_from_compact(&mut self, iter: ()) {
        self.set(iter);
        todo!();
    }
}

pub trait Polygon<C: CoordinateConcept>: PolygonWithHoles<C> + IntoPoints<C> {

    /// Set points from an iterator.
    fn set(&mut self, iter: ());
}

pub trait Polygon90WithHoles<C: CoordinateConcept>:  PolygonWithHoles<C> + Polygon90Set<C> {

}

pub trait PolygonWithHoles<C: CoordinateConcept>: PolygonSet<C> {
    fn num_holes(&self) -> usize;
}

pub trait Polygon90Set<C: CoordinateConcept>: PolygonSet<C> {}

pub trait PolygonSet<C: CoordinateConcept>: IntoSegments<C> {
    type Point: PointConcept<C>;
    type Segment: Segment<C, Point=Self::Point>;
    // type Polygon: PolygonWithHoles;
    // type PolygonIter: Iterator<Item=Self::Polygon>;
    type AllPoints: Iterator<Item=Self::Point>;

    fn num_polygons(&self) -> usize;

    fn convolved(self, p: &Self::Point) -> Self;

    fn convolve(&mut self, p: &Self::Point);

    fn scaled(self, scale: C::Coord) -> Self;

    fn scale(&mut self, scale: C::Coord);

    fn all_points(&self) -> Self::AllPoints;
}