// SPDX-FileCopyrightText: 2022 Thomas Kramer
//
// SPDX-License-Identifier: AGPL-3.0-or-later

//! Implement the geometry concepts for the default data types.

use crate::concepts::*;
use crate::isotropy::*;
use crate::prelude as types;
use crate::prelude::{IntoEdges, Translate, Scale, CoordinateType};


/// Coordinate types based on `i32` as main coordinate type.
pub struct I32Coordinates;

impl CoordinateBase for I32Coordinates {
    type Coord = i32;
}

impl CoordinateConcept for I32Coordinates {
    type Area = f64;
    type ManhattanArea = i64;
    type UnsignedArea = f64;
    type CoordinateDifference = i32;
    type CoordinateDistance = f64;
}

impl<C> PointBase<C> for types::Point<C::Coord>
    where C: CoordinateBase {
    fn new(x: C::Coord, y: C::Coord) -> Self {
        types::Point::new(x, y)
    }

    fn get(&self, orient: Orientation2D) -> C::Coord {
        match orient {
            Orientation2D::Horizontal => self.x,
            Orientation2D::Vertical => self.y
        }
    }

    fn set(&mut self, orient: Orientation2D, value: C::Coord) {
        let c = match orient {
            Orientation2D::Horizontal => &mut self.x,
            Orientation2D::Vertical => &mut self.y,
        };
        *c = value
    }
}

impl<C> PointConcept<C> for types::Point<C::Coord>
    where C: CoordinateConcept {}


impl<C> Segment<C> for types::REdge<C::Coord>
    where C: CoordinateConcept {
    type Point = types::Point<C::Coord>;

    fn get_point(&self, dir: Direction1D) -> Self::Point {
        match dir {
            Direction1D::Low => self.start(),
            Direction1D::High => self.end()
        }
    }
}

impl<C> Interval<C> for types::Interval<C>
    where C: CoordinateType {
    fn get(&self, d: Direction1D) -> C {
        match d {
            Direction1D::Low => self.start(),
            Direction1D::High => self.end()
        }
    }
}

impl<C> Rectangle<C> for types::Rect<C::Coord>
    where C: CoordinateConcept {
    type Interval = types::Interval<C::Coord>;

    fn get(&self, orientation: Orientation2D) -> Self::Interval {
        // let start = self.lower_left.get(orientation);
        // let end = self.upper_right.get(orientation);
        let (start, end) = match orientation {
            Orientation2D::Horizontal => (self.lower_left().x, self.upper_right().x),
            Orientation2D::Vertical => (self.lower_left().y, self.upper_right().y)
        };
        types::Interval::new(start, end)
    }
}

impl<C> Polygon90<C> for types::Rect<C::Coord>
    where C: CoordinateConcept {
    type CompactIterator = std::vec::IntoIter<C::Coord>;

    fn compact_iter(&self) -> Self::CompactIterator {
        vec![self.lower_left.x, self.lower_left.y, self.upper_right.x, self.upper_right.y]
            .into_iter()
    }
}

impl<C> Polygon<C> for types::Rect<C::Coord>
    where C: CoordinateConcept {
    fn set(&mut self, iter: ()) {
        unimplemented!()
    }
}


impl<C> PolygonWithHoles<C> for types::Rect<C::Coord>
    where C: CoordinateConcept {
    fn num_holes(&self) -> usize {
        0
    }
}

impl<C> Polygon90WithHoles<C> for types::Rect<C::Coord>
    where C: CoordinateConcept {}

impl<C> Polygon90Set<C> for types::Rect<C::Coord>
    where C: CoordinateConcept {}


impl<C> IntoSegments<C> for types::Rect<C::Coord>
    where C: CoordinateConcept {
    type Segment = types::REdge<C::Coord>;
    type SegmentIter = types::RectEdgeIterator<C::Coord>;

    fn into_segments(self) -> Self::SegmentIter {
        self.into_edges()
    }
}

impl<C> IntoPoints<C> for types::Rect<C::Coord>
    where C: CoordinateConcept {
    type Point = types::Point<C::Coord>;
    type PointIter = <Self as IntoIterator>::IntoIter;

    fn into_points(self) -> Self::PointIter {
        self.into_iter()
    }
}

impl<C> PolygonSet<C> for types::Rect<C::Coord>
    where C: CoordinateConcept {
    type Point = types::Point<C::Coord>;
    type Segment = types::REdge<C::Coord>;
    type AllPoints = <Self as IntoIterator>::IntoIter;

    fn num_polygons(&self) -> usize {
        1
    }

    fn convolved(mut self, p: &Self::Point) -> Self {
        // self.convolve(p);
        // self
        todo!()
    }

    fn convolve(&mut self, p: &Self::Point) {
        // let v: types::Vector<C::Coord> = p.into();
        // self.upper_right = self.upper_right + v;
        // self.lower_left = self.lower_left + v;
        todo!()
    }

    fn scaled(mut self, scale: C::Coord) -> Self {
        self.scale(scale);
        self
    }

    fn scale(&mut self, scale: C::Coord) {
        todo!()
    }

    fn all_points(&self) -> Self::AllPoints {
        self.into_iter()
    }
}

#[test]
fn test_point() {

    fn some_point_function<P: PointBase<C>, C: CoordinateBase>(p: P) -> C::Coord {
        p.x() + p.y()
    }

    let p = types::Point::new(1, 2);
    assert_eq!(some_point_function::<_, i32>(p), 3);

}