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

use crate::isotropy::*;
use num_traits::{Float, Num, Signed};

/// Define a type used for coordinates.
pub trait CoordinateBase {
    /// Base coordinate type used for x and y coordinates.
    type Coord: Num + Signed + Copy + PartialEq + PartialOrd;
}

impl<T> CoordinateBase for T
where
    T: Num + Signed + Copy + PartialEq + PartialOrd,
{
    type Coord = T;
}

/// Define the coordinate concept.
/// This will be used to parametrize all other geometric concepts.
pub trait CoordinateConcept: CoordinateBase {
    /// Type used for area. This is typically a floating point number or a rational number.
    type Area: Num + Signed + Copy + PartialEq + PartialOrd;
    /// Type used for area which can be expressed without fractions. The datatype usually has a bigger range than `Coord` to avoid overflows during multiplications.
    /// For example when using `i32` as `Coord`, a `i64` is recommended as area type.
    type ManhattanArea: Num + Copy + PartialEq + PartialOrd;
    /// Type for unsigned area.
    type UnsignedArea: Num + Copy + PartialEq + PartialOrd;
    /// Type for difference between coordinates. Typically the same type as `Coord` when `Coord` is signed.
    type CoordinateDifference: Num + Signed + Copy + PartialEq + From<Self::Coord> + PartialOrd;
    /// Type for distances. Typically a floating point type because distances cannot be represented in integers nor rationals.
    type CoordinateDistance: Num
        + Copy
        + PartialEq
        + From<Self::CoordinateDifference>
        + Float
        + PartialOrd;
}

// macro_rules! crd {
//  ($t:ty) => {<<Self as $t>::Coord as CoordinateConcept>}
// }

/// Basic traits to get and set Kartesian coordinates of a point in the two-dimensional plane.
pub trait PointBase<C: CoordinateBase> {
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
}

/// Concept of a point in the Euclidean plane.
pub trait PointConcept<C: CoordinateConcept>: PointBase<C> {
    /// Compute the x or y component of the vector from the point to the `other` point.
    fn projected_distance(&self, other: &Self, orient: Orientation2D) -> C::CoordinateDifference {
        let a: C::CoordinateDifference = self.get(orient).into();
        let b: C::CoordinateDifference = other.get(orient).into();
        b - a
    }

    /// Compute the 1-norm of the vector pointing from the point to the other.
    fn manhattan_distance(&self, other: &Self) -> C::CoordinateDifference {
        self.projected_distance(other, Orientation2D::Horizontal)
            + self.projected_distance(other, Orientation2D::Vertical)
    }

    /// Squared euclidean distance.
    fn distance_squared(&self, other: &Self) -> C::CoordinateDistance {
        let d_horiz: C::CoordinateDistance = self
            .projected_distance(other, Orientation2D::Horizontal)
            .into();
        let d_vert: C::CoordinateDistance = self
            .projected_distance(other, Orientation2D::Vertical)
            .into();

        d_horiz * d_horiz + d_vert * d_vert
    }

    /// Euclidean distance, i.e. 2-norm of the vector from the point to the other.
    fn euclidian_distance(&self, other: &Self) -> C::CoordinateDistance {
        self.distance_squared(other).sqrt()
    }
}

/// A polygon edge.
pub trait Segment<C: CoordinateConcept> {
    /// Type used to represent the end points of the segment.
    type Point: PointConcept<C>;

    /// Get the start (LOW) or end (HIGH) point of the segment.
    fn get_point(&self, dir: Direction1D) -> Self::Point;

    /// Shortcut to get the 'low' end of the segment.
    fn start(&self) -> Self::Point {
        self.get_point(Direction1D::Low)
    }

    /// Shortcut to get the 'high' end of the segment.
    fn end(&self) -> Self::Point {
        self.get_point(Direction1D::High)
    }
}

/// Convert a polygon into an iterator over polygon segments.
pub trait IntoSegments<C: CoordinateConcept> {
    /// Type which represents the segments.
    type Segment: Segment<C>;
    /// Iterator over segments.
    type SegmentIter: Iterator<Item = Self::Segment>;

    /// Iterate over segments/edges of a polygon.
    fn into_segments(self) -> Self::SegmentIter;
}

/// Loop over points/vertices.
pub trait IntoPoints<C: CoordinateConcept> {
    /// Type of the points.
    type Point: PointConcept<C>;
    /// Iterator over points.
    type PointIter: Iterator<Item = Self::Point>;

    /// Iterate over points.
    fn into_points(self) -> Self::PointIter;
}

/// Describe an interval of coordinates by a start value and an end value.
pub trait Interval<Coord>
where
    Coord: Copy,
{
    /// Get the low or high end.
    fn get(&self, d: Direction1D) -> Coord;

    /// Get the low end.
    fn start(&self) -> Coord {
        self.get(Direction1D::Low)
    }

    /// Get the high end.
    fn end(&self) -> Coord {
        self.get(Direction1D::High)
    }
}

/// Concept of an axis-aligned rectangle.
pub trait Rectangle<C: CoordinateConcept>: Polygon90<C> {
    /// Type used for representing a one-dimensional interval.
    type Interval: Interval<C::Coord>;

    /// Get the interval which is spanned by the rectangle in the given orientation.
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

/// Concept of a polygon with axis-aligned edges. The polygon consists of a single closed loop of vertices, i.e. has no holes.
pub trait Polygon90<C: CoordinateConcept>: Polygon<C> + Polygon90WithHoles<C> {
    /// Iterator over alternating x/y coordinates of the points. Starts with an x coordinate.
    type CompactIterator: Iterator<Item = C::Coord>;

    /// Iterate over alternating x/y coordinates of the polygon vertices. Start with an x coordinate.
    fn compact_iter(&self) -> Self::CompactIterator;

    // /// Set the points from a compact iterator.
    // fn set_from_compact(&mut self, iter: Self::CompactIterator) {
    //     self.set(iter);
    //     todo!();
    // }
}

/// Concept of a polygon which consists of a single closed loop of vertices.
pub trait Polygon<C: CoordinateConcept>: PolygonWithHoles<C> + IntoPoints<C> {
    /// Set points from an iterator.
    fn set(&mut self, iter: impl Iterator<Item = <Self as PolygonSet<C>>::Point>);
}

/// Concept of a polygon with axis-aligned edges and holes.
pub trait Polygon90WithHoles<C: CoordinateConcept>: PolygonWithHoles<C> + Polygon90Set<C> {}

/// Concept of a polygon with holes.
pub trait PolygonWithHoles<C: CoordinateConcept>: PolygonSet<C> {
    /// Get the number of holes.
    fn num_holes(&self) -> usize;
}

/// Set of multiple polygons with axis-aligned edges.
pub trait Polygon90Set<C: CoordinateConcept>: PolygonSet<C> {}

/// Set of multiple polygons edges.
pub trait PolygonSet<C: CoordinateConcept>: IntoSegments<C> {
    /// Point type used for the vertices.
    type Point: PointConcept<C>;
    /// Type used for the polygon segments.
    type Segment: Segment<C, Point = Self::Point>;
    // type Polygon: PolygonWithHoles;
    // type PolygonIter: Iterator<Item=Self::Polygon>;
    /// Iterator over all points.
    type AllPoints: Iterator<Item = Self::Point>;

    /// Get number of polygons.
    fn num_polygons(&self) -> usize;

    /// Add the point `p` to all vertices.
    fn convolved(self, p: &Self::Point) -> Self;

    /// Add the point `p` to all vertices.
    fn convolve(&mut self, p: &Self::Point);

    /// Multiply all vertices with a factor.
    fn scaled(self, scale: C::Coord) -> Self;

    /// Multiply all vertices with a factor.
    fn scale(&mut self, scale: C::Coord);

    /// Iterate over all vertices.
    fn all_points(&self) -> Self::AllPoints;
}
