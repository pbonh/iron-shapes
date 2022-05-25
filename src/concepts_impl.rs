use crate::concepts::{PointConcept, IntCoordinates};
use crate::isotropy::Orientation2D;
use crate::prelude::Point;

impl PointConcept<IntCoordinates> for Point<i32> {

    fn new(x: i32, y: i32) -> Self {
        Point::new(x, y)
    }

    fn get(&self, orient: Orientation2D) -> i32 {
        match orient {
            Orientation2D::Horizontal => self.x,
            Orientation2D::Vertical => self.y
        }
    }

    fn set(&mut self, orient: Orientation2D, value: i32) {
        let c = match orient {
            Orientation2D::Horizontal => &mut self.x,
            Orientation2D::Vertical => &mut self.y,
        };
        *c = value
    }
}