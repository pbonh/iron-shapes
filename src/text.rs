// Copyright (c) 2018-2020 Thomas Kramer.
// SPDX-FileCopyrightText: 2018-2022 Thomas Kramer
//
// SPDX-License-Identifier: AGPL-3.0-or-later

//! `Text` is used as labels associated with a point.

use crate::point::Point;
use crate::traits::*;
use std::ops::Deref;
use crate::rect::Rect;
use num_traits::NumCast;
use std::fmt;
use std::hash::Hash;

/// Trait for types that can be used as the text of this label.
/// The most simple solution is to use `String`. However, in many cases
/// where the same text is used in many labels it might be desirable to use 'string interning'
/// for more efficient memory usage. Then an `Rc<String>` could be used for instance.
pub trait TextType: Eq + Hash + Clone + fmt::Debug {}

impl<T: Eq + Hash + Clone + fmt::Debug> TextType for T {}

/// A text is a point associated with a string.
/// This struct does not define how the text should be rendered on screen.
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Text<T, S = String> {
    /// Location of the label.
    location: Point<T>,
    /// Text content.
    text: S,
}

/// Display format of the text label.
impl<T, S> fmt::Display for Text<T, S>
    where T: fmt::Display,
          S: TextType + fmt::Display {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Text({}, {})", self.text, self.location)
    }
}

impl<T, S> Deref for Text<T, S>
    where S: Deref<Target=String> {
    type Target = String;

    /// Dereference to String.
    fn deref(&self) -> &Self::Target {
        self.text.deref()
    }
}

impl<T, S: TextType> Text<T, S> {
    /// Create a new text object.
    pub fn new(text: S, location: Point<T>) -> Self {
        Text {
            location,
            text,
        }
    }

    /// Get a reference to the text string.
    pub fn text(&self) -> &S {
        &self.text
    }
}

impl<T: Copy, S: TextType> Text<T, S> {
    /// Get location of the text label.
    #[inline]
    pub fn location(&self) -> Point<T> {
        self.location
    }

    /// Get x-coordinate of the label location.
    #[inline]
    pub fn x(&self) -> T {
        self.location.x
    }

    /// Get y-coordinate of the label location.
    #[inline]
    pub fn y(&self) -> T {
        self.location.y
    }
}


impl<T: Copy + PartialOrd, S> TryBoundingBox<T> for Text<T, S> {
    fn try_bounding_box(&self) -> Option<Rect<T>> {
        Some(
            Rect::new(self.location, self.location)
        )
    }
}


impl<T, Dst, S> TryCastCoord<T, Dst> for Text<T, S>
    where T: Copy + NumCast,
          Dst: Copy + NumCast,
          S: Clone {
    type Output = Text<Dst, S>;

    fn try_cast(&self) -> Option<Self::Output> {
        self.location.try_cast()
            .map(|loc| Text {
                location: loc,
                text: self.text.clone(),
            })
    }
}

/// Point wise transformation for a single point.
impl<T, S> MapPointwise<T> for Text<T, S>
    where T: Copy,
          S: Clone
{
    /// Point wise transformation.
    #[inline]
    fn transform<F>(&self, transformation: F) -> Self
        where F: Fn(Point<T>) -> Point<T> {
        Text {
            location: transformation(self.location),
            text: self.text.clone(),
        }
    }
}