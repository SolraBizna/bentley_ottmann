//! This crate implements the [Bentley­Ottmann algorithm] for finding all
//! intersection points within a set of line segments. See the documentation
//! for [`find_intersections`].
//!
//! [Bentley­Ottmann algorithm]: https://en.wikipedia.org/wiki/Bentley%E2%80%93Ottmann_algorithm
//! [`find_intersections`]: fn.find_intersections.html

// The following references were the main ones that were used:
// - http://geomalgorithms.com/a09-_intersect-3.html
// - https://en.wikipedia.org/wiki/Bentley%E2%80%93Ottmann_algorithm

extern crate num_rational;

use std::fmt::Debug;

mod pqueue;
use pqueue::PriorityQueue;

/// An input coordinate, implemented as an `i32`. Note: If your input
/// coordinates get too large, they may still overflow.
pub type Fixed = i32;
/// An exact intersection coordinate, implemented as a ratio of two `i64`s.
pub type Frac = num_rational::Ratio<i64>;

fn twice_triangle_area(a: (Fixed, Fixed), b: (Fixed, Fixed), c: (Fixed, Fixed))
                       -> Fixed
{
    (b.0 - a.0)*(c.1 - a.1) - (c.0 - a.0)*(b.1 - a.1)
}

/// Returns which side of line a→b that c is on. -1 = left, 1 = right, 0 =
/// coincident.
fn side_of_line(a: (Fixed, Fixed), b: (Fixed, Fixed), c: (Fixed, Fixed))
    -> i32 {
    let area = twice_triangle_area(a, b, c);
    if area < 0 { -1 }
    else if area > 0 { 1 }
    else { 0 }
}

/// Any struct that implements this trait can be used as an input line segment.
/// This trait is automatically implemented by the 4-tuple of [`Fixed`]s, and
/// references thereto.
///
/// [`Fixed`]: type.Fixed.html
pub trait InputLineSegment : Copy + Eq {
    /// Return the exact coordinates of the endpoints of this line segment:
    /// `(x1, y1, x2, y2)`.
    fn get_coords(&self) -> (Fixed, Fixed, Fixed, Fixed);
    /// Return true if this line segment intersects with another. You only need
    /// to implement this if your internal line segment representation has a
    /// natural way to quickly check intersection. This method must be exact!
    ///
    /// The default implementation assumes that two lines that share an
    /// endpoint are consecutive parts of the same "path" and should not
    /// intersect. If your implementation has another method of preventing
    /// these spurious intersections, it need not check endpoint coordinates.
    fn has_intersection(&self, other: &Self) -> bool {
        let (a1x, a1y, a2x, a2y) = self.get_coords();
        let (b1x, b1y, b2x, b2y) = other.get_coords();
        if (a1x == b1x && a1y == b1y) || (a1x == b2x && a1y == b2y)
        || (a2x == b1x && a2y == b1y) || (a2x == b2x && a2y == b2y) {
            // treating any shared endpoint as non-intersection
            // sure hope this doesn't bite us in the wild...
            // (if it does, we can fix this with some careful preprocessing)
            return false
        }
        let b1side = side_of_line((a1x, a1y), (a2x, a2y), (b1x, b1y));
        let b2side = side_of_line((a1x, a1y), (a2x, a2y), (b2x, b2y));
        if b1side == b2side {
            // b1 and b2 are on the same side of a1→a2, therefore intersection
            // is not possible
            return false
        }
        let a1side = side_of_line((b1x, b1y), (b2x, b2y), (a1x, a1y));
        let a2side = side_of_line((b1x, b1y), (b2x, b2y), (a2x, a2y));
        if a1side == a2side {
            // a1 and a2 are on the same side of b1→b2, therefore intersection
            // is not possible
            return false
        }
        true
    }
    /// Returns the exact intersection point of two line segments, assuming it
    /// exists. You only need to implement this if your internal line segment
    /// representation has a natural way to quickly calculate *exact*
    /// intersection points. The resulting coordinates must be exact!
    fn get_intersection(&self, other: &Self)
                        -> (Frac, Frac) {
        let (a1x, a1y, a2x, a2y) = self.get_coords();
        let (b1x, b1y, b2x, b2y) = other.get_coords();
        let (a1x, a1y) = (Frac::from(a1x as i64), Frac::from(a1y as i64));
        let (a2x, a2y) = (Frac::from(a2x as i64), Frac::from(a2y as i64));
        let (b1x, b1y) = (Frac::from(b1x as i64), Frac::from(b1y as i64));
        let (b2x, b2y) = (Frac::from(b2x as i64), Frac::from(b2y as i64));
        // the normal of a1→a2
        let anx = a2y - a1y;
        let any = a1x - a2x; // = -(a2x-a1x)
        let b1dot = anx * (b1x - a1x) + any * (b1y - a1y);
        let b2dot = anx * (b2x - a1x) + any * (b2y - a1y);
        debug_assert!((*b1dot.numer() < 0) != (*b2dot.numer() < 0));
        let b1dot = if *b1dot.numer() < 0 { -b1dot } else { b1dot };
        let b2dot = if *b2dot.numer() < 0 { -b2dot } else { b2dot };
        let totaldot = b1dot + b2dot;
        (b1x + (b2x-b1x) * b1dot / totaldot,
         b1y + (b2y-b1y) * b1dot / totaldot)
    }
    /// Returns `None` if the two lines do not intersect, or `Some` exact
    /// intersection point if they do. You only need to implement this method
    /// if your internal line segment type has a natural way of simultaneously
    /// determining the existence of and exact location of an intersection
    /// point; otherwise, will be implemented in terms of [`has_intersection`]
    /// and [`get_intersection`]. This method's calculations must be exact!
    ///
    /// [`has_intersection`]: #method.has_intersection
    /// [`get_intersection`]: #method.get_intersection
    fn maybe_get_intersection(&self, other: &Self)
                              -> Option<(Frac, Frac)> {
        if self.has_intersection(other) { Some(self.get_intersection(other)) }
        else { None }
    }
}

impl InputLineSegment for (Fixed, Fixed, Fixed, Fixed) {
    fn get_coords(&self) -> (Fixed, Fixed, Fixed, Fixed) { *self }
}

impl<'a> InputLineSegment for &'a(Fixed, Fixed, Fixed, Fixed) {
    fn get_coords(&self) -> (Fixed, Fixed, Fixed, Fixed) { **self }
}

#[derive(Debug,PartialEq,Eq)]
enum CandidateEvent<LineSegmentType: Copy + PartialEq + Debug> {
    LeftEndpoint(LineSegmentType),
    RightEndpoint(LineSegmentType),
    // part of the implementation assumes that the first segment is *below* the
    // second
    Intersection(LineSegmentType, LineSegmentType),
}
use CandidateEvent::*;

#[derive(Debug)]
struct ActiveLine<LineSegmentType: Copy + Debug> {
    cur_y: Frac,
    slope: Frac,
    seg: LineSegmentType,
}

/// Uses the Bentley-Ottmann algorithm to find all intersections of the given
/// line segments.
///
/// Finds the beginning (left) and end (right) of each line, and position of
/// each intersection, in a deterministic order: left-to-right,
/// low-Y-to-high-Y. For vertical lines, the "left" endpoint is the one with
/// the lower Y coordinate.
///
/// Calls the respective handler for each thing found.
///
/// The default implementation of [`InputLineSegment`] assumes that any two
/// lines which share an endpoint are consecutive parts of the same "path" and
/// therefore should not intersect. If your input consists of multiple "paths",
/// you should preprocess the input so that no line segments from different
/// paths share an endpoint. (Or, implement [`InputLineSegment`] with knowledge
/// of path membership of segments.)
///
/// [`InputLineSegment`]: trait.InputLineSegment.html
pub fn find_intersections<LineIterator, LineSegmentType,
                          IntersectionHandler, LineBeginHandler,
                          LineEndHandler>
(lines: LineIterator,
 mut intersection: IntersectionHandler,
 mut line_begin: LineBeginHandler,
 mut line_end: LineEndHandler) -> ()
where LineSegmentType: InputLineSegment + PartialEq + Debug,
      LineIterator: Iterator<Item=LineSegmentType>,
      IntersectionHandler: FnMut(LineSegmentType, LineSegmentType, Frac, Frac),
      LineBeginHandler: FnMut(LineSegmentType, Fixed, Fixed),
      LineEndHandler: FnMut(LineSegmentType, Fixed, Fixed),
{
    let mut queue = PriorityQueue::new();
    for line in lines {
        let (x1, y1, x2, y2) = line.get_coords();
        if x1 > x2 || (x1 == x2 && y1 > y2) {
            queue.insert((Frac::from(x2 as i64), Frac::from(y2 as i64)),
                          LeftEndpoint(line));
            queue.insert((Frac::from(x1 as i64), Frac::from(y1 as i64)),
                          RightEndpoint(line));
        }
        else {
            queue.insert((Frac::from(x1 as i64), Frac::from(y1 as i64)),
                          LeftEndpoint(line));
            queue.insert((Frac::from(x2 as i64), Frac::from(y2 as i64)),
                          RightEndpoint(line));
        }
    }
    // we're using a Vec instead of a binary tree, on the grounds that this
    // vector will usually not contain more than a few elements, and that we'll
    // get better performance in practice by linear-searching and sorting a
    // small array than a small tree
    let mut active: Vec<ActiveLine<LineSegmentType>> = Vec::new();
    while let Some((key, event)) = queue.pop() {
        let x = key.0;
        let y = key.1;
        for line in active.iter_mut() {
            let (x1, y1, x2, _y2) = line.seg.get_coords();
            if x1 != x2 {
                line.cur_y = Frac::from(y1 as i64) + (x -Frac::from(x1 as i64))
                    * line.slope;
            }
        }
        active.sort_by_key(|x| x.cur_y);
        match event {
            LeftEndpoint(seg) => {
                line_begin(seg,
                           x.to_integer() as Fixed,
                           y.to_integer() as Fixed);
                let (x1, y1, x2, y2) = seg.get_coords();
                let slope = if x2 == x1 { Frac::from(0) }
                else { Frac::new((y2 - y1) as i64, (x2 - x1) as i64) };
                let new_line = ActiveLine {
                    cur_y: y,
                    slope,
                    seg: seg,
                };
                let mut target_index = active.len();
                for n in 0 .. active.len() {
                    if active[n].cur_y > new_line.cur_y {
                        target_index = n;
                        break;
                    }
                }
                if target_index != 0 && target_index != active.len() {
                    let line_below = &active[target_index-1];
                    let line_above = &active[target_index];
                    if let Some(intersection)
                    = line_below.seg.maybe_get_intersection(&line_above.seg) {
                        // if that intersection is a candidate event, remove it
                        let ev = Intersection(line_below.seg, line_above.seg);
                        queue.remove(&intersection, |x| x == &ev);
                    }
                }
                if target_index != 0 {
                    let line_below = &active[target_index-1];
                    if let Some(intersection)
                    = line_below.seg.maybe_get_intersection(&seg) {
                        queue.insert(intersection,
                                     Intersection(line_below.seg, seg));
                    }
                }
                if target_index != active.len() {
                    let line_above = &active[target_index];
                    if let Some(intersection)
                    = seg.maybe_get_intersection(&line_above.seg) {
                        queue.insert(intersection,
                                     Intersection(seg, line_above.seg));
                    }
                }
                active.insert(target_index, new_line);
            },
            RightEndpoint(seg) => {
                line_end(seg,
                         x.to_integer() as Fixed,
                         y.to_integer() as Fixed);
                let mut index = None;
                for n in 0 .. active.len() {
                    if active[n].seg == seg { index = Some(n); break }
                }
                if let Some(index) = index {
                    active.remove(index);
                    if index != 0 && index != active.len() {
                        let line_below = &active[index-1];
                        let line_above = &active[index];
                        if let Some(intersection)
                        = line_below.seg.maybe_get_intersection(&line_above.seg) {
                            queue.insert(intersection,
                                         Intersection(line_below.seg,
                                                      line_above.seg));
                        }
                    }
                }
            },
            Intersection(lower_seg, higher_seg) => {
                intersection(lower_seg, higher_seg, x, y);
                let mut lower_index = None;
                for n in 0 .. active.len() {
                    if active[n].seg == lower_seg {
                        lower_index = Some(n);
                        break
                    }
                }
                let lower_index = lower_index.expect("segment slipped out!");
                assert!(lower_index + 1 < active.len());
                let higher_index = lower_index + 1;
                if lower_index != 0 {
                    let line_below = &active[lower_index-1];
                    if let Some(intersection)
                    = line_below.seg.maybe_get_intersection(&lower_seg) {
                        // if that intersection is a candidate event, remove it
                        let ev = Intersection(line_below.seg, lower_seg);
                        queue.remove(&intersection, |x| x == &ev);
                    }
                    if let Some(intersection)
                    = line_below.seg.maybe_get_intersection(&higher_seg) {
                        queue.insert(intersection,
                                     Intersection(line_below.seg,
                                                  higher_seg));
                    }
                }
                if higher_index+1 != active.len() {
                    let line_above = &active[higher_index+1];
                    if let Some(intersection)
                    = higher_seg.maybe_get_intersection(&line_above.seg) {
                        // if that intersection is a candidate event, remove it
                        let ev = Intersection(higher_seg, line_above.seg);
                        queue.remove(&intersection, |x| x == &ev);
                    }
                    if let Some(intersection)
                    = lower_seg.maybe_get_intersection(&line_above.seg) {
                        queue.insert(intersection,
                                     Intersection(lower_seg,
                                                  line_above.seg));
                    }
                }
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::RefCell;
    #[test]
    fn test_twice_triangle_area() {
        assert_eq!(twice_triangle_area((10, 10), (20, 10), (10, 20)),
                   100);
        assert_eq!(twice_triangle_area((10, 10), (10, 20), (20, 10)),
                   -100);
        assert_eq!(twice_triangle_area((15, 10), (10, 20), (20, 10)),
                   -50);
    }
    #[test]
    fn some_intersections() {
        let line1 = (20, 10, 20, 30);
        let line2 = (10, 20, 30, 20);
        let line3 = (10, 10, 50, 30);
        let line4 = (10, 10, 30, 10);
        let line5 = (10, 9, 30, 9);
        assert_eq!(line1.maybe_get_intersection(&line2),
                   Some((20.into(), 20.into())));
        assert_eq!(line1.maybe_get_intersection(&line3),
                   Some((20.into(), 15.into())));
        assert_eq!(line1.maybe_get_intersection(&line4),
                   Some((20.into(), 10.into())));
        assert_eq!(line1.maybe_get_intersection(&line5),
                   None);
    }
    type TestSeg = (Fixed, Fixed, Fixed, Fixed);
    #[derive(Debug,PartialEq)]
    enum Event {
        Begin(TestSeg, Fixed, Fixed),
        End(TestSeg, Fixed, Fixed),
        Intersection(TestSeg, TestSeg, Frac, Frac),
    }
    use self::Event::{Begin,End,Intersection};
    #[test]
    fn pentagram() {
        let lines: &[TestSeg] = &[
            (-3, -5,  0,  5),
            ( 0,  5,  3, -5),
            ( 3, -5, -5,  1),
            (-5,  1,  5,  1),
            ( 5,  1, -3, -5),
        ];
        let events = RefCell::new(Vec::new());
        find_intersections(lines.iter(),
                           |a,b,x,y| events.borrow_mut()
                           .push(Event::Intersection(*a, *b, x, y)),
                           |a,x,y| events.borrow_mut()
                           .push(Event::Begin(*a, x, y)),
                           |a,x,y| events.borrow_mut()
                           .push(Event::End(*a, x, y)));
        assert_eq!(&events.into_inner()[..],
                   &[
                       Begin((-5, 1, 5, 1), -5, 1),
                       Begin((3, -5, -5, 1), -5, 1),
                       Begin((5, 1, -3, -5), -3, -5),
                       Begin((-3, -5, 0, 5), -3, -5),
                       Intersection((-3, -5, 0, 5), (3, -5, -5, 1),
                                    Frac::new(-93, 49), Frac::new(-65, 49)),
                       Intersection((-3, -5, 0, 5), (-5, 1, 5, 1),
                                    Frac::new(-6, 5), Frac::new(1, 1)),
                       Intersection((5, 1, -3, -5), (3, -5, -5, 1),
                                    Frac::new(0, 1), Frac::new(-11, 4)),
                       Begin((0, 5, 3, -5), 0, 5),
                       End((-3, -5, 0, 5), 0, 5),
                       Intersection((-5, 1, 5, 1), (0, 5, 3, -5),
                                    Frac::new(6, 5), Frac::new(1, 1)),
                       Intersection((5, 1, -3, -5), (0, 5, 3, -5),
                                    Frac::new(93, 49), Frac::new(-65, 49)),
                       End((3, -5, -5, 1), 3, -5),
                       End((0, 5, 3, -5), 3, -5),
                       End((5, 1, -3, -5), 5, 1),
                       End((-5, 1, 5, 1), 5, 1),
                   ]);
    }
}
