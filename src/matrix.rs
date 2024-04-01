use std::fmt;
use std::fmt::Debug;
use std::ops::{Index, IndexMut};

#[derive(Clone, RustcDecodable, RustcEncodable)]
pub struct Matrix<T> {
    pub width: usize,
    pub height: usize,
    pub data: Vec<T>,
}

impl<T: Clone + Ord> Matrix<T> {
    pub fn new(init: T, w: usize, h: usize) -> Matrix<T> {
        Matrix {
            width: w,
            height: h,
            data: {
                // FIXME: why doesn't the vec! macro work here?  this seems inefficient
                let mut v: Vec<T> = Vec::new();
                for _ in 0..(w * h) {
                    v.push(init.clone());
                }
                v
            },
        }
    }

    pub fn max(&self) -> (/*val*/ T, (/*x*/ i32, /*y*/ i32)) {
        let idx = self
            .data
            .iter()
            .enumerate()
            .max_by_key(|&(_, item)| item)
            .unwrap()
            .0;
        (
            self.data[idx].clone(),
            ((idx % self.width) as i32, (idx / self.height) as i32),
        )
    }

    /*
    pub fn map<T,F>(&'a self, f: F) -> Matrix<B>
    where F: &T -> B {
        Matrix {
            width:  self.width,
            height: self.height,
            data:   self.data.map(f)
        }
    }
    */
}

impl<T> Index<(i32, i32)> for Matrix<T> {
    type Output = T;
    fn index<'a>(&'a self, _index: (i32, i32)) -> &'a T {
        let (a, b) = _index;
        let a2: usize = if a < 0 {
            self.width - (a.abs() as usize)
        } else {
            a as usize
        };
        let b2: usize = if b < 0 {
            self.height - (b.abs() as usize)
        } else {
            b as usize
        };
        return &self.data[b2 * self.width + a2];
    }
}

impl<T> IndexMut<(i32, i32)> for Matrix<T> {
    fn index_mut<'a>(&'a mut self, _index: (i32, i32)) -> &'a mut T {
        let (a, b) = _index;
        let a2: usize = if a < 0 {
            self.width - (a.abs() as usize)
        } else {
            a as usize
        };
        let b2: usize = if b < 0 {
            self.height - (b.abs() as usize)
        } else {
            b as usize
        };
        return &mut self.data[b2 * self.width + a2];
    }
}

impl<T: Debug> fmt::Debug for Matrix<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut s: String = String::new();
        for i in 0..self.height {
            let x = &self.data[(i * self.width)..((i + 1) * self.width)];
            s.push_str(&format!("{:?}\n", x));
        }
        write!(f, "{}", s)
    }
}
