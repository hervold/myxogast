use std::fmt;
use std::fmt::Debug;
use std::ops::Add;
use std::ops::RangeTo;
use std::ops::{Index, IndexMut};

pub type Mmer = u8;
pub const HYPHEN: Mmer = 4;

#[derive(Clone, RustcDecodable, RustcEncodable)]
pub struct Sequence(pub Vec<Mmer>);

pub fn base_to_char(ch: Mmer) -> char {
    match ch {
        0 => 'A',
        1 => 'T',
        2 => 'G',
        3 => 'C',
        HYPHEN => '-',
        _ => 'X',
    }
}

impl Sequence {
    pub fn from_str(seq: &str) -> Result<Sequence, String> {
        let x = seq
            .to_uppercase()
            .chars()
            .map(|ch| match ch {
                'A' => Ok(0),
                'T' => Ok(1),
                'G' => Ok(2),
                'C' => Ok(3),
                '-' => Ok(HYPHEN),
                _ => Err(format!("unrecognized base: {}", ch)),
            })
            .collect();
        match x {
            Ok(arr) => Ok(Sequence(arr)),
            Err(msg) => Err(msg),
        }
    }
    pub fn mmer_to_str(v: &[Mmer]) -> String {
        let viz = |ch: &Mmer| -> char { base_to_char(*ch) };

        v.iter().map(viz).collect()
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }
    pub fn reverse(&self) -> Sequence {
        let mut x: Vec<Mmer> = self.0.clone();
        x.reverse();
        Sequence(x)
    }
}

impl Add for Sequence {
    type Output = Sequence;

    fn add(self, rhs: Sequence) -> Sequence {
        let mut a = self.0.clone();
        a.extend(rhs.0.iter().cloned());
        Sequence(a)
    }
}

impl fmt::Display for Sequence {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", Sequence::mmer_to_str(&self.0[..]))
    }
}

impl fmt::Debug for Sequence {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", Sequence::mmer_to_str(&self.0[..]))
    }
}

impl Index<i32> for Sequence {
    type Output = Mmer;
    fn index<'a>(&'a self, _index: i32) -> &'a Mmer {
        let a: usize = if _index < 0 {
            self.0.len() - (_index.abs() as usize)
        } else {
            _index as usize
        };
        &self.0[a]
    }
}

impl PartialEq for Sequence {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}
impl Eq for Sequence {}
