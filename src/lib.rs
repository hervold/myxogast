extern crate rustc_serialize;
use rustc_serialize::json;

use std::ops::{Index,IndexMut};
use std::fmt;
use std::fmt::Debug;

pub type Mmer = u8;

#[derive(Debug,RustcDecodable,RustcEncodable)]
pub struct Sequence( Vec<Mmer> );

impl Sequence {
    pub fn new(sequence: Vec<Mmer>) -> Sequence {
        Sequence(sequence)
    }

    pub fn from_str( seq : &str ) -> Result<Sequence, String> {
        let x = seq.to_uppercase()
            .chars()
            .map( |ch| match ch {
                'A' => Ok(0),
                'T' => Ok(1),
                'G' => Ok(2),
                'C' => Ok(3),
                _   => Err(format!("unrecognized base: {}", ch))
            } )
            .collect();
        match x {
            Ok(arr) => Ok( Sequence(arr) ),
            Err(msg) => Err(msg)
        }
    }
    pub fn len(&self) -> usize { self.0.len() }
}

impl Index<i32> for Sequence {
    type Output = Mmer;
    fn index<'a>(&'a self, _index: i32) -> &'a Mmer {
        let a : usize = if _index < 0 {self.0.len() - (_index.abs() as usize)} else {_index as usize};
        return &self.0[ a ];
    }
}

#[derive(RustcDecodable,RustcEncodable)]
pub struct Matrix<T:Clone> {
    pub width:   usize,
    pub height:  usize,
    pub data :   Vec<T>,
}

impl<T:Clone> Matrix<T> {
    pub fn new( init: T, w: usize, h: usize ) -> Matrix<T> {
        Matrix {
            width:  w,
            height: h,
            data : {
                // FIXME: why doesn't the vec! macro work here?  this seems inefficient
                let mut v : Vec<T> = Vec::new();
                for _ in 0..(w * h) {
                    v.push(init.clone());
                }
                v
            }
        }
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

impl<T:Clone> Index<(i32, i32)> for Matrix<T> {
    type Output = T;
    fn index<'a>(&'a self, _index: (i32, i32)) -> &'a T {
        let (a, b) = _index;
        let a2 : usize = if a < 0 {self.width - (a.abs() as usize)} else {a as usize};
        let b2 : usize = if b < 0 {self.height - (b.abs() as usize)} else {b as usize};
        return &self.data[ b2 * self.width + a2 ];
    }
}

impl<T:Clone> IndexMut<(i32, i32)> for Matrix<T> {
    fn index_mut<'a>(&'a mut self, _index: (i32, i32)) -> &'a mut T {
        let (a, b) = _index;
        let a2 : usize = if a < 0 {self.width - (a.abs() as usize)} else {a as usize};
        let b2 : usize = if b < 0 {self.height - (b.abs() as usize)} else {b as usize};
        return &mut self.data[ b2 * self.width + a2 ];
    }
}


impl<T:Clone + Debug> fmt::Debug for Matrix<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut s : String = String::new();
        for i in 0 .. self.height {
            let x = &self.data[ (i * self.width) .. ((i+1) * self.width) ];
            s.push_str( &format!("{:?}\n", x));
        }
        write!(f, "{}", s)
    }
}





// used to store HMMer-style probabilities; dimensions: BASE-COUNT x REF-LENGTH
pub type ProbMatr = Matrix<f32>;


#[derive(Debug,RustcDecodable,RustcEncodable)]
pub enum SeqNode {
    Frag { id: u32, val: Sequence },
    Splat,
    Dist { id: u32, scores: ProbMatr },
    List { id: u32, members: Vec<SeqNode> },
    Branch { id: u32, members: Vec<SeqNode> }
}


pub struct AlnParams {
    pub llocal:     bool,        // global or local on the left?
    pub rlocal:     bool,        // global or local on the right?
    pub max_indel:  Option<u8>,  // maximum number of indels before short-circuiting
    pub gap_open:   i16,         // penalty for opening a gap
    pub gap_ext:    i16,         // gap extention penalty
    pub mismatch:   i16,         // mismatch penalty
    pub equal:      i16,         // match score, but "match" is a keyword
}

impl AlnParams {
    pub fn new(  _llocal: Option<bool>, _rlocal: Option<bool>, _max_indel: Option<u8>,
                 _gap_open: Option<i16>, _gap_ext: Option<i16>, _mismatch: i16, _equal: Option<i16> ) -> AlnParams {
        AlnParams {
            llocal:     match _llocal { Some(b) => b, None => true },
            rlocal:     match _rlocal { Some(b) => b, None => true },
            max_indel:  _max_indel,
            gap_open:   match _gap_open { Some(g) => g, None => -1 },
            gap_ext:    match _gap_ext { Some(g) => g, None => -1 },
            mismatch:   _mismatch,
            equal:      match _equal { Some(g) => g, None => 1 }, // "match" is a keyword
        }
    }
}

#[derive(Debug,Clone)]
pub enum AlnState {
    Nil,
    Match,
    Mismatch,
    Ins,
    Del
}

// NOTE: I'm a little shocked I can't just #derive this...
impl PartialEq for AlnState {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (&AlnState::Nil, &AlnState::Nil) => true,
            (&AlnState::Match, &AlnState::Match) => true,
            (&AlnState::Mismatch, &AlnState::Mismatch) => true,
            (&AlnState::Ins, &AlnState::Ins) => true,
            (&AlnState::Del, &AlnState::Del) => true,
            _ => false }
    }
}
impl Eq for AlnState {}

#[derive(Clone)]
pub struct Cell(pub i32);

impl Cell {
    // bit-packing operations
    // FIXME: i16 might not be enough, so I should investigate a different division, eg
    //        4 bits for AlnState and 28 bits get munged into an i32
    pub fn pack( state : AlnState, score : i16 ) -> Cell {
        let _state : i32 = match state {
            AlnState::Nil =>  0xffff,
            AlnState::Match => (1 << 16) | 0xffff,
            AlnState::Mismatch => (2 << 16) | 0xffff,
            AlnState::Ins => (3 << 16) | 0xffff,
            AlnState::Del => (4 << 16) | 0xffff };
        let _score : i32 = (score as i32) | 0xffff0000;
        Cell(_score & _state)
    }

    pub fn unpack( c : Cell ) -> Result<(AlnState, i16), ()> {
        let Cell(packed) = c;
        let state = try!(
            match packed >> 16 {
                0 => Ok(AlnState::Nil),
                1 => Ok(AlnState::Match),
                2 => Ok(AlnState::Mismatch),
                3 => Ok(AlnState::Ins),
                4 => Ok(AlnState::Del),
                _ => Err(()) }
        );
        Ok( (state, packed as i16) )
    }
}

impl fmt::Debug for Cell {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let (a, b) = Cell::unpack( self.clone() ).unwrap();
        write!(f, "{:?}({})", a, b)
    }
}



/// Dynamic-programming alignment
/// params can specify a max_indel, in which case this can return None if a
///   solution can't be found with fewer indels
///
pub fn align( reference: Sequence, query: Sequence, params: AlnParams ) -> Option<Matrix<Cell>> {
    let mut m = Matrix::<Cell>::new( Cell(0), reference.len() + 1, query.len() + 1 );
    let ref_len : i32 = reference.len() as i32;
    let query_len : i32 = query.len() as i32;
    let mut trail : Vec<AlnState> = vec![];

    // initialize edges
    if !params.llocal {
        for i in 0 .. ref_len + 1 { m[ (i as i32, 0) ] = Cell::pack( AlnState::Nil, -i as i16); }
    }
    for j in 0i16 .. (query_len + 1) as i16 { m[ (0, j as i32) ] = Cell::pack( AlnState::Nil, -j as i16); }

    println!("-- initialized:\n");
    println!("{:?}\n", Matrix { width: m.width, height: m.height,
                                data: m.data.iter().map( |c| Cell::unpack(c.clone()).unwrap().1 ).collect() });

    for i in 1 .. ref_len {
        for j in 1 .. query_len {

            let (dstate, del) = Cell::unpack( m[ (i-1, j) ].clone() ).unwrap();
            let del_score = del + match (params.rlocal && i == ref_len, dstate) {
                (true, _) => 0,
                (_, AlnState::Del) => params.gap_ext,
                _ => params.gap_open
            };

            let (istate, ins) = Cell::unpack( m[ (i, j-1) ].clone() ).unwrap();
            let ins_score = ins + if istate == AlnState::Ins { params.gap_ext } else { params.gap_open };

            let (_, diag) = Cell::unpack( m[ (i-1, j-1) ].clone() ).unwrap();
            let diag_score = diag + if reference[i-1] == query[j-1] { params.equal } else { params.mismatch };

            println!("@@@ diag={}, diag_score+{}", diag, if reference[i-1] == query[j-1] { params.equal } else { params.mismatch });

            let (a,b) = {
                if diag_score >= del_score && diag_score >= ins_score {
                    (if reference[i-1] == query[j-1] {AlnState::Match} else {AlnState::Mismatch}, diag_score)
                } else if del_score > diag_score && del_score >= ins_score {
                    (AlnState::Del, del_score)
                } else {
                    (AlnState::Ins, ins_score)
                }};
            println!("@@ [{},{}]: del={}, ins={}, diag={} -> {}", i, j, del_score, ins_score, diag_score, b);

            m[ (i, j) ] = Cell::pack( a.clone(), b );
            trail.push( a.clone() );
        }
    };
    println!("trail: {:?}", trail);
    Some(m)
}

pub fn align_hmm( reference: ProbMatr, query: Sequence, params: AlnParams ) -> Option<()> {
    //let mut m = Matrix::<f32>::new( 0., reference.width + 1, query.len() + 1 );
    None
}

#[test]
fn it_has_unit_tests() {
    assert!(true);
}

#[test]
fn it_correctly_calculates_sequence_length() {
    let s = Sequence::new( vec![0,1,2,3] );
    assert!(s.len() == 4);

    let s = Sequence::from_str("AAAAATGCTCGAAAAAAAA").unwrap();
    assert!(s.len() == 19);
}

#[test]
fn it_rejects_bad_sequences() {
    let bad_seq = Sequence::from_str("ACZXCV");
    match bad_seq {
        Ok(v) => panic!("Invalid characters accepted as a sequence"),
        Err(e) => ()
    }
}
