use matrix::Matrix;
use seq::*;
use std::cmp::min;
use std::fmt;
use std::fmt::Debug;

pub type AlnScore = i32;

const SCOREBITS: usize = 28;

/// Cell: bit-pack maximal alignment state and score into 32 bits
///
/// layout: 1bit for sign; 3 bits for state (Ins, Del, etc); 28 bits for score
///
/// note that a simple Smith-Waterman doesn't require storing the "direction" of the
/// previous max in the cell, however, we need it to implement separate gap-opening
/// and gap-extension penalties
///
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Cell(pub u32);

impl Cell {
    // bit-packing operations
    pub fn pack(state: &AlnState, score: &AlnScore) -> Cell {
        // 0x1fffffff = 29 1's
        let _state: u32 = if *score < 0 { 0x8fffffff } else { 0xfffffff }
            | (match *state {
                AlnState::Nil => 0,
                AlnState::Match => 1 << SCOREBITS,
                AlnState::Mismatch => 2 << SCOREBITS,
                AlnState::Ins => 3 << SCOREBITS,
                AlnState::Del => 4 << SCOREBITS,
            });

        // ensure the first 4 bits are either 0111 or 1111, depending on sign
        let _score: u32 = (*score).abs() as u32 | {
            if *score < 0 {
                0xf0000000
            } else {
                0x70000000
            }
        };
        Cell(_score & _state)
    }

    pub fn unpack(c: &Cell) -> Result<(AlnState, AlnScore), String> {
        let Cell(packed) = *c;
        let state = try!(match (packed >> SCOREBITS) & 7 {
            0 => Ok(AlnState::Nil),
            1 => Ok(AlnState::Match),
            2 => Ok(AlnState::Mismatch),
            3 => Ok(AlnState::Ins),
            4 => Ok(AlnState::Del),
            x => Err(format!("unknown value: {}", x)),
        });

        // the first bit is a sign bit
        let score: AlnScore = if (packed & 0x80000000) == 0 {
            0xfffffff & packed as i32
        } else {
            // score was negative
            -1 * (0xfffffff & packed as i32)
        };
        Ok((state, score))
    }
}

impl fmt::Debug for Cell {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let (a, b) = Cell::unpack(&self).unwrap();
        write!(f, "{:?}({})", a, b)
    }
}

// used to store HMMer-style probabilities; dimensions: BASE-COUNT x REF-LENGTH
pub type ProbMatr = Matrix<f32>;

#[derive(Clone, Debug)]
pub struct AlnParams {
    pub llocal: bool,          // global or local on the left?
    pub rlocal: bool,          // global or local on the right?
    pub max_indel: Option<u8>, // maximum number of indels before short-circuiting
    pub gap_open: AlnScore,    // penalty for opening a gap
    pub gap_ext: AlnScore,     // gap extention penalty
    pub mismatch: AlnScore,    // mismatch penalty
    pub equal: AlnScore,       // match score, but "match" is a keyword
}

impl AlnParams {
    pub fn new(
        _llocal: Option<bool>,
        _rlocal: Option<bool>,
        _max_indel: Option<u8>,
        _gap_open: Option<AlnScore>,
        _gap_ext: Option<AlnScore>,
        _mismatch: AlnScore,
        _equal: Option<AlnScore>,
    ) -> AlnParams {
        AlnParams {
            llocal: match _llocal {
                Some(b) => b,
                None => true,
            },
            rlocal: match _rlocal {
                Some(b) => b,
                None => true,
            },
            max_indel: _max_indel,
            gap_open: match _gap_open {
                Some(g) => g,
                None => -2,
            },
            gap_ext: match _gap_ext {
                Some(g) => g,
                None => -1,
            },
            mismatch: _mismatch,
            equal: match _equal {
                Some(g) => g,
                None => 1,
            }, // "match" is a keyword
        }
    }
    pub fn copy_but_llocal(tgt: &AlnParams, ll: bool) -> AlnParams {
        AlnParams {
            llocal: ll,
            rlocal: tgt.rlocal,
            max_indel: tgt.max_indel,
            gap_open: tgt.gap_open,
            gap_ext: tgt.gap_ext,
            mismatch: tgt.mismatch,
            equal: tgt.equal,
        }
    }
}

#[derive(Debug, Clone)]
pub enum AlnState {
    Nil,
    Match,
    Mismatch,
    Ins,
    Del,
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
            _ => false,
        }
    }
}
impl Eq for AlnState {}

/// Dynamic-programming alignment
/// params can specify a max_indel, in which case this can return None if a
///   solution can't be found with fewer indels
///
/// returns Option(m, i, j), where m is the alignment matrix, and i,j is the
///   the location of the cell with the highest value
///
pub fn align_matrix(
    reference: &Sequence,
    query: &Sequence,
    params: &AlnParams,
    _ref_offset: Option<i32>,
    m: &mut Matrix<Cell>,
) -> Option<(/*x*/ i32, /*y*/ i32)> {
    let ref_len: i32 = reference.len() as i32;
    let ref_offset = match _ref_offset {
        Some(x) => x,
        None => 0i32,
    };
    let query_len: i32 = query.len() as i32;

    let mut best_val: AlnScore = AlnScore::min_value();
    let mut best_loc: (i32, i32) = (0, 0);

    // initialize edges
    if !params.llocal {
        for i in 0..ref_len + 1 {
            m[(i as i32, 0)] = Cell::pack(&AlnState::Nil, &(-i as AlnScore));
        }
    }
    for j in 0i16..(query_len + 1) as i16 {
        m[(0, j as i32)] = Cell::pack(&AlnState::Nil, &(-j as AlnScore));
    }

    for i in 1..ref_len + 1 {
        for j in 1..query_len + 1 {
            let m_i = ref_offset + i;

            let (dstate, del) = Cell::unpack(&m[(m_i - 1, j)]).unwrap();
            let del_score = del
                + match (params.rlocal && i == ref_len, dstate) {
                    (true, _) => 0,
                    (_, AlnState::Del) => params.gap_ext,
                    _ => params.gap_open,
                };

            let (istate, ins) = Cell::unpack(&m[(m_i, j - 1)]).unwrap();
            let ins_score = ins
                + if istate == AlnState::Ins {
                    params.gap_ext
                } else {
                    params.gap_open
                };

            let (_, diag) = Cell::unpack(&m[(m_i - 1, j - 1)]).unwrap();
            let diag_score = diag
                + if reference[i - 1] == query[j - 1] {
                    params.equal
                } else {
                    params.mismatch
                };

            let (a, b) = {
                if diag_score >= del_score && diag_score >= ins_score {
                    (
                        if reference[i - 1] == query[j - 1] {
                            AlnState::Match
                        } else {
                            AlnState::Mismatch
                        },
                        diag_score,
                    )
                } else if del_score > diag_score && del_score >= ins_score {
                    (AlnState::Del, del_score)
                } else {
                    (AlnState::Ins, ins_score)
                }
            };

            if b > best_val {
                best_val = b;
                best_loc = (m_i, j);
            }
            m[(m_i, j)] = Cell::pack(&a, &b);
        }
    }
    Some((best_loc.0, best_loc.1))
}

pub fn align_hmm(reference: ProbMatr, query: Sequence, params: AlnParams) -> Option<()> {
    //let mut m = Matrix::<f32>::new( 0., reference.width + 1, query.len() + 1 );
    None
}

pub fn aln_from_coord(
    st_i: &i32,
    st_j: &i32,
    _inc: &i32,
    reference: &Sequence,
    query: &Sequence,
    alignment: &Matrix<Cell>,
) -> (Sequence, Sequence) {
    let ref_len: i32 = reference.len() as i32;
    let query_len: i32 = query.len() as i32;

    let mut padded_ref: Vec<Mmer> = Vec::with_capacity(reference.len());
    let mut padded_query: Vec<Mmer> = Vec::with_capacity(query.len());
    let inc = *_inc;
    let mut i = *st_i;
    let mut j = *st_j;

    let max_i = if inc < 0 { ref_len + 1 } else { ref_len };
    let max_j = if inc < 0 { query_len + 1 } else { query_len };

    while i < max_i && i > 1 && j < max_j && j > 1 {
        let (_, m_sc) = Cell::unpack(&alignment[(i + inc, j + inc)]).unwrap();
        let (_, r_sc) = Cell::unpack(&alignment[(i + inc, j)]).unwrap();
        let (_, d_sc) = Cell::unpack(&alignment[(i, j + inc)]).unwrap();

        if m_sc >= r_sc && m_sc >= d_sc {
            i += inc;
            j += inc;
            padded_ref.push(reference[i - 1]);
            padded_query.push(query[j - 1]);
        } else if r_sc >= m_sc && r_sc >= d_sc {
            i += inc;
            padded_ref.push(reference[i - 1]);
            padded_query.push(HYPHEN);
        } else {
            j += inc;
            padded_ref.push(HYPHEN);
            padded_query.push(query[j - 1]);
        }
    }

    while i < ref_len && i > 1 {
        i += inc;
        padded_ref.push(reference[i - 1]);
        padded_query.push(HYPHEN);
    }

    while j < query_len && j > 1 {
        j += inc;
        padded_ref.push(HYPHEN);
        padded_query.push(query[j - 1]);
    }

    (Sequence(padded_ref), Sequence(padded_query))
}

pub fn align(
    reference: &Sequence,
    query: &Sequence,
    params: &AlnParams,
) -> Option<(Sequence, Sequence)> {
    let mut m = Matrix::<Cell>::new(Cell(0), reference.len() + 2, query.len() + 2);

    match align_matrix(reference, query, &params, None, &mut m) {
        Some((st_i, st_j)) => {
            let (fwd_r, fwd_q) = aln_from_coord(&st_i, &st_j, &1, &reference, &query, &m);
            let (rev_r1, rev_q1) = aln_from_coord(&st_i, &st_j, &(-1), &reference, &query, &m);

            let mut rev_r2 = rev_r1.reverse().0;
            if st_i <= reference.len() as i32 {
                rev_r2.push(reference[st_i - 1]);
            }

            let mut rev_q2 = rev_q1.reverse().0;
            if st_j <= query.len() as i32 {
                rev_q2.push(query[st_j - 1]);
            }

            Some((Sequence(rev_r2) + fwd_r, Sequence(rev_q2) + fwd_q))
        }

        None => None,
    }
}
