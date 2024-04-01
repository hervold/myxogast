extern crate rustc_serialize;

pub mod align;
pub mod matrix;
pub mod seq;
pub mod tree;

use align::*;
use matrix::*;
use seq::*;
use tree::*;

const params: AlnParams = AlnParams {
    llocal: false,
    rlocal: false,
    max_indel: None,
    gap_open: -1,
    gap_ext: -1,
    mismatch: -1,
    equal: 1,
};

#[test]
fn test_basic() {
    assert_eq!(
        Cell::unpack(&Cell::pack(&AlnState::Ins, &-10)).unwrap(),
        (AlnState::Ins, -10)
    );
    assert_eq!(
        Cell::unpack(&Cell::pack(&AlnState::Match, &222)).unwrap(),
        (AlnState::Match, 222)
    );

    let s = "ATGCATGC";
    assert_eq!(format!("{}", Sequence::from_str(s).unwrap()), s);
}

#[test]
fn test_aln() {
    let reference = Sequence::from_str("AAAAATGCTCGAAAAAAAA").unwrap();
    let query = Sequence::from_str("TGCTCG").unwrap();

    let (r, q) = align(&reference, &query, &params).unwrap();
    assert_eq!(r, reference);
    assert_eq!(q, Sequence::from_str("-----TGCTCG--------").unwrap());

    let ref2 = Sequence::from_str("ATGCAT").unwrap();
    let query2 = Sequence::from_str("ATGCA").unwrap();
    let (r2, q2) = align(&ref2, &query2, &params).unwrap();
    assert_eq!(r2, ref2);
    assert_eq!(q2, Sequence::from_str("ATGCA-").unwrap());
}

#[test]
fn test_graph() {
    let loc_g = SeqGraph::from_json(
        r#"[{"branch": [["ATCG",{"branch":["TTGG","AAAA"]}],  ["ATGC","TTTT"]]}]"#,
    )
    .unwrap();
    let loc_q = Sequence::from_str("ATGCAAAA").unwrap();

    assert_eq!(
        loc_g.align__global_max(&loc_q, &params),
        Some((
            vec![24, 16, 15, 14],
            Sequence::from_str("A-TGCAAAA").unwrap(),
            Sequence::from_str("ATC-GAAAA").unwrap()
        ))
    );

    assert_eq!(
        loc_g.align__local_max(&loc_q, &params),
        Some((
            vec![24, 23, 22],
            Sequence::from_str("ATGCAAAA").unwrap(),
            Sequence::from_str("ATGCTTTT").unwrap()
        ))
    );
}
