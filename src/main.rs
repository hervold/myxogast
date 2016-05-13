extern crate rustc_serialize;
extern crate myxogast;

use rustc_serialize::json;
use std::ops::{Index,IndexMut};
use std::clone::Clone;

use myxogast::*;

fn main() {
    let x = SeqNode::Frag { id: 0, val: Sequence::new( vec![0,1,2,3] ) };
    println!("Hello World: {:?}", x );
    println!("{:?}", json::encode(&x));
    let mut y = Matrix::<u32>::new(0, 10, 10);
    y[(1,0)] = 5u32;
    println!("5 == {}?", y[(1,0)]);

    let q = Sequence::from_str("atgc");
    println!("{:?}", q);
    let r = Sequence::from_str("atgcn");
    println!("{:?}", r);

    assert_eq!( Cell::unpack( Cell::pack( AlnState::Ins, -10 ) ).unwrap(), (AlnState::Ins, -10) );
    assert_eq!( Cell::unpack( Cell::pack( AlnState::Match, 222 ) ).unwrap(), (AlnState::Match, 222) );

    let params = AlnParams {
        llocal:    false,
        rlocal:    false,
        max_indel: None,
        gap_open:  -1,
        gap_ext:   -1,
        mismatch:  -1,
        equal:     1 };

    let x = align( Sequence::from_str("AAAAATGCTCGAAAAAAAA").unwrap(),
                   Sequence::from_str("TGCTCG").unwrap(),
                   params ).unwrap();

    /*
    let x = align( Sequence::from_str("TGCCG").unwrap(),
                   Sequence::from_str("TGCTCG").unwrap(),
                   params ).unwrap();
     */
    println!("\n{:?}\n\n", x );

    let x_scores = Matrix { width: x.width, height: x.height, data: x.data.iter().map( |c| Cell::unpack(c.clone()).unwrap().1 ).collect() };
    println!("\n{:?}", x_scores );
}
