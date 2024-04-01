extern crate bio;
extern crate myxogast;
extern crate rustc_serialize;

extern crate argparse;

use bio::io::{fasta, fastq};
use rustc_serialize::json;
use rustc_serialize::json::Json;
use std::clone::Clone;
use std::cmp::{max, Ordering, PartialOrd};
use std::fmt;
use std::fmt::Debug;
use std::fs::File;
use std::io;
use std::io::Read;
use std::ops::{Index, IndexMut};
use std::str;

use myxogast::align::*;
use myxogast::matrix::*;
use myxogast::seq::*;
use myxogast::tree::*;

use argparse::{ArgumentParser, Store, StoreTrue};

// FIXME:
const params: AlnParams = AlnParams {
    llocal: true,
    rlocal: true,
    max_indel: None,
    gap_open: -1,
    gap_ext: -1,
    mismatch: -1,
    equal: 1,
};

fn main() {
    let mut ref_fname: String = String::new();
    let mut query_fname: String = String::new();
    {
        // scope block?
        let mut parser = ArgumentParser::new();
        parser
            .refer(&mut ref_fname)
            .add_argument("ref", Store, "refernce JSON file")
            .required();
        parser
            .refer(&mut query_fname)
            .add_argument("query", Store, "query file")
            .required();
        parser.parse_args_or_exit();
    }

    let mut graph_s = String::new();
    File::open(ref_fname).unwrap().read_to_string(&mut graph_s);
    let graph = SeqGraph::from_json(&graph_s).unwrap();
    println!("graph: {:?}", graph);

    // FIXME: it seems like Rust-Bio implements nearly IDENTICAL structures
    //   for records from FASTA and FASTQ files, but they're not the same types!
    let query_lc = query_fname.to_lowercase();
    let reader = if query_lc.ends_with(".fasta") || query_lc.ends_with(".fa") {
        for gene in fasta::Reader::from_file(query_fname).unwrap().records() {
            let seq = gene.unwrap();
            let gene_name = seq.id();

            // FIXME: ugh ..
            let seq_str: String = String::from_utf8_lossy(&seq.seq()).into_owned();
            println!("name:{}", gene_name);
            let query = Sequence::from_str(seq_str.as_str()).unwrap();
            let (_, tgt, refe) = graph.align__global_max(&query, &params).unwrap();
            println!("name:{} - {}, {}", gene_name, tgt, refe);
        }
    } else if query_lc.ends_with(".fastq") || query_lc.ends_with(".fq") {
        for gene in fastq::Reader::from_file(query_fname).unwrap().records() {
            let seq = gene.unwrap();

            // FIXME: ugh ..
            let seq_str: String = String::from_utf8_lossy(&seq.seq()).into_owned();

            let query = Sequence::from_str(seq_str.as_str()).unwrap();
        }
    } else {
        panic!("don't recognize file type: {}", query_fname);
    };
}
