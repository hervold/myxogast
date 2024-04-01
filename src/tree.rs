extern crate serde_json;

//use std::fmt::Write;
use self::serde_json::Value as JSON_Val;
use align::*;
use matrix::*;
use seq::*;
use std::cmp::max;
use std::collections::BTreeMap;
use std::io::Write;
use std::iter::Iterator;
use std::rc::Rc;

#[derive(Debug, Clone)]
pub enum SeqNode {
    Nil,

    Frag {
        id: u32,
        val: Sequence,
        llocal: bool,
        rlocal: bool,

        next: Rc<SeqNode>,
    },

    //Dist { id: u32, scores: ProbMatr },
    Branch {
        id: u32,
        llocal: bool,
        rlocal: bool,
        members: Vec<SeqNode>,
    },
}

#[derive(Debug)]
pub struct SeqGraph {
    root: SeqNode,
    names: BTreeMap<u32, String>,
}

/// GraphPath - used for Iterator trait on SeqGraph
#[derive(Debug)]
pub struct GraphPath<'a> {
    curr_node: &'a SeqNode,
    path: Vec<u32>,
    pos: usize,
}

#[derive(Debug)]
pub enum SeqErr {
    BadJsonElement,
    Ambiguous,
    StringExpected,
    Unsupported,
    Mismatch,
}

impl SeqNode {
    fn dispatch(
        idx: &mut u32,
        names: &mut BTreeMap<u32, String>,
        elem: &JSON_Val,
        next: Rc<SeqNode>,
    ) -> Result<SeqNode, SeqErr> {
        *idx += 1;

        match elem {
            &JSON_Val::String(ref s) => SeqNode::read_str(idx, names, s, next),
            &JSON_Val::Object(ref map) => {
                // map is serde_json::Map<std::string::String, Value>
                SeqNode::read_obj(idx, names, map, next)
            }
            &JSON_Val::Array(ref l) => SeqNode::read_list(idx, names, l, 0),
            _ => Err(SeqErr::BadJsonElement),
        }
    }

    fn read_str(
        idx: &mut u32,
        names: &mut BTreeMap<u32, String>,
        s: &String,
        next: Rc<SeqNode>,
    ) -> Result<SeqNode, SeqErr> {
        Ok(SeqNode::Frag {
            id: *idx,
            val: Sequence::from_str(&s).unwrap(),
            llocal: false,
            rlocal: false,
            next: next.clone(),
        })
    }

    fn read_list(
        idx: &mut u32,
        names: &mut BTreeMap<u32, String>,
        l: &Vec<JSON_Val>,
        pos: usize,
    ) -> Result<SeqNode, SeqErr> {
        *idx += 1;
        if pos < l.len() {
            let next = SeqNode::read_list(idx, names, l, pos + 1).unwrap();
            SeqNode::dispatch(idx, names, &l[pos], Rc::new(next))
        } else {
            Ok(SeqNode::Nil)
        }
    }

    fn read_obj(
        idx: &mut u32,
        names: &mut BTreeMap<u32, String>,
        map: &serde_json::Map<std::string::String, JSON_Val>,
        next: Rc<SeqNode>,
    ) -> Result<SeqNode, SeqErr> {
        let _ = match map.get("id") {
            Some(s) => match s {
                &JSON_Val::String(ref s2) => {
                    names.insert(*idx, s2.clone());
                    ()
                }
                _ => return Err(SeqErr::StringExpected),
            },
            None => (),
        };

        if map.contains_key("seq") {
            if map.contains_key("dist") || map.contains_key("branch") {
                Err(SeqErr::Ambiguous)
            } else {
                match map.get("seq").unwrap() {
                    &JSON_Val::String(ref s2) => Ok(SeqNode::Frag {
                        id: *idx,
                        val: Sequence::from_str(&s2).unwrap(),
                        llocal: false,
                        rlocal: false,
                        next: next.clone(),
                    }),
                    _ => Err(SeqErr::BadJsonElement),
                }
            }
        } else if map.contains_key("branch") {
            if map.contains_key("dist") {
                Err(SeqErr::Ambiguous)
            } else {
                let members = match map.get("branch").unwrap() {
                    &JSON_Val::Array(ref _l) => {
                        let l: &Vec<JSON_Val> = _l;
                        let mut m = Vec::new();
                        for x in l {
                            *idx += 1;
                            m.push(SeqNode::dispatch(idx, names, x, next.clone()).unwrap());
                        }
                        Ok(m)
                        //Ok( l.iter().map(|x| SeqNode::dispatch( idx, names, x, next ).unwrap() ).collect() )
                    }
                    _ => Err(SeqErr::BadJsonElement),
                }
                .unwrap();
                *idx += 1;
                Ok(SeqNode::Branch {
                    id: *idx,
                    members: members,
                    llocal: false,
                    rlocal: false,
                })
            }
        } else if map.contains_key("dist") {
            Err(SeqErr::Unsupported)
        } else {
            Err(SeqErr::BadJsonElement)
        }
    }

    pub fn iden(&self) -> Option<u32> {
        match self {
            &SeqNode::Nil => None,
            &SeqNode::Frag { id: ref id, .. } => Some(id.clone()),
            &SeqNode::Branch { id: ref id, .. } => Some(id.clone()),
        }
    }
}

enum GraphAlnMode {
    Global,
    LocalTest,
    LocalFollow,
}

impl SeqGraph {
    pub fn from_json(serialized: &str) -> Result<SeqGraph, SeqErr> {
        let mut names: BTreeMap<u32, String> = BTreeMap::new();
        let mut idx = 0;

        let value = serde_json::from_str(serialized).unwrap();
        let mut idx = 0u32;
        let tree = SeqNode::dispatch(&mut idx, &mut names, &value, Rc::new(SeqNode::Nil)).unwrap();
        Ok(SeqGraph {
            root: tree,
            names: names,
        })
    }

    fn _max_len(n: &SeqNode) -> usize {
        match n {
            &SeqNode::Nil => 0,
            &SeqNode::Frag {
                val: ref val,
                next: ref next,
                ..
            } => val.len() + SeqGraph::_max_len(next),
            &SeqNode::Branch {
                members: ref members,
                ..
            } => members.iter().map(|n| SeqGraph::_max_len(n)).max().unwrap(),
        }
    }
    /// max_len - the maximum sequence length encoded by this graph
    pub fn max_len(&self) -> usize {
        SeqGraph::_max_len(&self.root)
    }

    fn _align(
        node: &SeqNode,
        query: &Sequence,
        m: &mut Matrix<Cell>,
        base_params: &AlnParams,
        start: i32,
        path: &mut Vec<u32>,
        pos: usize,
        mode: GraphAlnMode,
    ) -> Option<(/*score*/ i32, /*path*/ Vec<u32>)> {
        match node {
            &SeqNode::Nil => Some((0, path.clone())),
            &SeqNode::Frag {
                id: ref id,
                val: ref val,
                next: ref next,
                ..
            } => {
                path.push(id.clone());
                let t = align_matrix(val, query, &base_params, Some(start), m).unwrap();

                match mode {
                    GraphAlnMode::Global => {
                        let score = Cell::unpack(&m[t]).unwrap().1;
                        let (next_score, next_path) = SeqGraph::_align(
                            next,
                            query,
                            m,
                            base_params,
                            start + val.len() as i32,
                            path,
                            pos + 1,
                            mode,
                        )
                        .unwrap();
                        Some((max(score, next_score), next_path.to_vec()))
                    }
                    GraphAlnMode::LocalFollow => SeqGraph::_align(
                        next,
                        query,
                        m,
                        base_params,
                        start + val.len() as i32,
                        path,
                        pos + 1,
                        GraphAlnMode::LocalFollow,
                    ),
                    GraphAlnMode::LocalTest => {
                        Some((Cell::unpack(&m[t]).unwrap().1, path.to_vec()))
                    }
                }
            }
            &SeqNode::Branch {
                members: ref members,
                ..
            } => {
                path.push(node.iden().unwrap());

                let mut best_node = &SeqNode::Nil;
                let mut best_score = i32::min_value();

                for n in members {
                    path.truncate(pos + 1);
                    let (next_score, _) = SeqGraph::_align(
                        n,
                        query,
                        m,
                        base_params,
                        start,
                        path,
                        pos + 1,
                        match mode {
                            GraphAlnMode::Global => GraphAlnMode::Global,
                            _ => GraphAlnMode::LocalTest,
                        },
                    )
                    .unwrap();
                    if next_score > best_score {
                        best_score = next_score;
                        best_node = &n;
                    }
                }
                path.truncate(pos + 1);
                SeqGraph::_align(best_node, query, m, base_params, start, path, pos + 1, mode)
            }
        }
    }

    /// SeqGraph::align__global_max -- align query to graph, testing every possible branch to
    ///   find the global maximum.
    pub fn align__global_max(
        &self,
        query: &Sequence,
        base_params: &AlnParams,
    ) -> Option<(
        /*path*/ Vec<u32>,
        /*padded_ref*/ Sequence,
        /*padded_query*/ Sequence,
    )> {
        let mut _path = Vec::new();
        let ref_len = self.max_len();
        let mut m = Matrix::<Cell>::new(Cell(0), ref_len + 2, query.len() + 2);
        match SeqGraph::_align(
            &self.root,
            query,
            &mut m,
            base_params,
            0,
            &mut _path,
            0,
            GraphAlnMode::Global,
        ) {
            Some((score, path)) => {
                let mut full_ref_v = Vec::new();
                for s in GraphPath::from_graph(self, path.to_vec()) {
                    full_ref_v.extend(s.0);
                }
                let full_ref = Sequence(full_ref_v);

                // FIXME - we can't use the existing alignment matrix, because we don't
                //   know what path was tested last, and it's unlikely to be the maximum.
                //   however, we can probably do better than this
                let padded = align(&full_ref, query, base_params).unwrap();
                Some((path, padded.0, padded.1))
            }
            None => None,
        }
    }

    /// SeqGraph::align__local_max -- align query to graph, testing each branch to a depth of 1
    ///   to quickly find a maximum
    pub fn align__local_max(
        &self,
        query: &Sequence,
        base_params: &AlnParams,
    ) -> Option<(
        /*path*/ Vec<u32>,
        /*padded_ref*/ Sequence,
        /*padded_query*/ Sequence,
    )> {
        let mut _path = Vec::new();
        let ref_len = self.max_len();
        let mut m = Matrix::<Cell>::new(Cell(0), ref_len + 2, query.len() + 2);
        match SeqGraph::_align(
            &self.root,
            query,
            &mut m,
            base_params,
            0,
            &mut _path,
            0,
            GraphAlnMode::LocalFollow,
        ) {
            Some((score, path)) => {
                let mut full_ref_v = Vec::new();
                for s in GraphPath::from_graph(self, path.to_vec()) {
                    full_ref_v.extend(s.0);
                }
                let full_ref = Sequence(full_ref_v);

                // FIXME - we can't use the existing alignment matrix, because we don't
                //   know what path was tested last, and it's unlikely to be the maximum.
                //   however, we can probably do better than this
                let padded = align(&full_ref, query, base_params).unwrap();
                Some((path, padded.0, padded.1))
            }
            None => None,
        }
    }
}

impl<'a> GraphPath<'a> {
    pub fn from_graph(graph: &'a SeqGraph, path: Vec<u32>) -> GraphPath {
        GraphPath {
            curr_node: &graph.root,
            path: path,
            pos: 0,
        }
    }

    fn _next(&mut self) -> Option<Sequence> {
        if self.pos < self.path.len() {
            assert_eq!(self.path[self.pos], self.curr_node.iden().unwrap());

            match self.curr_node {
                &SeqNode::Frag {
                    val: ref val,
                    next: ref next,
                    ..
                } => {
                    self.pos += 1;
                    self.curr_node = next;
                    Some(val.clone())
                }
                &SeqNode::Branch {
                    members: ref members,
                    ..
                } => {
                    assert_eq!(self.curr_node.iden(), Some(self.path[self.pos]));
                    self.pos += 1;
                    for n in members {
                        if n.iden() == Some(self.path[self.pos]) {
                            self.curr_node = n;
                            return self._next();
                        }
                    }
                    panic!("couldn't find thing");
                }
                &SeqNode::Nil => None,
            }
        } else {
            None
        }
    }
}

impl<'a> Iterator for GraphPath<'a> {
    type Item = Sequence;

    fn next(&mut self) -> Option<Sequence> {
        if self.pos < self.path.len() {
            self._next()
        } else {
            None
        }
    }
}
