![Slime mold growth](https://tedideas.files.wordpress.com/2014/06/gif4b.gif?w=1000&h=563)

# Myxogast
## A graph aligner in Rust

Myxogast will be a network aligner based on the [classic dynamic programming](http://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm) [algorithms](https://en.wikipedia.org/wiki/Viterbi_algorithm), and using a JSON-based reference sequence format supporting:
 * wildcards ('N' and '*')
 * branch points
 * HMM-style probabilities

For example, this is a valid graph-JSON reference:
```JSON
    [{"id":    "domain_1",
     "seq":    "*ATGCATGC"},
     {"id":    "domain_2",
      "branch": [{"id":   "A",
                  "seq":  "GGCGGC"},
                 {"id":   "B",
                  "seq":  "TTATAG"}]
     },
     {"id":    "domain_3",
      "dist":  [{"A": 0.5, "T": 0.5}, {"A": 0.1, "T": 0.1, "G": 0.8}]},
     "TATTATA",

     {"__META__": {"version": 0.1,
                   "gap-penalty": -1}}]
```

Broadly, there are 3 node types: sequence fragments, branches, and probability distributions.  The simplest valid graph-JSON reference sequence is a simple JSON string, eg `"AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTG"`, and represents a "fragment" node without an id.  Functionally, this is equivalent to `{"id": "BA000007.2", "seq": "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTG"}`.

## Status
Currently, basic Needleman-Wunsch alignment works, basic graphs are parsed correctly, and some basic form of graph alignment seems to work.


## References
Myxogast's core graph-alignment algorithm is based on an algorithm discussed in section 12.4.2 (p283) of Gusfield's ["Algorithms on Strings, Trees and Sequences: Computer Science and Computational Biology"](https://www.amazon.com/Algorithms-Strings-Trees-Sequences-Computational/dp/0521585198), where an alignment of two suffix trees is discussed.  Another variation combining trees and dynamic programming is discussed later, in section 17.7 (p474).

Internally, many existing tools align a tree structure to a reference, including the most popular aligner used in NGS analysis, [bwa](http://bio-bwa.sourceforge.net/).  Sequence-graph representations include [FASTG](http://fastg.sourceforge.net/).


## Internals
Reference sequence graphs are represented like so:
```Rust
enum SeqNode {
    Nil,

    Frag {
        id: u32,
        val: Sequence,
        llocal: bool,
        rlocal: bool,

        next: Rc<SeqNode>
    },

    //Dist { id: u32, scores: ProbMatr },

    Branch {
        id: u32,
        llocal: bool,
        rlocal: bool,
        members: Vec<SeqNode>
    }
}

struct SeqGraph {
    tree: SeqNode,
    names: BTreeMap<u32, String>
}
```
Nodes are identified by a unique u32 integer.  Names from the "id" field in the graph-JSON input are stored in a BTreeMap, and not with the nodes, because (1) names are optional, and (2) not guaranteed to be unique.


## TODO
* set Frag llocal/rlocal fields based on splats
* add __PARAMS__ block to JSON format
* enforce short-circuiting of alignment (eg, max 3 indels)
* HMM/Viterbi alignmed
* modify Sequence to bit-pack bases (ATGCN- -> 3 bits, or 10 bases per u32) (?? is this worth the effort?)

