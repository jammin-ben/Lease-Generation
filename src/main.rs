use clap::{Arg, App};
use std::collections::HashMap;
use csv::{ReaderBuilder};
use serde::{Serialize, Deserialize};
use std::fmt;
use std::i64;
use cshel;

fn main(){
    let matches = App::new("cshel")
        .version("1.0")
        .author("B. Reber <breber@cs.rochester.edu>")
        .about("Lease assignment generator for phased traces")
        .arg(Arg::new("INPUT")
            .about("Sets the input file name")
            .required(true)
            .index(1))
        .arg(Arg::new("OUTPUT")
            .short('o')
            .multiple(true)
            .takes_value(true)
            .about("Sets the output file name"))
        .get_matches();

    let ri_hists = cshel::build_ri_hists(matches.value_of("INPUT").unwrap());
}
