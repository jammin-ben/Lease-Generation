use clap::{Arg, App};
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
        .arg(Arg::new("CACHE_SIZE")
            .short('s')
            .multiple(false)
            .takes_value(true)
            .required(true))
        .get_matches();

    let s  = matches.value_of("CACHE_SIZE").unwrap().to_string();
    let cache_size=s.parse::<u64>().unwrap();

    let sample_rate = 256;


    let verbose = true;
    let debug   = true;
    let cshel   = false;

    let (ri_hists,samples_per_phase) = cshel::io::build_ri_hists(matches.value_of("INPUT").unwrap());
    let (leases, dual_leases, predicted_misses) = cshel::lease_gen::shel_cshel(cshel,&ri_hists,cache_size,sample_rate,samples_per_phase,verbose,debug).unwrap();
    println!("Dump predicted miss count (no contention misses): 
{}",predicted_misses);
    println!("Dump leases");
    cshel::io::dump_leases(leases,dual_leases);
}
