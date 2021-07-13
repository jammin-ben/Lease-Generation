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
        .get_matches();

    let cache_size = 128;
    let sample_rate = 256;

    let (ri_hists,samples_per_phase) = cshel::build_ri_hists(matches.value_of("INPUT").unwrap());
    let (leases, dual_leases) = cshel::gen_leases_c_shel(&ri_hists,cache_size,sample_rate,samples_per_phase,true).unwrap();
    cshel::dump_leases(leases,dual_leases);


}
