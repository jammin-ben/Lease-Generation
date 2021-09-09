use clap::{Arg, App};
use cshel;
use std::collections::HashMap;
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
    let leases:HashMap<u64,u64>;
    let dual_leases:HashMap<u64,(f32,u64)>;
    let predicted_misses:u64;
    let sample_rate = 256;
    let perl_bin_num =5;
   
    let verbose = true;
    let debug   = true;
    let cshel   = false;
    let PRL =true;
    let (binned_ri_distributions,binned_freqs,bin_width) = cshel::io::get_binned_hists(matches.value_of("INPUT").unwrap(),perl_bin_num);
      let (ri_hists,samples_per_phase) = cshel::io::build_ri_hists(matches.value_of("INPUT").unwrap(),cshel);
   if PRL{
    let (x, y, z) = cshel::lease_gen::PRL(bin_width,&ri_hists,&binned_ri_distributions,&binned_freqs,256,cache_size,samples_per_phase).unwrap();
    leases=x;
    dual_leases=y;
    predicted_misses=z;
   }
   else{
   let (x, y, z) = cshel::lease_gen::shel_cshel(cshel,&ri_hists,cache_size,sample_rate,samples_per_phase,verbose,debug).unwrap();
   leases=x;
    dual_leases=y;
    predicted_misses=z;
}
    println!("Dump predicted miss count (no contention misses): 
{}",predicted_misses);
    cshel::io::dump_leases(leases,dual_leases);
}






