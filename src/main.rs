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
        .arg(Arg::new("CACHE_SIZE"))
        .get_matches();

    let s  = matches.value_of("CACHE_SIZE").unwrap().to_string();
    let cache_size=s.parse::<u64>().unwrap();
    let leases:HashMap<u64,u64>;
    let dual_leases:HashMap<u64,(f32,u64)>;
    let lease_hits:HashMap<u64,HashMap<u64,u64>>;
    let sample_rate = 256;
    let perl_bin_num =5;
   let trace_length: u64;
    let verbose = true;
    let debug   = true;
    let cshel   = false;
    let PRL =true;
    let (binned_ri_distributions,binned_freqs,bin_width) = cshel::io::get_binned_hists(matches.value_of("INPUT").unwrap(),perl_bin_num);
      let (ri_hists,samples_per_phase) = cshel::io::build_ri_hists(matches.value_of("INPUT").unwrap());
   if PRL{
    let (leases_temp, dual_leases_temp, lease_hits_temp,trace_length_temp) = cshel::lease_gen::PRL(bin_width,&ri_hists,&binned_ri_distributions,&binned_freqs,256,cache_size,samples_per_phase).unwrap();
    leases=leases_temp;
    dual_leases=dual_leases_temp;
    lease_hits=lease_hits_temp;
    trace_length=trace_length_temp;
   }
   else{
   let (leases_temp, dual_leases_temp, lease_hits_temp,trace_length_temp) = cshel::lease_gen::shel_cshel(cshel,&ri_hists,cache_size,sample_rate,samples_per_phase,verbose,debug).unwrap();
   leases=leases_temp;
    dual_leases=dual_leases_temp;
    lease_hits=lease_hits_temp;
    trace_length=trace_length_temp;
}
   
    cshel::io::dump_leases(leases,dual_leases,lease_hits,trace_length);
}






