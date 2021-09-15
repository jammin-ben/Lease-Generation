use clap::{Arg, App};
use cshel;

use regex::Regex;
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
            .index(2)
            .takes_value(true)
            .about("Sets the output file Location")
            .required(true))
        .arg(Arg::new("CACHE_SIZE")
            .short('s')
            .takes_value(true)
            .required(true)
            .about("target cache size for algorithms"))
        .arg(Arg::new("PRL")
            .short('p')
            .default_missing_value("5")
            .default_value("5")
            .required(false)
            .about("calculate leases for prl (only for non_phased sampling files)"))
        .arg(Arg::new("CSHEL")
            .short('c')
            .required(false)
            .about("calculate leases for CSHEL"))
        .arg(Arg::new("VERBOSE")
            .short('V')
            .takes_value(false)
            .required(false)
            .about("output information about lease assignment"))
        .arg(Arg::new("SAMPLING_RATE")
            .short('S')
            .default_value("256")
            .about("benchmark sampling rate")
            .required(false))
        .arg(Arg::new("DEBUG")
            .short('d')
            .takes_value(false)
            .required(false)
            .about("enable even more information about lease assignment")).get_matches();

    let cache_size=matches.value_of("CACHE_SIZE").unwrap().parse::<u64>().unwrap();
    let sample_rate = matches.value_of("SAMPLING_RATE").unwrap().parse::<u64>().unwrap();
    let perl_bin_num =matches.value_of("PRL").unwrap().parse::<u64>().unwrap();
    let verbose = matches.is_present("VERBOSE");
    let debug   = matches.is_present("DEBUG");
    let cshel   = matches.is_present("CSHEL");
    let prl =matches.occurrences_of("PRL")>0;
    let mut output_file_name:String;
    let re= Regex::new(r"/(clam|shel).*/(.*?)\.txt$").unwrap();
    let search_string=matches.value_of("INPUT").unwrap().to_lowercase();
    let cap= re.captures(&*search_string).unwrap();

//generate distributions
       let (ri_hists,samples_per_phase,misses_from_first_access) = cshel::io::build_ri_hists(matches.value_of("INPUT").unwrap(),cshel);
   //generates PRL
   if prl {
    //generate bins
    let (binned_ri_distributions,binned_freqs,bin_width) = cshel::io::get_binned_hists(matches.value_of("INPUT").unwrap(),perl_bin_num);
    //compose output file name
    //this panic here avoids the almost certain panic that will result from running PRL on multi phase sampling files
    if &cap[1]=="shel"{
        println!("Error! You can only use prl on sampling files with a single phase!");
        panic!();
    }
    output_file_name=format!("{}/{}_{}_{}",matches.value_of("OUTPUT").unwrap(),&cap[2],"prl","leases");
    //generate prl leases
    let (leases, dual_leases, lease_hits,trace_length) = cshel::lease_gen::prl(bin_width,
        &ri_hists,&binned_ri_distributions,&binned_freqs,256,cache_size,&samples_per_phase,verbose).unwrap();
    
    cshel::io::dump_leases(leases,dual_leases,lease_hits,trace_length,&output_file_name[..],misses_from_first_access);
   }
 
   
    output_file_name=format!("{}/{}_{}_{}",matches.value_of("OUTPUT").unwrap(),&cap[2],&cap[1],"leases");
       //generates based on input file phases, CLAM or SHEL 
   let (leases,dual_leases, lease_hits,trace_length) = cshel::lease_gen::shel_cshel(false,&ri_hists,cache_size,sample_rate,&samples_per_phase,verbose,debug).unwrap();
  cshel::io::dump_leases(leases,dual_leases,lease_hits,trace_length,&output_file_name[..],misses_from_first_access);

  //generate CSHEL if option specified
   if cshel {
    //generate leases
     let (leases, dual_leases, lease_hits,trace_length) = cshel::lease_gen::shel_cshel(true,&ri_hists,cache_size,sample_rate,&samples_per_phase,verbose,debug).unwrap();
       //compose output file name
     output_file_name=format!("{}/{}_{}_{}",matches.value_of("OUTPUT").unwrap(),&cap[2],"c-shel","leases");
      //output to file
    cshel::io::dump_leases(leases,dual_leases,lease_hits,trace_length,&output_file_name[..],misses_from_first_access);
    }

}






