use clap::{Arg, App};
use clam;

use regex::Regex;
fn main(){
    let matches = App::new("clam")
        .version("1.0")
        .author("B. Reber <breber@cs.rochester.edu> 
M. Gould <mdg2838@rit.edu>")
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
        .arg(Arg::new("SET ASSOCIATIVITY")
            .short('S')
            .takes_value(true)
            .default_value("0")
            .required(false)
            .about("set associativity of the cache being targeted"))
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
        .arg(Arg::new("LLT_SIZE")
            .short('L')
            .default_value("128")
            .about("Number of elements in the lease lookup table")
            .required(false))
        .arg(Arg::new("MEM_SIZE")
            .short('M')
            .default_value("65536")
            .about("total memory allocated for lease information"))
        .arg(Arg::new("DISCRETIZE_WIDTH")
        .short('D')
        .default_value("9")
        .about("bit width avaiable for discretized short lease probability")
            .required(false))
        .arg(Arg::new("DEBUG")
            .short('d')
            .takes_value(false)
            .required(false)
            .about("enable even more information about lease assignment"))
         .arg(Arg::new("SAMPLING_RATE")
             .short('S')
             .default_value("256")
             .about("benchmark sampling rate")
             .required(false))
         .arg(Arg::new("EMPIRICAL_SAMPLE_RATE")
             .short('E')
             .default_value("yes")
             .about("Use given or empirically derived sampling rate")
             .required(false)).get_matches();
    

    let cache_size=matches.value_of("CACHE_SIZE").unwrap().parse::<u64>().unwrap();
    let perl_bin_num =matches.value_of("PRL").unwrap().parse::<u64>().unwrap();
    let llt_size=matches.value_of("LLT_SIZE").unwrap().parse::<u64>().unwrap();
    let mem_size=matches.value_of("MEM_SIZE").unwrap().parse::<u64>().unwrap();
    //get maximum number of scopes that can fit in given memory size with given llt size
    let max_scopes=mem_size/((2*llt_size + 16)*4);
    let discretize_width=matches.value_of("DISCRETIZE_WIDTH").unwrap().parse::<u64>().unwrap();
    let verbose = matches.is_present("VERBOSE");
    let debug   = matches.is_present("DEBUG");
    let cshel   = matches.is_present("CSHEL");
    let prl =matches.occurrences_of("PRL")>0;
    let mut output_file_name:String;
    let re= Regex::new(r"/(clam|shel).*/(.*?)\.txt$").unwrap();
    let search_string=matches.value_of("INPUT").unwrap().to_lowercase();
    let cap= re.captures(&*search_string).unwrap();
    let empirical_rate=matches.value_of("EMPIRICAL_SAMPLE_RATE").unwrap().to_lowercase();

//if associativity not specified, set as fully associative
    let num_ways:u64;
    let num_ways_given=matches.value_of("SET ASSOCIATIVITY").unwrap().parse::<u64>().unwrap();
    if num_ways_given==0{
        num_ways=cache_size;
    }
    else if  num_ways_given >cache_size {
        println!("The number of ways exceeds number of blocks in cache");
        panic!();
    }
    else {
        num_ways=num_ways_given;
    }
    //get mask for set bits
    let set_mask=(cache_size as f64 /num_ways as f64) as u32-1;

//generate distributions
       let (ri_hists,samples_per_phase,misses_from_first_access,empirical_sample_rate) = clam::io::build_ri_hists(matches.value_of("INPUT").unwrap(),cshel,set_mask);

//if specified used empirical sampling rate else use given or default rate.
    let sample_rate=if empirical_rate=="no" {matches.value_of("SAMPLING_RATE").unwrap().parse::<u64>().unwrap()}
    else {empirical_sample_rate};
   
   //generates PRL
   if prl {
    //generate bins
    let (binned_ri_distributions,binned_freqs,bin_width) = clam::io::get_binned_hists(matches.value_of("INPUT").unwrap(),perl_bin_num,set_mask);
    //compose output file name
    //this panic here avoids the almost certain panic that will result from running PRL on multi phase sampling files
    if &cap[1]=="shel"{
        panic!("Error! You can only use prl on sampling files with a single phase!");
    }
    output_file_name=format!("{}/{}_{}_{}",matches.value_of("OUTPUT").unwrap(),&cap[2],"prl","leases");
    //generate prl leases
    let (leases, dual_leases, lease_hits,trace_length) = clam::lease_gen::prl(bin_width,
        &ri_hists,&binned_ri_distributions,&binned_freqs,sample_rate,cache_size,discretize_width,&samples_per_phase,verbose,debug,set_mask).unwrap();
    println!("running PRL");
    let lease_vectors=clam::io::dump_leases(leases,dual_leases,lease_hits,trace_length,&output_file_name[..],sample_rate,misses_from_first_access);
    let output_lease_file_name=format!("{}/{}_{}_{}",matches.value_of("OUTPUT").unwrap(),&cap[2],"prl","lease.c");
    clam::io::gen_lease_c_file(lease_vectors,llt_size,max_scopes,mem_size,output_lease_file_name,discretize_width);
   }
 
   println!("running {}",&cap[1]);
    output_file_name=format!("{}/{}_{}_{}",matches.value_of("OUTPUT").unwrap(),&cap[2],&cap[1],"leases");
       //generates based on input file phases, CLAM or SHEL 
   let (leases,dual_leases, lease_hits,trace_length) = clam::lease_gen::shel_cshel(false,&ri_hists,cache_size,sample_rate,&samples_per_phase,discretize_width, 
    verbose,debug,set_mask).unwrap();
  let lease_vectors=clam::io::dump_leases(leases,dual_leases,lease_hits,trace_length,&output_file_name[..],sample_rate,misses_from_first_access);
   let output_lease_file_name=format!("{}/{}_{}_{}",matches.value_of("OUTPUT").unwrap(),&cap[2],&cap[1],"lease.c");
   //generate lease file
 clam::io::gen_lease_c_file(lease_vectors,llt_size,max_scopes,mem_size,output_lease_file_name,discretize_width);
  //generate CSHEL if option specified
   if cshel {
    //generate leases
    println!("Running C-SHEL.");
     let (leases, dual_leases, lease_hits,trace_length) = clam::lease_gen::shel_cshel(true,&ri_hists,cache_size,sample_rate,&samples_per_phase,discretize_width,
        verbose,debug,set_mask).unwrap();
       //compose output file name
    output_file_name=format!("{}/{}_{}_{}",matches.value_of("OUTPUT").unwrap(),&cap[2],"c-shel","leases");
      //output to file
      let output_lease_file_name=format!("{}/{}_{}_{}",matches.value_of("OUTPUT").unwrap(),&cap[2],"c-shel","lease.c");
    let lease_vectors=clam::io::dump_leases(leases,dual_leases,lease_hits,trace_length,&output_file_name[..],sample_rate,misses_from_first_access);
     clam::io::gen_lease_c_file(lease_vectors,llt_size,max_scopes,mem_size,output_lease_file_name,discretize_width);
    }

}






