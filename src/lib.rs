



// Functions for parsing input files, debug prints, 
// and lease output
pub mod io {
    use serde::{Serialize, Deserialize};
    use std::collections::BinaryHeap;
    use std::collections::HashMap;
    use std::io::Write;
    



    #[derive(Deserialize)]
    #[derive(Debug)]
    struct Sample{
        phase_id_ref: String,
        ri: String,
        tag: String,
        time: u64,
    }

    #[derive(Serialize)]
    struct PhaseTime{
        phase_id: u16,
        start_time: u64,
    }
    pub fn get_binned_hists(input_file:&str, num_bins: u64,set_mask:  u32) -> (super::lease_gen::BinnedRIs,super::lease_gen::BinFreqs,u64){

        let mut curr_bin: u64= 0;
        let mut curr_bin_dict=HashMap::<u64,u64>::new();
        let mut bin_freqs=HashMap::<u64,HashMap<u64,u64>>::new();
        let mut bin_ri_distributions=HashMap::<u64,HashMap<u64,HashMap<u64,u64>>>::new();
        let mut curr_ri_distribution_dict=HashMap::<u64,HashMap<u64,u64>>::new();
        let mut last_address: u64 =0;
        let mut all_keys: Vec<u64>=Vec::new();
        
        bin_freqs.insert(0,curr_bin_dict.clone());
        bin_ri_distributions.insert(0,curr_ri_distribution_dict.clone());
       
        let mut rdr = csv::ReaderBuilder::new()
            .has_headers(false)
            .from_path(input_file)
            .unwrap();
        for result in rdr.deserialize() {
            let sample: Sample = result.unwrap();
            last_address= sample.time;
        }
        let bin_width = ((last_address as f64 / (num_bins as f64)) as f64).ceil() as u64;
        let mut rdr = csv::ReaderBuilder::new()
            .has_headers(false)
            .from_path(input_file)
            .unwrap();
        for result in rdr.deserialize(){
            let sample: Sample = result.unwrap();

            
            //if outside of current bin, moved to the next
            if sample.time>curr_bin+bin_width{
              //store the RI and frequency data for the old bin
                bin_freqs.insert(curr_bin,curr_bin_dict.clone());
                bin_ri_distributions.insert(curr_bin,curr_ri_distribution_dict.clone());
                //initalize storage for the new bin
                curr_bin_dict.clear();
                curr_ri_distribution_dict.clear();
                curr_bin+=bin_width;

            }
            let tag = u32::from_str_radix(&sample.tag,16).unwrap();
                let set = (tag&set_mask) as u64;
                
            let addr = u64::from_str_radix(&sample.phase_id_ref,16).unwrap()|(set<<32);
            let ri = u64::from_str_radix(&sample.ri,16).unwrap();
            
            if curr_bin_dict.contains_key(&addr){
                curr_bin_dict.insert(addr,curr_bin_dict.get(&addr).unwrap()+1);
            }
            else {
                curr_bin_dict.insert(addr,1);
            }
            //if the address isn't a key in the outer level hashmap, add a hashmap corresponding to that key; then add the key value pair (ri,1) to that hashmap
            //if the address is a key in the outer level hashmap, check to see if the current ri is a key in the corresponding inner hashmap.
            //if yes, then increment the count for that RI, if not, insert a new key value pair (ri,1).
           *curr_ri_distribution_dict.entry(addr).or_insert(HashMap::new()).entry(ri).or_insert(0) += 1;

           //create an array of all found references
            if !all_keys.iter().any(|&i| i==addr){
                all_keys.push(addr);
            }
                
        }
       
        //store the frequency and RI data for the last bin
        bin_freqs.insert(curr_bin,curr_bin_dict.clone());
        bin_ri_distributions.insert(curr_bin,curr_ri_distribution_dict.clone());
        //if a reference is not in a bin, add it with a frequency count of 0
        let temp= bin_freqs.clone();
        for (bin,_addrs) in &temp{
            let bin_freqs_temp=bin_freqs.entry(*bin).or_insert(HashMap::new());
            for key in &all_keys{
                bin_freqs_temp.entry(*key).or_insert(0);
            }
           
        }

        (super::lease_gen::BinnedRIs::new(bin_ri_distributions),super::lease_gen::BinFreqs::new(bin_freqs),bin_width)
    }



    pub fn build_phase_transitions(input_file:&str) -> (Vec<(u64,u64)>,usize,u64){
        println!("Reading input from: {}", input_file);


        let mut rdr = csv::ReaderBuilder::new()
            .has_headers(false)
            .from_path(input_file)
            .unwrap();
        let mut u_tags=HashMap::<u64,bool>::new();
        let mut sample_hash = HashMap::new();
        let mut last_sample_time: u64=0;
        let mut sample_num: u64=0;

        for result in rdr.deserialize() {
            let sample: Sample = result.unwrap();
            let ri = u64::from_str_radix(&sample.ri,16).unwrap();

            
            //don't use end of benchmark infinite RIs
           

            //store unique tags
             u_tags.insert(u64::from_str_radix(&sample.tag,16).unwrap(),false);

        
            let phase_id_ref = u64::from_str_radix(&sample.phase_id_ref,16).unwrap();

            let phase_id = (phase_id_ref & 0xFF000000)>>24;

           
            let reuse_time = sample.time;
            let use_time;
            if (ri as i32)<0 {
               use_time=reuse_time -(!ri+1)as u64;
            }
            else{
                use_time = reuse_time - ri as u64;
            }
            sample_hash.insert(use_time,phase_id);
            //get empircal sampling rate
            last_sample_time=sample.time;
            sample_num+=1;

        }
        //empircally calculate sampling rate
        let sampling_rate=(last_sample_time as f64/sample_num as f64).round() as u64;
        println!("sampling_rate:{}",sampling_rate);
        //every data block is associated with at least one miss in the absense of hardware prefetching.
        let first_misses=u_tags.len();
       

        let mut sorted_samples: Vec<_> = sample_hash.iter().collect();
        sorted_samples.sort_by_key(|a| a.0);

        //get phase transitions
        let mut phase_transitions = HashMap::new(); //(time,phase start)
        let mut phase = 0;
        phase_transitions.insert(0,phase);

        for s in sorted_samples.iter(){
            if *s.1 != phase{
                phase_transitions.insert(*s.0,*s.1);
                phase = *s.1;
            } 
        }

        let mut sorted_transitions: Vec<_> = phase_transitions.iter().collect();
        sorted_transitions.sort_by_key(|a| a.0);

        return(sorted_transitions.iter().map(|&x| (*(x.0),*(x.1))).collect(),first_misses,sampling_rate);
      
    }

    //Build ri hists in the following form
    //{ref_id,
    //  {ri,
    //    (count,{phase_id,(head_cost,tail_cost)})}}
    //
    // Head cost refers to the accumulation of 
    // cost from reuses with length ri, which may span 
    // phase boundaries.
    
    // Tail cost refers to the accumulation of cost from reuses greater than ri.
    // This cost may span a phase boundary 
    pub fn build_ri_hists(input_file:&str,cshel:bool,set_mask:u32) -> (super::lease_gen::RIHists,HashMap<u64,u64>,usize,u64){
        let (phase_transitions,first_misses,sampling_rate) = build_phase_transitions(input_file);
        let mut rdr = csv::ReaderBuilder::new()
            .has_headers(false)
            .from_path(input_file)
            .unwrap();

        let mut ri_hists = HashMap::new();
        let mut samples_per_phase = HashMap::new();
     
        //Don't need to calculate head or tail costs for PRL or CLAM or SHEL, 
     

        if cshel {
            println!("before first pass!");
            for result in rdr.deserialize() {
                let sample: Sample = result.unwrap();
                let ri = u32::from_str_radix(&sample.ri,16).unwrap();
                let reuse_time = sample.time;
                let use_time;

                let mut ri_signed = ri as i32;

                 if ri_signed < 0 {
                   use_time=reuse_time -(!ri_signed+1)as u64;
                    ri_signed = i32::MAX; //canonical value for negatives
                }
                else {
                    use_time = reuse_time - ri_signed as u64;
                }

                let phase_id_ref = u64::from_str_radix(&sample.phase_id_ref,16).unwrap();
                
                 let tag = u32::from_str_radix(&sample.tag,16).unwrap();
                let set = (tag&set_mask) as u64;
                let set_phase_id_ref=phase_id_ref | set<<32;

                let next_phase_tuple = match super::helpers::binary_search(&phase_transitions,use_time){
                    Some(v) => v,
                    None => (reuse_time+1,0),
                };
                super::lease_gen::process_sample_head_cost(&mut ri_hists,
                                                       set_phase_id_ref,
                                                       ri_signed as u64,
                                                       use_time,
                                                       next_phase_tuple);

                let phase_id = (phase_id_ref & 0xFF000000)>>24;
                *samples_per_phase.entry(phase_id).or_insert_with(|| 0)+=1;
            }
            
            let mut rdr = csv::ReaderBuilder::new()
                .has_headers(false)
                .from_path(input_file)
                .unwrap();
                   println!("before second pass!");
            //second pass for tail costs
            for result in rdr.deserialize() {
                let sample: Sample = result.unwrap();
                let ri = u32::from_str_radix(&sample.ri,16).unwrap();

                let reuse_time = sample.time;
                let use_time;

                let mut ri_signed = ri as i32;

                if ri_signed < 0 {
                   use_time=reuse_time -(!ri_signed+1)as u64;
                    ri_signed = i32::MAX; //canonical value for negatives
                }
                else {
                    use_time = reuse_time - ri_signed as u64;
                }

                let phase_id_ref = u64::from_str_radix(&sample.phase_id_ref,16).unwrap();
                let tag = u32::from_str_radix(&sample.tag,16).unwrap();
                let set = (tag&set_mask) as u64;
                let set_phase_id_ref=phase_id_ref | set<<32;
               
                let next_phase_tuple = match super::helpers::binary_search(&phase_transitions,use_time){
                    Some(v) => v,
                    None => (reuse_time+1,0),
                };

                super::lease_gen::process_sample_tail_cost(&mut ri_hists,
                                         set_phase_id_ref,
                                         ri_signed as u64,
                                         use_time,
                                         next_phase_tuple);
            }
        }
        //if not doing C-SHEL generates RI distribution with head and tail costs set to 0 (significantly faster)
        else{
            for result in rdr.deserialize() {
                let sample: Sample = result.unwrap();
                let mut ri = u64::from_str_radix(&sample.ri,16).unwrap();
                let _reuse_time = sample.time;
                //if sample is negative, there is no reuse 
                ri = std::cmp::min(ri,i32::MAX as u64);

                let phase_id_ref = u64::from_str_radix(&sample.phase_id_ref,16).unwrap();
                let tag = u32::from_str_radix(&sample.tag,16).unwrap();
                let set = (tag&set_mask) as u64;
                let set_phase_id_ref=phase_id_ref | set<<32;
                let phase_id = (phase_id_ref & 0xFF000000)>>24;
                *samples_per_phase.entry(phase_id).or_insert(0)+=1;
                //if the reference isn't in the distrubtion add it
                //if an ri for that reference isn't in the distrubtion add it as key with a value of (ri_count,{phaseID,(0,0)))
                //if an ri for that reference is in the distribuntion increment the ri count by 1
               ri_hists.entry(set_phase_id_ref).or_insert(HashMap::new()).entry(ri).and_modify(|e| {e.0+=1}).or_insert((1,HashMap::new())).1.entry(phase_id).or_insert((0,0));
            }
        }
        (super::lease_gen::RIHists::new(ri_hists),samples_per_phase,first_misses,sampling_rate)
    }
    pub fn dump_leases(leases: HashMap<u64,u64>, 
        dual_leases: HashMap<u64,(f64,u64)>,
        lease_hits:HashMap<u64,HashMap<u64,u64>>,
        trace_length:u64,
        output_file:&str,
        sampling_rate:u64,
        first_misses:usize) ->Vec<(u64,u64,u64,u64,f64)> {
        let mut num_hits=0;
        //create lease output vector
        let mut lease_vector: Vec<(u64,u64,u64,u64,f64)> = Vec::new();
        for (&phase_address,&lease) in leases.iter() {
            let lease = if lease >0 {lease} else {1}; 
            let phase   = (phase_address & 0xFF000000)>>24;
            let address =  phase_address & 0x00FFFFFF;
            if dual_leases.contains_key(&phase_address){
               lease_vector.push((phase,address,lease,dual_leases.get(&phase_address).unwrap().1,1.0-dual_leases.get(&phase_address).unwrap().0));
            }
            else {
                lease_vector.push((phase,address,lease,0, 1.0));
            }
        } 
        lease_vector.sort_by_key(|a| (a.0,a.1)); //sort by phase and then by reference
        //get number of predicted misses
        for (_phase, address, lease_short, lease_long, percentage) in lease_vector.iter(){
            
            //we are assuming that our sampling captures all RIS by assuming the distribution is normal
            //thus if an RI for a reference didn't occur during runtime (i.e., the base lease of 1 that all references get) 
            //we can assume the number of hits it gets is zero, and moreover, even if that reuse interval does happen, we have no way
            if lease_hits.get(address).unwrap().get(lease_short)!=None{
                num_hits+=(*lease_hits.get(address).unwrap().get(lease_short).unwrap() as f64 *(percentage)).round() as u64;
            }
            if lease_hits.get(address).unwrap().get(lease_long)!=None{
                num_hits+=(*lease_hits.get(address).unwrap().get(lease_long).unwrap() as f64 *(1.0-percentage)).round() as u64;
            }
            
         }
         println!("Writing output to: {}",output_file);
         let mut file = std::fs::File::create(output_file).expect("create failed");
         file.write_all(&format!("Dump predicted miss count (no contention misses): {}\n",trace_length-num_hits*sampling_rate+first_misses as u64)[..].as_bytes()).expect("write failed");
         file.write_all("Dump formated leases\n".as_bytes()).expect("write failed");

         for (phase, address, lease_short, lease_long, percentage) in lease_vector.iter(){
             file.write_all(&format!("{:x}, {:x}, {:x}, {:x}, {}\n",
                                     phase, 
                                     address, 
                                     lease_short, 
                                     lease_long, 
                                     percentage)[..].as_bytes()).expect("write failed");
        }
        return lease_vector;
    }
  #[derive(Debug)]
    struct LeaseItem{
        phase: String,
        reference:   String,
        long_lease: String,
        short_lease: String,
        short_lease_prob: f64,
    }
// function for generating c-files
    pub fn gen_lease_c_file(
        mut lease_vector:Vec<(u64,u64,u64,u64,f64)>,llt_size:usize,mem_size:usize,output_file:String,discretize_width:u64){
    
    let mut phase_lease_arr:HashMap<u64,HashMap<u64,(u64,u64,f64,bool)>>=HashMap::new();
    let mut phases:Vec<u64>=Vec::new();
    for lease in lease_vector.iter(){
        if !phases.contains(&lease.0){
            phases.push(lease.0);
        }
    }
    //due to the way we store leases in memory, we can't skip any phases so 
    //if there are phases with no leases, assign a dummy lease to the skipped phase
    for phase in 0..*phases.iter().max().unwrap(){
        if !phases.contains(&phase){
            lease_vector.push((phase,0,0,0,1.0));
        }
    }



    //convert lease vector to hashmap of leases per phase
      for (phase, address, lease_short, lease_long, percentage) in lease_vector.iter(){
         phase_lease_arr.entry(*phase).or_insert(HashMap::new()).entry(*address).or_insert((*lease_short,*lease_long,*percentage,if lease_long >&0 {true} else {false})); 
      }
      let default_lease=1;

      //make sure each phase can fit in the specified LLT
      for (phase,phase_leases) in phase_lease_arr.iter(){
        if phase_leases.len()>llt_size {
            println!("Leases for Phase {} don't fit in lease lookup table!",phase);
            panic!();
        }

      }

      //make sure that all phases can fit in the memory allocated
      if (4*(2*llt_size)+16*4)*phase_lease_arr.len()>mem_size{
        println!("Error: phases cannot fit in specified {} byte memory",mem_size);
            panic!();
      }
      //write header
      let mut file = std::fs::File::create(output_file).expect("create failed");
      file.write_all("#include \"stdint.h\"\n\n".as_bytes()).expect("write failed");
      file.write_all(format!("static uint32_t lease[{}] __attribute__((section (\".lease\"))) __attribute__ ((__used__)) = {{\n",mem_size/4).as_bytes()).expect("write failed");
      file.write_all("// lease header\n".as_bytes()).expect("write failed");
      let mut phase_index:u64=0;//len returns usize which can't directly substituted as u64
       for i in 0..phase_lease_arr.len(){
            let phase_leases=phase_lease_arr.get(&phase_index).unwrap();
            phase_index=phase_index+1; //increment to next phase
            file.write_all(format!("// phase {}\n",i).as_bytes()).expect("write failed");
     
        let mut dual_lease_ref=(0,0,1.0);
                let mut lease_phase:Vec<(u64,u64)>=Vec::new();
                let dual_lease_found=false;
       for(lease_ref,lease_data) in phase_leases.iter(){
           //convert hashmap of leases for phase to vector
         lease_phase.push((*lease_ref,lease_data.0));
         //get dual lease if it exists;
         if lease_data.3==true && dual_lease_found==false {
            dual_lease_ref=(*lease_ref,lease_data.1,lease_data.2);
        }          
       }
      
       //output config
       for j in 0..16{
                if j==0{
                 file.write_all(format!("\t0x{:08x},\t// default lease\n",default_lease).as_bytes()).expect("write failed");
                }
                else if j==1{
                  file.write_all(format!("\t0x{:08x},\t// long lease value\n",dual_lease_ref.1).as_bytes()).expect("write failed");
                }
                else if j==2{
                    file.write_all(format!("\t0x{:08x},\t// short lease probability\n",discretize(dual_lease_ref.2,discretize_width)).as_bytes()).expect("write failed");
                }
                else if j==3{
file.write_all(format!("\t0x{:08x},\t// num of references in phase\n",phase_leases.len()).as_bytes()).expect("write failed");
  
                }
                else if j==4{
file.write_all(format!("\t0x{:08x},\t// dual lease ref (word address)\n",dual_lease_ref.0>>2).as_bytes()).expect("write failed");
  
                }
                else {
file.write_all(format!("\t0x{:08x},\t // unused\n",0).as_bytes()).expect("write failed");
                }
        }
        let field_list=["reference address","lease0 value"];
       

            // loop through lease fields
        for k in 0..2{
            file.write_all(format!("\t//{}\n\t",field_list[k]).as_bytes()).expect("write failed");
            

            for j in 0..llt_size{
                if j<phase_leases.len(){
                    if k==0 {
                        file.write_all(format!("0x{:08x}",lease_phase[j].0).as_bytes()).expect("write failed");
                    }
                    else {
                        file.write_all(format!("0x{:08x}",lease_phase[j].1).as_bytes()).expect("write failed");
                    }
                }
                else {
                    file.write_all(format!("0x{:08x}",0).as_bytes()).expect("write failed");
                }
                //print delimiter
                if j+1==llt_size && k==1 && i+1 ==phase_lease_arr.len(){
                    file.write_all(format!("\n").as_bytes()).expect("write failed");
                }
                else if j+1==llt_size{
                    file.write_all(format!(",\n").as_bytes()).expect("write failed");
                }
                else if ((j+1)%10)==0{
                    file.write_all(format!(",\n\t").as_bytes()).expect("write failed");
                }
                else {
                    file.write_all(format!(", ").as_bytes()).expect("write failed");
                }
   
            }
        }

    }
  file.write_all(format!("}};").as_bytes()).expect("write failed");
}


        pub fn discretize(percentage: f64,discretization:u64)->u64{
            let percentage_binary =((percentage*((2<<(discretization-1)) as f64)-1.0)).round() as u64;
            return percentage_binary;
        }    
    


    
    pub mod debug {
        pub fn print_ri_hists(rihists: &super::super::lease_gen::RIHists){
            for (ref_id, ref_ri_hist) in &rihists.ri_hists{
                println!("({},0x{:x}):",(ref_id & 0xFF000000) >> 24, ref_id & 0x00FFFFFF);
                for (ri,tuple) in ref_ri_hist{
                    println!(" | ri 0x{:x}: count {}",ri, tuple.0);
                    for (phase_id,cost) in &tuple.1{
                        println!(" | | phase {} head_cost {} tail_cost {}",
                                 phase_id,
                                 cost.0,
                                 cost.1);
                    }
                }
            }
        }

        pub fn destructive_print_ppuc_tree(ppuc_tree: &mut super::BinaryHeap<super::super::lease_gen::PPUC>){
            while ppuc_tree.peek() != None{
                println!("ppuc: {:?}",ppuc_tree.pop().unwrap());
            }
        }
    }
}

//Small miscellaneous functions used 
mod helpers {
    pub fn float_min(a: f64, b:f64) -> f64{
        if a.lt(&b){ 
            return a;
        }
        b
    }

    pub fn binary_search(vector: &Vec<(u64,u64)>,value:u64) -> Option<(u64,u64)>{
        let mut min = 0;
        let mut max = vector.len() - 1;

        while max>=min{
            let guess = (max + min)/2;
            if vector[guess].0==value{ //on transition sample
                if guess < vector.len() - 1{
                    return Some(vector[guess+1]);
                }
                //transition of last phase
                return None;
            }

            if vector[guess].0 < value{
                min = guess+1;
            }
            if vector[guess].0 > value{
                if guess >0{
                    max = guess-1;
                }
                else{
                    return Some(vector[guess])
                }
            }
        }
        assert_eq!(max,min-1);
        if min==vector.len(){
            return None;
        }

        Some(vector[min])
    }
}

//Core Algorithms
pub mod lease_gen {
   
    use std::collections::BinaryHeap;
    use std::collections::HashMap;
    use core::cmp::Ordering;
    #[derive(Debug,Clone)]
    pub struct BinFreqs{
        pub bin_freqs: HashMap <u64,HashMap<u64,u64>>,
          
    }
      #[derive(Debug,Clone)]
    pub struct BinnedRIs{
        pub bin_ri_distribution: HashMap <u64,HashMap<u64,HashMap<u64,u64>>>,
       
    }
    
    impl BinFreqs{
        pub fn new(bin_freqs_input: HashMap<u64,HashMap<u64,u64>>) -> Self{
            BinFreqs{
                bin_freqs:bin_freqs_input,
            } 
        }
    }

    impl BinnedRIs{
        pub fn new(bin_ri_input: HashMap <u64,HashMap<u64,HashMap<u64,u64>>>) -> Self{
            BinnedRIs{
                bin_ri_distribution:bin_ri_input,
            } 
        }
    }



    pub struct RIHists{
        pub ri_hists: HashMap<u64,HashMap<u64,(u64,HashMap<u64,(u64,u64)>)>>,
    }

    impl RIHists {
        pub fn new(ri_hists_input: HashMap<u64,HashMap<u64,(u64,HashMap<u64,(u64,u64)>)>>) -> Self{
            RIHists{
                ri_hists: ri_hists_input,
            }
        }

        pub fn get_ref_hist(&self,ref_id: u64) -> &HashMap<u64,(u64,HashMap<u64,(u64,u64)>)>{
            return self.ri_hists.get(&ref_id).unwrap();
        }

        pub fn get_ref_ri_count(&self,ref_id: u64,ri:u64) -> u64 {
            return self.ri_hists.get(&ref_id).unwrap().get(&ri).unwrap().0;
        }

        pub fn get_ref_ri_cost(&self,ref_id: u64,ri:u64) -> &HashMap<u64,(u64,u64)> {
            return &self.ri_hists.get(&ref_id).unwrap().get(&ri).unwrap().1;
        }

        pub fn get_ref_ri_phase_cost(&self,ref_id: u64, ri:u64, phase: u64) -> (u64,u64) {
            return *self.ri_hists.get(&ref_id).unwrap().get(&ri).unwrap().1.get(&phase).unwrap();
        }
    }

    #[derive(Debug,Copy,Clone)]
    pub struct PPUC {
        ppuc: f64,
        lease: u64,
        old_lease: u64,
        ref_id: u64,
        new_hits: u64,
    }
    impl PartialOrd for PPUC {
        fn partial_cmp(&self, other: &Self) -> Option<Ordering>{
            self.ppuc.partial_cmp(&other.ppuc)
        }
    }
    impl Ord for PPUC {
        fn cmp(&self, other: &Self) -> Ordering {
            other.ppuc.partial_cmp(&self.ppuc).unwrap()
        }
    }
    impl PartialEq for PPUC {
        fn eq(&self, other: &Self) -> bool{
            other.ppuc.eq(&self.ppuc)
        }
    }
    impl Eq for PPUC {}

    pub fn process_sample_head_cost(ri_hists: &mut HashMap<u64,HashMap<u64,(u64,HashMap<u64,(u64,u64)>)>>,
                      phase_id_ref: u64,
                      ri: u64,
                      use_time: u64,
                      next_phase_tuple: (u64,u64)){

        let phase_id = (phase_id_ref & 0xFF000000)>>24;

        //increment count
        let ref_hist = ri_hists.entry(phase_id_ref).or_insert_with(|| HashMap::new());
        let ri_tuple = ref_hist.entry(ri).or_insert_with(|| (0,HashMap::new()));
        ri_tuple.0  += 1;

         //increment head costs
        let this_phase_head_cost= std::cmp::min(next_phase_tuple.0 - use_time,ri);
        let next_phase_head_cost= std::cmp::max(use_time as i64 + ri as i64 - next_phase_tuple.0 as i64,
                                                0) as u64;
       
        
        ri_tuple.1.entry(phase_id)
            .or_insert_with(|| (0,0)).0 += this_phase_head_cost;

        if next_phase_head_cost>0 {
            ri_tuple.1.entry(next_phase_tuple.1)
                .or_insert_with(|| (0,0)).0 += next_phase_head_cost;
        }
    }

    pub fn process_sample_tail_cost(ri_hists: &mut HashMap<u64,HashMap<u64,(u64,HashMap<u64,(u64,u64)>)>>,
                      phase_id_ref: u64,
                      ri: u64,
                      use_time: u64,
                      next_phase_tuple: (u64,u64)){

        let phase_id = (phase_id_ref & 0xFF000000)>>24;

        let ref_hist = ri_hists.entry(phase_id_ref)
            .or_insert_with(|| HashMap::new());

        //this heinous code exists so we can iterate through a HashMap while modifying it
        let ris: Vec<&u64> = ref_hist.keys().collect();
        let mut ris_keys :Vec<u64> = Vec::new();
        for ri_other in ris{
            ris_keys.push(*ri_other);
        }

        //increment tail costs
        for ri_other in ris_keys {

            //no tail cost if the other ri is greater
            if ri_other >= ri {
                continue;
            }
            let count_phase_cost_tuple = ref_hist.entry(ri_other)
                .or_insert_with(|| (0,HashMap::new()));

            let this_phase_tail_cost = std::cmp::min(next_phase_tuple.0 - use_time,
                                                     ri_other);
            let next_phase_tail_cost = std::cmp::max(0,
                                                     use_time as i64 + 
                                                       ri_other as i64 - 
                                                       next_phase_tuple.0 as i64) as u64;

            count_phase_cost_tuple.1.entry(phase_id)
                .or_insert_with(|| (0,0)).1 += this_phase_tail_cost;

            if next_phase_tail_cost > 0 {
                count_phase_cost_tuple.1.entry(next_phase_tuple.1)
                    .or_insert_with(|| (0,0)).1 += next_phase_tail_cost;
            }
        }
    }

    fn cshel_phase_ref_cost(sample_rate: u64, phase: u64,ref_id: u64,old_lease: u64,new_lease: u64,ri_hists: &RIHists) -> u64 {
        let mut old_cost = 0;
        let mut new_cost = 0;
        //if a set in a phase does not have this reference
        if ri_hists.ri_hists.get(&ref_id)==None {
            return 0;
        }
        let ri_hist = ri_hists.ri_hists.get(&ref_id).unwrap();

        for (&ri,(_,phase_cost_hashmap)) in ri_hist.iter(){
            let (phase_head_cost,phase_tail_cost) = match phase_cost_hashmap.get(&phase) {
                Some((a,b)) => (*a,*b),
                None        => (0,0), 
            };
              if ri <= old_lease {
                old_cost += phase_head_cost;
            }
            if ri ==old_lease {
                old_cost += phase_tail_cost;
            }

            if ri <= new_lease {
                new_cost += phase_head_cost;
            }
            if ri ==new_lease {
                new_cost += phase_tail_cost;
            }
          
        }
      
        (new_cost - old_cost) * sample_rate   
    }

    fn shel_phase_ref_cost(sample_rate: u64, 
                           phase:      u64,
                           ref_id:      u64,
                           old_lease:   u64,
                           new_lease:   u64,
                           ri_hists: &RIHists) -> u64 {
        //if a set in a phase does not have this reference
        if ri_hists.ri_hists.get(&ref_id)==None {
            return 0;
        }
        let ref_ri_hist : &HashMap<u64,(u64,HashMap<u64,(u64,u64)>)> = 
            ri_hists.ri_hists.get(&ref_id).unwrap(); 
        let ri_hist: Vec<(u64,u64)> = ref_ri_hist.iter().map(|(k,v)|(*k,v.0)).collect();
        let mut old_cost = 0;
        let mut new_cost = 0;
        if phase != (ref_id & 0xFF000000) >> 24 {
            return 0;
        }
        for (ri,count) in ri_hist.iter(){
            if *ri <= old_lease {
                old_cost += *count * *ri;
            }
            else {
                old_cost += *count * old_lease;
            }

            if *ri <= new_lease {
                new_cost += *count * *ri;
            }
            else {
                new_cost += *count * new_lease;
            }
            
        }
        (new_cost - old_cost ) * sample_rate
    }

    pub fn get_ppuc(ref_id: u64, 
                base_lease: u64, 
                ref_ri_hist: &HashMap<u64,(u64,HashMap<u64,(u64,u64)>)>) -> Vec<PPUC>{

        let ri_hist: Vec<(u64,u64)> = ref_ri_hist.iter().map(|(k,v)|(*k,v.0)).collect();
        let total_count = ri_hist.iter().fold(0,|acc,(_k,v)| acc+v);
        let mut hits = 0;
        let mut head_cost = 0;

        let mut lease_hit_table = HashMap::new();
        let mut lease_cost_table = HashMap::new();

        //prevent kernel panic for a base lease that doesn't correspond to sampled ri for a reference
        lease_hit_table.insert(base_lease,0);
        lease_cost_table.insert(base_lease,0);

        let mut ri_hist_clone = ri_hist.clone();
        ri_hist_clone.sort_by(|a,b| a.0.cmp(&b.0));
        
        for (ri,count) in ri_hist_clone.iter(){
            hits += *count;
            head_cost += *count * *ri ;
            let tail_cost = (total_count - hits) * *ri ;

            lease_hit_table.insert(*ri,hits);
            lease_cost_table.insert(*ri,head_cost + tail_cost);
        }

        ri_hist_clone.iter().map(|(k,_v)| k).filter(|k| **k > base_lease).map(
            |k| PPUC {
                ppuc:((*lease_hit_table.get(k).unwrap() -*lease_hit_table.get(&base_lease).unwrap()) as f64/
                      (*lease_cost_table.get(k).unwrap()-*lease_cost_table.get(&base_lease).unwrap())as f64),
                lease: *k,
                old_lease: base_lease,
                ref_id: ref_id, 
                new_hits: *lease_hit_table.get(k).unwrap() - *lease_hit_table.get(&base_lease).unwrap(),
            }
        ).collect()

    }
     pub fn get_avg_lease(distribution:&BinnedRIs
                          , addr: &u64, bin: u64,lease: u64) ->u64{ 
        let mut total=0;
        for (ri,freq) in distribution.bin_ri_distribution.get(&bin).unwrap().get(&addr).unwrap(){
            if *ri<=lease && *ri>0{
                total+=ri*freq;
            }
            else {
                total+=lease*freq;
            }
        }
        return total;

     }
 pub fn prl(bin_width : u64,
                ri_hists : &RIHists,
                 binned_ris: &BinnedRIs,
                binned_freqs: &BinFreqs,
                sample_rate : u64,
                cache_size : u64,
                discretize: u64,
                samples_per_phase : &HashMap<u64,u64>,
                verbose: bool,
                debug: bool,
                set_mask : u32)-> Option<(HashMap<u64,u64>,HashMap<u64,(f64,u64)>,HashMap<u64,HashMap<u64,u64>>,u64)>{
        let mut new_lease: PPUC;
        let mut dual_leases : HashMap<u64,(f64,u64)>= HashMap::new(); //{ref_id, (alpha, long_lease)}
        let mut trace_length : u64=0;
        let mut bin_endpoints:Vec<u64>=Vec::new();
         let mut num_unsuitable:u64;
        let mut ppuc_tree = BinaryHeap::new();
       let mut impact_dict:HashMap<u64,HashMap<u64,f64>>=HashMap::new();
        let mut bin_ranks:HashMap<u64,f64>=HashMap::new();
        let mut sorted_bins:Vec<(u64,f64)>=Vec::new();
        let mut acceptable_ratio:f64;
        let mut neg_impact;
        let mut lease_hits=HashMap::new();
        let mut num_full_bins;
        let mut bin_saturation:HashMap<u64,HashMap<u64,f64>>=HashMap::new();
         let mut leases:HashMap<u64,u64>=HashMap::new();

        let num_sets=set_mask as u64+1;
         let bin_target:u64=bin_width*cache_size/num_sets;
       let min_alpha=1.0-(((2<<(discretize-1)) as f64)-1.5 as f64)/(((2<<(discretize-1)) as f64)-1.0 as f64); //threshold for meaningful dual lease

        if verbose{
        println!("bin_width:  {}",bin_width);
         }
       for key in binned_freqs.bin_freqs.keys(){
            bin_endpoints.push(*key);
       }
       //each bin will have all addresses although freq may be 0
       let mut addrs:Vec<u64>=Vec::new();
       for key in binned_freqs.bin_freqs.get(&0).unwrap().keys(){
            addrs.push(*key);
       }
     

        for endpoint in bin_endpoints.iter(){
            for set in 0..num_sets{
            bin_saturation.entry(*endpoint).or_insert(HashMap::new()).entry(set).or_insert(0.0);
            }
        }
        //make all references have lease of 1
        for addr in addrs{
            leases.insert(addr&0x00FFFFFF,1);
           // update saturation to take into account each reference having a lease of 1
            for (bin,_sat) in &bin_saturation.clone(){
                for set in 0..num_sets{
                    if  binned_ris.bin_ri_distribution.get(bin).unwrap().contains_key(&addr){
                        let  old_avg_lease=get_avg_lease(binned_ris,&addr,*bin,0);
                        let avg_lease =get_avg_lease(binned_ris,&addr,*bin,1);
                        let impact= (avg_lease as f64-old_avg_lease as f64)*&(sample_rate as f64);
                        bin_saturation.get_mut(&bin).unwrap().insert(set,impact);
                    }
                     //init impact dict for later
                        impact_dict.entry(*bin).or_insert(HashMap::new()).entry(set).or_insert(0.0);
                }
            }
       }
        
       
        for (_phase,&num) in samples_per_phase.iter(){
            trace_length += num * sample_rate;
        }
       
        for (&ref_id, ri_hist) in ri_hists.ri_hists.iter(){
            let ppuc_vec = get_ppuc(ref_id,0,ri_hist);
            for ppuc in ppuc_vec.iter(){
                ppuc_tree.push(*ppuc);
            }
        }
       // get lease hits assuming a base lease of 0
        for _r in ppuc_tree.clone(){
            let lease= ppuc_tree.pop().unwrap();
            lease_hits.entry(lease.ref_id&0x00FFFFFF).or_insert(HashMap::new()).entry(lease.lease).or_insert(lease.new_hits);
      }
        //reinitallize ppuc tree, assuming a base lease of 1
        for (&ref_id, ri_hist) in ri_hists.ri_hists.iter(){
            let ppuc_vec = get_ppuc(ref_id,1,ri_hist);
            for ppuc in ppuc_vec.iter(){
                ppuc_tree.push(*ppuc);
            }
        }
        
    
        loop {
             new_lease = match ppuc_tree.pop(){
                Some(i) => i,
                None => return Some((leases,dual_leases, lease_hits,trace_length)),
            };
           
             //continue to pop until we have a ppuc with the right base_lease
            if new_lease.old_lease != *leases.get(&(new_lease.ref_id&0xFFFFFFFF)).unwrap(){
                continue;
            }
            //won't assign a reference to a reference that has recieved a dual lease
            if dual_leases.contains_key(&(new_lease.ref_id&0xFFFFFFFF)){
                continue;
            }
            neg_impact=false;
            num_unsuitable=0;
            let addr= new_lease.ref_id;
           for (bin,bin_sat_set) in &bin_saturation{
            for (set,_sat) in bin_sat_set{
                let mut impact:f64 =0.0;
                let set_addr=(addr&0xFFFFFFFF)|(set<<32);
                if  binned_ris.bin_ri_distribution.get(&bin).unwrap().contains_key(&set_addr){
                    let  old_avg_lease=get_avg_lease(binned_ris,&set_addr,*bin,*leases.get(&(addr&0xFFFFFFFF)).unwrap());
                    let avg_lease =get_avg_lease(binned_ris,&set_addr,*bin,new_lease.lease);
                    impact= (avg_lease as f64-old_avg_lease as f64)*&(sample_rate as f64);
                    //don't assign leases that decrease bin saturation
                    neg_impact = if impact>=0.0 {false} else {true};
                    impact_dict.get_mut(&bin).unwrap().insert(*set,impact as f64);
                }
                else{
                    impact_dict.get_mut(&bin).unwrap().insert(*set,0 as f64);
                }

                if (bin_saturation.get(&bin).unwrap().get(&set).unwrap()+impact)>bin_target as f64{
                   
                    num_unsuitable+=1;
                }
            }
        }
            for (bin,sat_set) in &bin_saturation.clone(){
                for (set,sat) in sat_set{
                    if verbose && debug {
                        println!("bin: {} set: {} current capacity: {:.7} impact: {:.7}", bin/bin_width, set, sat/bin_width as f64,impact_dict.get(bin).unwrap().get(set).unwrap()/bin_width as f64); 
                    }
                }
             }
             if verbose{
                println!("addr:{:x} ri:{:x}",addr&0xFFFFFFFF,new_lease.lease);
             }
            
            //skip lease, if it makes it worse
           
            if neg_impact{
                if verbose{
                    println!("lease value would be negative");
                }
                continue;
            }
            if num_unsuitable<1{
                leases.insert(addr&0xFFFFFFFF,new_lease.lease);
                 //push new ppucs
                let ppuc_vec = get_ppuc(new_lease.ref_id,
                                        new_lease.lease,
                                        ri_hists.ri_hists.get(&new_lease.ref_id).unwrap());

                for ppuc in ppuc_vec.iter(){
                    ppuc_tree.push(*ppuc);
                }
                 let mut print_string:String=String::new();
                for (bin,sat_set) in &bin_saturation.clone(){
                    for (set,sat) in sat_set{ 
                        let set_addr=(addr&0xFFFFFFFF)|(set<<32);
                     if  binned_ris.bin_ri_distribution.get(bin).unwrap().contains_key(&set_addr){
                        bin_saturation.get_mut(bin).unwrap().insert(*set,sat+impact_dict.get(bin).unwrap().get(set).unwrap());
                    }
                    print_string=format!("{:} {1:.7}",print_string,&(bin_saturation.get(bin).unwrap().get(set).unwrap()/bin_width as f64));
                   }
                }
                if verbose{
                  println!("assigning lease: {:x} to reference {:x}",new_lease.lease, addr);
                 println!("Average cache occupancy per bin: [{}]",print_string);
             }
            }
            else {
                num_full_bins=0;
                let mut alpha=1.0;
                for (bin,sat_set) in &bin_saturation{
                    for (set,sat) in sat_set{
                        if sat>=&(bin_target as f64){
                            num_full_bins+=1
                        }
                        let new_capacity=sat+impact_dict.get(bin).unwrap().get(set).unwrap();
                        
                        if&new_capacity>=&(bin_target as f64){
                            if *impact_dict.get(bin).unwrap().get(set).unwrap()!=0.0{
                                //get minimum set alpha for bin
                                let set_alpha = ((bin_target as f64) -sat)/ *impact_dict.get(bin).unwrap().get(set).unwrap();
                                if set_alpha<alpha {
                                    alpha=set_alpha;
                                }
                            }
                          
                        }
                    }
                      bin_ranks.insert(*bin,alpha);
                }
                for (num,key_val_pair) in bin_ranks.iter().enumerate(){
                    if sorted_bins.get_mut(num)==None{
                        sorted_bins.push((*key_val_pair.0,*key_val_pair.1));
                    }
                    else{
                        sorted_bins[num]=(*key_val_pair.0,*key_val_pair.1);
                    }
                }
                sorted_bins.sort_by(|a,b| a.1.partial_cmp(&b.1).unwrap());
                acceptable_ratio= if num_full_bins==0 {sorted_bins[0].1} else {0.0};
               
                if acceptable_ratio>min_alpha{
                    dual_leases.insert(addr,(acceptable_ratio as f64,new_lease.lease));
                        let mut print_string:String=String::new();

                    for (bin,sat_set) in &bin_saturation.clone(){
                        for (set,sat) in sat_set{ 
                            let set_addr=(addr&0xFFFFFFFF)|(set<<32);
                            if  binned_ris.bin_ri_distribution.get(bin).unwrap().contains_key(&set_addr){
                              bin_saturation.get_mut(bin).unwrap().insert(*set,sat+impact_dict.get(bin).unwrap().get(set).unwrap()*acceptable_ratio);
                             }
                            print_string=format!("{:} {1:.7}",print_string,&(bin_saturation.get(bin).unwrap().get(set).unwrap()/bin_width as f64));
                        }
                    }
                    if verbose{
                            println!("Assigning dual lease {} to address {} with percentage: {}",new_lease.lease,addr,acceptable_ratio);
                            println!("Average cache occupancy per bin: [{}]",print_string);
                        }
                    }
                }
            }

}


    pub fn shel_cshel(cshel: bool,
                  ri_hists : &RIHists, 
                  cache_size : u64, 
                  sample_rate : u64, 
                  samples_per_phase : &HashMap<u64,u64>,
                  discretize : u64,
                  verbose: bool,
                  debug: bool,
                  set_mask: u32) -> Option<(HashMap<u64,u64>,HashMap<u64,(f64,u64)>,HashMap<u64,HashMap<u64,u64>>,u64)> {

        let mut new_lease: PPUC;
        let mut cost_per_phase:HashMap<u64,HashMap<u64,u64>>=HashMap::new();
        let mut budget_per_phase:HashMap<u64,u64>=HashMap::new();
        let mut leases = HashMap::new(); //{ri, lease}
        let mut dual_leases : HashMap<u64,(f64,u64)>= HashMap::new(); //{ref_id, (alpha, long_lease)}
        let mut trace_length : u64 = 0;
        let mut lease_hits=HashMap::new();
        let mut dual_lease_phases :Vec<u64>=Vec::new();
        //{phase,(cost with alpha, cost if alpha was 1, ref ID)}
        let mut past_lease_values :HashMap<u64,(u64,u64)>=HashMap::new();
        let mut last_lease_cost: HashMap<u64,HashMap<u64,(u64,u64,u64)>>=HashMap::new();
        
        let num_sets=set_mask as u64+1;
        let phase_ids: Vec<&u64> = samples_per_phase.keys().collect();
        //since we can't run CSHEL without also running SHEL, don't output RI history twice
        if !cshel {
            if verbose {
            println!("---------Dump RI Hists------------");
            super::io::debug::print_ri_hists(&ri_hists);
            println!("---------Dump Samples Per Phase---");
            println!("{:?}",&samples_per_phase);
            }
        }
        
        
        let min_alpha:f64=1.0-(((2<<(discretize-1)) as f64)-1.5 as f64)/(((2<<(discretize-1)) as f64)-1.0 as f64); //threshold for meaningful dual lease
         //initialize ppucs
        let mut ppuc_tree = BinaryHeap::new();

        for (&ref_id, ri_hist) in ri_hists.ri_hists.iter(){
           
            let ppuc_vec = get_ppuc(ref_id,0,ri_hist);
            for ppuc in ppuc_vec.iter(){
                ppuc_tree.push(*ppuc);
            }
        }
       // get lease hits assuming a base lease of 0
        for _r in ppuc_tree.clone(){
            let lease= ppuc_tree.pop().unwrap();
            //sum hits for reference over all sets
            *lease_hits.entry(lease.ref_id&0x00FFFFFF).or_insert(HashMap::new()).entry(lease.lease).or_insert(0)+=lease.new_hits;

      }
        //reinitallize ppuc tree, assuming a base lease of 1
        for (&ref_id, ri_hist) in ri_hists.ri_hists.iter(){
            let ppuc_vec = get_ppuc(ref_id,1,ri_hist);
            for ppuc in ppuc_vec.iter(){
                ppuc_tree.push(*ppuc);
            }
        }

        //initialize cost + budget
        for (&phase,&num) in samples_per_phase.iter(){
            
                budget_per_phase.entry(phase).or_insert(num * cache_size / num_sets * sample_rate);
             trace_length += num * sample_rate;
        }

        if verbose{
            println!("
            ---------------------
            Initial budget per phase: 
            {:?} 
            ---------------------",
            budget_per_phase);
        }
        //initialize leases to a default value of 1
         for (&ref_id, _ ) in ri_hists.ri_hists.iter(){
            leases.insert(ref_id&0xFFFFFFFF,1);
              let phase=(ref_id&0xFF000000)>>24;
              let phase_id_ref=ref_id&0xFFFFFFFF;
            // get cost of assigning a lease of 1 for each set
            for set in 0..num_sets{
                let set_phase_id_ref=phase_id_ref|(set<<32);
                 let new_cost = match cshel{
                          true => cshel_phase_ref_cost(sample_rate,
                                                            phase,
                                                            set_phase_id_ref,
                                                            0,
                                                            1,
                                                            &ri_hists),
                          false => shel_phase_ref_cost(sample_rate,
                                                            phase,
                                                            set_phase_id_ref,
                                                            0,
                                                            1,
                                                            &ri_hists),
                  };
                  if cost_per_phase.get(&phase)==None{
                    cost_per_phase.insert(phase,HashMap::new());
                  }
                  else {
                   cost_per_phase.get_mut(&phase).unwrap().insert(set,new_cost);
                 }
             }
            
         }
         if verbose {
             println!("costs per phase{:?}",cost_per_phase);
        }

        loop {
            new_lease = match ppuc_tree.pop(){
                Some(i) => i,
                None => return Some((leases,dual_leases,lease_hits,trace_length)),
            };
            let phase=(new_lease.ref_id&0xFFFFFFFF)>>24;
            let ref_id=new_lease.ref_id&0xFFFFFFFF;
           
             //continue to pop until we have a ppuc with the right base_lease
            if new_lease.old_lease != *leases.get(&ref_id).unwrap(){
                continue;
            }
            
             let mut set_full=false;
                for set in 0..num_sets{
                    if cost_per_phase.get(&phase).unwrap().get(&set).unwrap()==budget_per_phase.get(&phase).unwrap(){
                        set_full=true;
                        break;
                    }
                }
            //if any set in phase is full, skip
            // if   dual_lease_phases.contains(&(new_lease.ref_id>>24)){  
                if set_full{
                    continue;
                }
                //if we've already assigned dual leases to all phases, end
            if dual_lease_phases.len()==cost_per_phase.len(){
              return Some((leases,dual_leases,lease_hits,trace_length)) ;
            }
            //if we've already assigned a dual lease for the phase
            if dual_lease_phases.contains(&phase){
                continue;
            }

            let old_lease = *leases.get(&ref_id).unwrap();
            //check for capacity
            let mut acceptable_lease = true;
            let mut new_phase_ref_cost :HashMap<u64,HashMap<u64,u64>>=HashMap::new(); 
            for (&phase,current_cost) in cost_per_phase.iter(){
                
            // get cost of assigning a lease of 1 for each set
                for set in 0..num_sets{
                    let set_phase_id_ref=ref_id|(set<<32);
                    let new_cost = match cshel{
                        true => cshel_phase_ref_cost(sample_rate,
                                                          phase,
                                                          set_phase_id_ref,
                                                          old_lease,
                                                          new_lease.lease,
                                                          &ri_hists),
                        false => shel_phase_ref_cost(sample_rate,
                                                          phase,
                                                          set_phase_id_ref,
                                                          old_lease,
                                                          new_lease.lease,
                                                          &ri_hists),
                    };
                  
                     new_phase_ref_cost.entry(phase).or_insert(HashMap::new()).entry(set).or_insert(new_cost); 
                    if (new_cost + current_cost.get(&set).unwrap()) > *budget_per_phase.get(&phase).unwrap() {
                        acceptable_lease = false;
                    }
                }
            }
            if verbose & debug {
                 println!("\nDebug: budgets per phase {:?}",&budget_per_phase);
                println!("Debug: Current cost budgets {:?}",&cost_per_phase);
                println!("Debug: NEW_PHASE_REF_COST {:?}",&new_phase_ref_cost);
               
             }

            if acceptable_lease {
           
           
                //update cache use
                for (phase,phase_set_costs) in cost_per_phase.iter_mut() {
                    for (set,set_costs) in phase_set_costs.iter_mut(){
                            *set_costs+=new_phase_ref_cost.get(phase).unwrap().get( set).unwrap();
                    }
                }
                let phase=(new_lease.ref_id & 0xFF000000)>>24;
                //store lease value we assign to the reference and the value of the previously assigned lease for that reference
                past_lease_values.insert(new_lease.ref_id&0xFFFFFFFF,(new_lease.lease,*leases.get(&(&new_lease.ref_id&0xFFFFFFFF)).unwrap()));
                
                if last_lease_cost.get_mut(&phase)==None{
                    last_lease_cost.insert(phase,HashMap::new());
                }
                for set in 0..num_sets{
                    last_lease_cost.get_mut(&phase).unwrap().insert(set,(*new_phase_ref_cost.get(&phase).unwrap().get(&set).unwrap(),*new_phase_ref_cost.get(&phase).unwrap().get(&set).unwrap(),ref_id&0xFFFFFFFF));
                }  
                //update leases
                leases.insert(new_lease.ref_id&0xFFFFFFFF,new_lease.lease);
                //push new ppucs
                let ppuc_vec = get_ppuc(new_lease.ref_id,
                                        new_lease.lease,
                                        ri_hists.ri_hists.get(&new_lease.ref_id).unwrap());

                for ppuc in ppuc_vec.iter(){
                    ppuc_tree.push(*ppuc);
                }
                if verbose {
                    println!("Assigned lease {:x} to reference ({},{:x}).", 
                             new_lease.lease, (new_lease.ref_id & 0xFF000000) >> 24, 
                             new_lease.ref_id & 0x00FFFFFF);
                }
            }

            else {
                //unacceptable lease, must assign a dual lease
                let mut alpha = 1.0;
                let mut current_phase_alpha = 1.0;
                for (&phase,phase_set_current_cost) in cost_per_phase.iter(){
                     let set_budget=*budget_per_phase.get(&phase).unwrap();
                    for(&set,&current_set_cost) in phase_set_current_cost.iter(){
                        let &set_phase_ref_cost   = new_phase_ref_cost.get(&phase).unwrap().get(&set).unwrap(); 
                        if set_phase_ref_cost > 0 {
                            if set_budget< current_set_cost{
                                println!("
                                ERROR: current cost exceeds budget
                                *budget_per_phase.get(&phase)=.unwrap():  {}
                                currenc_cost:                            {}
                                ",
                                set_budget,
                                current_set_cost
                                );
                                panic!();
                            }

                            let remaining_budget = set_budget - current_set_cost; 
                            //get the best alpha for any set  (ignoring other phases) that we want for the current reference
                            if phase==(new_lease.ref_id&0xFF000000)>>24 {
                                current_phase_alpha=super::helpers::float_min(current_phase_alpha,remaining_budget as f64/ 
                                                                                                   set_phase_ref_cost as f64);
                            }
                            alpha = super::helpers::float_min(alpha, 
                                                              remaining_budget as f64 / 
                                                                set_phase_ref_cost as f64);
                        }
                    }

                }
                //if the alpha we wish to assign would result in a long lease that is never used because the short lease probabiliy will be 1 after descretizing
                //don't assign dual lease.
                if current_phase_alpha<min_alpha{
                     println!("Assigning lease {:x} with percentage {} to reference ({},{:x}) would not be meaningful.", 
                                 new_lease.lease,current_phase_alpha,(new_lease.ref_id & 0xFF000000) >> 24, 
                                 new_lease.ref_id & 0x00FFFFFF);
                    continue;
                }
               
                if alpha > min_alpha{
                     //update cache use
                    for (phase,phase_set_costs) in cost_per_phase.iter_mut() {
                        let mut set_budget=*budget_per_phase.get(phase).unwrap();
                        for (set,set_costs) in phase_set_costs.iter_mut(){
                                *set_costs+=(*new_phase_ref_cost.get(phase).unwrap().get( set).unwrap() as f64 *alpha).round() as u64;
                                //fix floating point precision error leading to "overallocation" or underallocation
                                if set_costs>&mut set_budget{
                                    *set_costs=set_budget;
                                }
                        }

                    }
                   
                }
                
              if cshel{
                    //if there's no alpha that would assign a meaningful dual lease that wouldn't put other phases over budget
                    if alpha <=min_alpha{
                        let mut new_costs=HashMap::new();
                        let mut new_alpha=HashMap::new();
                        let mut adjust_lease=true;
                        let mut phase_alpha=1.0;
                        for phase in &phase_ids{
                            for set in 0..num_sets{
                             let set_phase_ref_cost   = new_phase_ref_cost.get(&phase).unwrap().get(&set).unwrap(); 
                             //if the phase would be effected by the lease assignment
                                if set_phase_ref_cost > &0 {
                                     //get phases that would be over budgeted by assigning the current lease.
                                     //then subtract the cost of their prior dual lease (which may be, due to the default, a non-dual lease)
                                    //and then add the spillover cost from the new leases
                                    let past_cost_actual=if last_lease_cost.get(&phase)==None {0 } else if last_lease_cost.get(&phase).unwrap().get(&set)==None{0}
                                    else {last_lease_cost.get(&phase).unwrap().get(&set).unwrap().0};
                                    if new_costs.get(&phase)==None{
                                        new_costs.insert(phase,HashMap::new());
                                    }
                                    let new_cost=cost_per_phase.get(&phase).unwrap().get(&set).unwrap()-past_cost_actual+(*set_phase_ref_cost as f64 * current_phase_alpha).round() as u64;
                                new_costs.get_mut(&phase).unwrap().insert(set,new_cost); 
                                  //if no lease adjustment can be made to keep the phase from being over budget 
                                    if new_costs.get(&phase).unwrap().get(&set).unwrap()>budget_per_phase.get(&phase).unwrap(){
                                        adjust_lease=false;
                                        break;
                                    }
                                     let remaining_budget = *budget_per_phase.get(&phase).unwrap()  - new_costs.get(&phase).unwrap().get(&set).unwrap(); 
                                           //if cost of last lease was zero i.e., no prior lease for phase, then alpha will be 1 and will not be adjusted
                                    let past_cost_max=if past_cost_actual!=0 {last_lease_cost.get(phase).unwrap().get(&set).unwrap().1} else {0};
                                     if past_cost_max!=0{
                                         //if previous long lease didn't fill phase, could be greater than one
                                        let set_phase_alpha=super::helpers::float_min(1.0,remaining_budget as f64/past_cost_max as f64);
                                        if set_phase_alpha<=min_alpha{
                                            let old_phase_ref=last_lease_cost.get(phase).unwrap().get(&set).unwrap().2;
                                            dual_leases.get(&old_phase_ref).unwrap().1;
                                              println!("Assigning adjusted dual lease {:x} with percentage {} to reference ({},{:x}) would not be meaningful.", 
                                        new_lease.lease,set_phase_alpha,phase,old_phase_ref);
                                               adjust_lease=false;
                                              break;
                                        }
                                        //need the minimum alpha of any set in the phase
                                        else if set_phase_alpha<phase_alpha{
                                            phase_alpha=set_phase_alpha;
                                        }

                                        new_alpha.insert(phase,phase_alpha);

                                    }
                                     
                                    
                                }
                                //new costs is equal to old cost
                                else{
                                    if new_costs.get(&phase)==None{
                                        new_costs.insert(&phase,HashMap::new());
                                    }
                                    new_costs.get_mut(&phase).unwrap().insert(set,*cost_per_phase.get(&phase).unwrap().get(&set).unwrap());
                                }
                            }
                        }
                        if adjust_lease==true {
                            for phase in &phase_ids{
                                //if adjusting lease
                                for set in 0..num_sets{
                                    if new_alpha.get(&phase)!=None{
                                        let old_phase_cost_max=last_lease_cost.get(phase).unwrap().get(&set).unwrap().1;
                                        let old_phase_ref=last_lease_cost.get(phase).unwrap().get(&set).unwrap().2;
                                        let new_phase_cost=(old_phase_cost_max as f64*new_alpha.get(phase).unwrap()) as u64;
                                       
                                        //if phase had a dual lease
                                        if dual_lease_phases.contains(phase){
                                           dual_leases.insert(old_phase_ref,(*new_alpha.get(&phase).unwrap(),dual_leases.get(&old_phase_ref).unwrap().1));
                                        }
                                        //if we are not currently assigning a dual lease to this phase
                                        //generate dual lease from the past two lease values of the last reference assigned in this phase
                                        else if **phase!=new_lease.ref_id>>24{
                                            //set prior single lease as long lease value with new alpha
                                            dual_leases.insert(old_phase_ref,(*new_alpha.get(&phase).unwrap(),past_lease_values.get(&old_phase_ref).unwrap().0));
                                            //set the lease two references back as the short lease value
                                            leases.insert(old_phase_ref,past_lease_values.get(&old_phase_ref).unwrap().1);
                                            dual_lease_phases.push(**phase);
                                        } 

                                        last_lease_cost.get_mut(&phase).unwrap().insert(set,(new_phase_cost,old_phase_cost_max,old_phase_ref));
                                        //update phase costs
                                        cost_per_phase.get_mut(&phase).unwrap().insert(set,*new_costs.get(&phase).unwrap().get(&set).unwrap()+new_phase_cost);
                                    }
                                    //if not adjusting the lease
                                    else {
                                            //update phase costs
                                        cost_per_phase.get_mut(&phase).unwrap().insert(set,*new_costs.get(&phase).unwrap().get(&set).unwrap());
                                       //fix floating point precision error leading to "overallocation"
                                        if cost_per_phase.get(&phase).unwrap().get(&set).unwrap()>budget_per_phase.get(&phase).unwrap(){
                                            cost_per_phase.get_mut(&phase).unwrap().insert(set,*budget_per_phase.get(&phase).unwrap());
                                        }

                                    }
                                }
                            }
                            alpha=current_phase_alpha;
                        }
                        else {
                            //if we can't assign a dual lease without overflowing a phase
                            //without adjustment of past dual leases, with adjustment of past dual leases, 
                            //or in the the unlikely case a phase is full with no dual lease
                           
                              println!("Unable to assign lease {:x} with percentage {} to reference ({},{:x})", 
                                 new_lease.lease,current_phase_alpha,(new_lease.ref_id & 0xFF000000) >> 24, 
                                 new_lease.ref_id & 0x00FFFFFF);
                               continue; 
                        }

                    }
                }
              

                
            
                let phase=(new_lease.ref_id & 0xFF000000)>>24;

                //detect if set full
                let mut set_full=false;
                for set in 0..num_sets{
                    if cost_per_phase.get(&phase).unwrap().get(&set).unwrap()==budget_per_phase.get(&phase).unwrap(){
                        set_full=true;
                        break;
                    }
                }
                //if last lease was a dual lease with alpha of 1 that didn't fill the budget, then it is actually a short lease and adjustments can be made to ensure 
                //there is only 1 dual lease per phase.
                if alpha==1.0 && set_full==false{
                    //update leases
                    leases.insert(new_lease.ref_id&0xFFFFFFFF,new_lease.lease);
                    //push new ppucs
                    let ppuc_vec = get_ppuc(new_lease.ref_id,
                                            new_lease.lease,
                                            ri_hists.ri_hists.get(&new_lease.ref_id).unwrap());

                    for ppuc in ppuc_vec.iter(){
                        ppuc_tree.push(*ppuc);
                    }
                    if verbose {
                        println!("Assigned lease {:x} to reference ({},{:x}).", 
                                 new_lease.lease, (new_lease.ref_id & 0xFF000000) >> 24, 
                                 new_lease.ref_id & 0x00FFFFFF);
                    }

                }
                //add dual lease
                else{

                    //store cost of dual lease and store cost of lease with no dual lease and the reference for that lease
                     for set in 0..num_sets{
                      last_lease_cost.get_mut(&phase).unwrap().insert(set,((*new_phase_ref_cost.get(&phase).unwrap().get(&set).unwrap() as f64 *alpha).round() as u64,
                        *new_phase_ref_cost.get(&phase).unwrap().get(&set).unwrap(),ref_id&0xFFFFFFFF));
                    }  
                    
                   dual_lease_phases.push(phase);
                   //update dual lease HashMap
                    dual_leases.insert(new_lease.ref_id&0xFFFFFFFF,(alpha,new_lease.lease));
                    

                    
                    if verbose {
                        println!("Assigned dual lease ({:x},{}) to reference ({},{:x}).", 
                                  new_lease.lease, 
                                  alpha, 
                                  (new_lease.ref_id & 0xFF000000) >> 24,
                                  new_lease.ref_id & 0x00FFFFFF);
                    }
                   
                }
              
            }//unacceptable lease
       
            if verbose & debug { 
                
                for (phase,num) in samples_per_phase.iter(){
                    for set in 0..num_sets{


                    println!("Debug: phase: {}. set: {} Allocation: {}",phase,set,
                        cost_per_phase.get(&phase).unwrap().get(&set).unwrap()  / (num*sample_rate));
                    }
                }
                /*
                    println!("Debug: phase: {}",phase);
                    println!("Debug:    cost_per_phase:   {:?}",
                             cost_per_phase.get(&phase).unwrap());
                    println!("Debug:    budget_per_phase: {:?}",
                             budget_per_phase.get(&phase).unwrap());
                }*/
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use super::io::*;
    use super::io::debug::*;
    use super::helpers::*;
    use super::lease_gen::*;

    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }

    #[test]
    fn test_bin_search(){
        let a = vec![(0,0),(5,1),(120,2),(288,0),(1025,1)];
        assert_eq!(binary_search(&a,2),Some((5,1)));
        assert_eq!(binary_search(&a,120),Some((288,0)));
        assert_eq!(binary_search(&a,1024),Some((1025,1)));
        assert_eq!(binary_search(&a,2048),None);
        assert_eq!(binary_search(&a,1025),None);
    }

    #[test]
    fn test_process_sample_head_cost(){
        let mut ri_hists = HashMap::new();
        process_sample_head_cost(&mut ri_hists,1,12,20,(100000,1));
        let hist_struct = RIHists::new(ri_hists);

        assert_eq!(hist_struct.get_ref_ri_count(1,12),1);
        assert_eq!(hist_struct.get_ref_ri_phase_cost(1,12,0),(12,0));

    }

    #[test]
    fn test_cross_phase_head_cost(){

        let mut ri_hists = HashMap::new();
        process_sample_head_cost(&mut ri_hists,1,12,8,(10,1)); //reference 1, phase 0
        let hist_struct = RIHists::new(ri_hists);

        assert_eq!(hist_struct.get_ref_ri_count(1,12),1);
        assert_eq!(hist_struct.get_ref_ri_phase_cost(1,12,0),(2,0));
        assert_eq!(hist_struct.get_ref_ri_phase_cost(1,12,1),(10,0));
    }

    #[test]
    fn test_get_ppuc(){
        let mut ri_hist = HashMap::new();

        ri_hist.insert(3,(1,HashMap::new()));
        ri_hist.insert(5,(7,HashMap::new()));
        ri_hist.insert(17,(4,HashMap::new()));
        ri_hist.insert(19,(3,HashMap::new()));
        let ppucs = get_ppuc(1,0,&ri_hist);
        println!("{:?}",ppucs);//According to Asplos19, the numbers are [.02,.11,.08,.09];
    }


    #[test]
    fn test_process_sample_tail_cost(){
        let mut ri_hists = HashMap::new();
        let ri_short = 10;
        let ri_long = 100;
        let ri_very_long = 1000;
        process_sample_head_cost(&mut ri_hists,1,ri_short,10,(10000,1));
        process_sample_head_cost(&mut ri_hists,1,ri_short,20,(10000,1));
        process_sample_head_cost(&mut ri_hists,1,ri_long,200,(10000,1));
        process_sample_head_cost(&mut ri_hists,1,ri_very_long,2200,(10000,1));

        process_sample_tail_cost(&mut ri_hists,1,ri_short,10,(10000,1));
        process_sample_tail_cost(&mut ri_hists,1,ri_short,20,(10000,1));
        process_sample_tail_cost(&mut ri_hists,1,ri_long,200,(10000,1));
        process_sample_tail_cost(&mut ri_hists,1,ri_very_long,2200,(10000,1));

        let hist_struct = RIHists::new(ri_hists);

       // print_ri_hists(&hist_struct);
        assert_eq!(hist_struct.get_ref_ri_phase_cost(1,ri_short,0).1,20);
        assert_eq!(hist_struct.get_ref_ri_phase_cost(1,ri_long,0).1,100);
        assert_eq!(hist_struct.get_ref_ri_phase_cost(1,ri_very_long,0).1,0);

    }

    //pub fn process_sample_head_cost(ri_hists: &mut HashMap<u64,HashMap<u64,(u64,HashMap<u64,(u64,u64)>)>>,
     //p                 phase_id_ref: u64,
        //p              ri: u64,
           //p           use_time: u64,
              //p        next_phase_tuple: (u64,u64)){
    #[test]
    fn tail_cost_cross_phase(){
        let mut ri_hists = HashMap::new();
        process_sample_head_cost(&mut ri_hists,1,100,70,(100,1));
        process_sample_head_cost(&mut ri_hists,1,50,70,(100,1));

        process_sample_tail_cost(&mut ri_hists,1,100,70,(100,1));
        process_sample_tail_cost(&mut ri_hists,1,50,70,(100,1));

        let hist_struct = RIHists::new(ri_hists);

        print_ri_hists(&hist_struct);
        assert_eq!(hist_struct.get_ref_ri_phase_cost(1,50,0).1,30);
        assert_eq!(hist_struct.get_ref_ri_phase_cost(1,50,1).1,20);
    }

    #[test]
    fn negative_ri(){
        let mut ri_hists = HashMap::new();
        process_sample_head_cost(&mut ri_hists,1, 50, 80, (100,1));
        process_sample_head_cost(&mut ri_hists,1,i32::max as u64, 90, (100,1));

        process_sample_tail_cost(&mut ri_hists,1, 50, 80, (100,1));
        process_sample_tail_cost(&mut ri_hists,1,i32::max as u64, 90, (100,1));
        let hist_struct = RIHists::new(ri_hists);

        assert_eq!(hist_struct.get_ref_ri_phase_cost(1,50,0).0,20);
        assert_eq!(hist_struct.get_ref_ri_phase_cost(1,50,1).0,30);
        assert_eq!(hist_struct.get_ref_ri_phase_cost(1,50,0).1,10);
        assert_eq!(hist_struct.get_ref_ri_phase_cost(1,50,1).1,40);

    }
}
