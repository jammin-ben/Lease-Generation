//use std::collections::HashMap;
//use std::collections::BinaryHeap;
//use std::u64;


// Functions for parsing input files, debug prints, 
// and lease output
pub mod io {
    use serde::{Serialize, Deserialize};
    use std::collections::BinaryHeap;
    use std::collections::HashMap;


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
    pub fn build_phase_transitions(input_file:&str) -> Vec<(u64,u64)>{
        println!("Reading input from: {}", input_file);

        let mut rdr = csv::ReaderBuilder::new()
            .has_headers(false)
            .from_path(input_file)
            .unwrap();

        let mut sample_hash = HashMap::new();

        for result in rdr.deserialize() {
            let sample: Sample = result.unwrap();
            let ri = u64::from_str_radix(&sample.ri,16).unwrap();
            //if sample is negative, there is no reuse, so ignore
            if ri>2147483647 {
                continue;
            }
            let mut time = sample.time;
            let phase_id_ref = u64::from_str_radix(&sample.phase_id_ref,16).unwrap();

            let phase_id = (phase_id_ref & 0xFF000000)>>24;
            time = time - ri;
            sample_hash.insert(time,phase_id);
        }

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

        sorted_transitions.iter().map(|&x| (*(x.0),*(x.1))).collect()
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
    pub fn build_ri_hists(input_file:&str) -> (super::lease_gen::RIHists,HashMap<u64,u64>){
        let phase_transitions = build_phase_transitions(input_file);
        let mut rdr = csv::ReaderBuilder::new()
            .has_headers(false)
            .from_path(input_file)
            .unwrap();

        let mut ri_hists = HashMap::new();
        let mut samples_per_phase = HashMap::new();

        //first pass for head costs
        for result in rdr.deserialize() {
            let sample: Sample = result.unwrap();
            let ri = u64::from_str_radix(&sample.ri,16).unwrap();
            //if sample is negative, there is no reuse, so ignore
            if ri>2147483647 {
                continue;
            }
            let time = sample.time;
            let phase_id_ref = u64::from_str_radix(&sample.phase_id_ref,16).unwrap();
            let next_phase_tuple = match super::helpers::binary_search(&phase_transitions,time-ri){
                Some(v) => v,
                None => (time+1,0),
            };

            super::lease_gen::process_sample_head_cost(&mut ri_hists,phase_id_ref,ri,time,next_phase_tuple);

            let phase_id = (phase_id_ref & 0xFF000000)>>24;
            *samples_per_phase.entry(phase_id).or_insert_with(|| 0)+=1;
        }

        let mut rdr = csv::ReaderBuilder::new()
            .has_headers(false)
            .from_path(input_file)
            .unwrap();
        //second pass for tail costs
        for result in rdr.deserialize() {
            let sample: Sample = result.unwrap();
            let ri = u64::from_str_radix(&sample.ri,16).unwrap();
            //if sample is negative, there is no reuse, so ignore
            if ri>2147483647 {
                continue;
            }
            let time = sample.time;
            let phase_id_ref = u64::from_str_radix(&sample.phase_id_ref,16).unwrap();
            let next_phase_tuple = match super::helpers::binary_search(&phase_transitions,time-ri){
                Some(v) => v,
                None => (time+1,0),
            };
            super::lease_gen::process_sample_tail_cost(&mut ri_hists,
                                     phase_id_ref,
                                     ri,
                                     time,
                                     next_phase_tuple);
        }

        (super::lease_gen::RIHists::new(ri_hists),samples_per_phase)
    }
    pub fn dump_leases(leases: HashMap<u64,u64>, dual_leases: HashMap<u64,(f32,u64)>) {
        println!("Dump formated leases");
              let mut lease_vector: Vec<(u64,u64,u64,u64,f32)> = Vec::new();
        for (&phase_address,&lease) in leases.iter(){
            let lease = if lease >0 {lease} else {1}; 
            let phase   = (phase_address & 0xFF000000)>>24;
            let address =  phase_address & 0x00FFFFFF;
            if dual_leases.contains_key(&phase_address){
               lease_vector.push((phase,address,lease,dual_leases.get(&phase_address).unwrap().1,1.0-dual_leases.get(&phase_address).unwrap().0));
            }
            else{
                lease_vector.push((phase,address,lease,0, 1.0));
            }
        }
        lease_vector.sort_by_key(|a| a.0); //sort by phase
        lease_vector.sort_by_key(|a| a.1); //sort by reference
        for (phase, address, lease_short, lease_long, percentage) in lease_vector.iter(){
            println!("{:x}, {:x}, {:x}, {:x}, {}",phase, address, lease_short, lease_long, percentage);
        }
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
    pub fn float_min(a: f32, b:f32) -> f32{
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
        ppuc: f32,
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
                      time: u64,
                      next_phase_tuple: (u64,u64)){

        let phase_id = (phase_id_ref & 0xFF000000)>>24;
        let start_time = time - ri;

        //increment count
        let ref_hist = ri_hists.entry(phase_id_ref).or_insert_with(|| HashMap::new());
        let ri_tuple = ref_hist.entry(ri).or_insert_with(|| (0,HashMap::new()));
        ri_tuple.0  += 1;

        //increment head costs
        let this_phase_head_cost= std::cmp::min(next_phase_tuple.0 - start_time,ri);
        let next_phase_head_cost= std::cmp::max(time as i64 - next_phase_tuple.0 as i64,
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
                      time: u64,
                      next_phase_tuple: (u64,u64)){

        let phase_id = (phase_id_ref & 0xFF000000)>>24;
        let start_time = time - ri;

        let ref_hist = ri_hists.entry(phase_id_ref)
            .or_insert_with(|| HashMap::new());

        //this heinous code exists so we can iterate through a hashmap while modifying it
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

            let this_phase_tail_cost = std::cmp::min(next_phase_tuple.0 - start_time,
                                                     ri_other);
            let next_phase_tail_cost = std::cmp::max(0,
                                                     start_time as i64 + 
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
        let ri_hist = ri_hists.ri_hists.get(&ref_id).unwrap();

        for (&ri,(_,phase_cost_hashmap)) in ri_hist.iter(){
            let (phase_head_cost,phase_tail_cost) = match phase_cost_hashmap.get(&phase) {
                Some((a,b)) => (*a,*b),
                None        => (0,0), 
            };

            if ri <= old_lease {
                old_cost += phase_head_cost;
            }
            if ri == old_lease {
                old_cost += phase_tail_cost;
            }

            if ri <= new_lease {
                new_cost += phase_head_cost;
            }
            if ri == new_lease {
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

        if phase != (ref_id & 0xFF000000) >> 24 {
            return 0;
        }

        let ref_ri_hist : &HashMap<u64,(u64,HashMap<u64,(u64,u64)>)> = 
            ri_hists.ri_hists.get(&ref_id).unwrap(); 
        let ri_hist: Vec<(u64,u64)> = ref_ri_hist.iter().map(|(k,v)|(*k,v.0)).collect();
        let mut old_cost = 0;
        let mut new_cost = 0;

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

    fn get_ppuc(ref_id: u64, 
                base_lease: u64, 
                ref_ri_hist: &HashMap<u64,(u64,HashMap<u64,(u64,u64)>)>) -> Vec<PPUC>{

        let ri_hist: Vec<(u64,u64)> = ref_ri_hist.iter().map(|(k,v)|(*k,v.0)).collect();
        let total_count = ri_hist.iter().fold(0,|acc,(_k,v)| acc+v);
        let mut hits = 0;
        let mut head_cost = 0;

        let mut lease_hit_table = HashMap::new();
        let mut lease_cost_table = HashMap::new();

        lease_hit_table.insert(0,0);
        lease_cost_table.insert(0,0);

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
                ppuc:((*lease_hit_table.get(k).unwrap() -*lease_hit_table.get(&base_lease).unwrap()) as f32/
                      (*lease_cost_table.get(k).unwrap()-*lease_cost_table.get(&base_lease).unwrap())as f32),
                lease: *k,
                old_lease: base_lease,
                ref_id: ref_id, 
                new_hits: *lease_hit_table.get(k).unwrap() - *lease_hit_table.get(&base_lease).unwrap(),
            }
        ).collect()
    }

    pub fn shel_cshel(cshel: bool,
                  ri_hists : &RIHists, 
                  cache_size : u64, 
                  sample_rate : u64, 
                  samples_per_phase : HashMap<u64,u64>,
                  verbose: bool,
                  debug: bool) -> Option<(HashMap<u64,u64>,HashMap<u64,(f32,u64)>,u64)>{

        let mut new_lease: PPUC;
        let mut cost_per_phase = HashMap::new();
        let mut budget_per_phase = HashMap::new();
        let mut leases = HashMap::new(); //{ri, lease}
        let mut dual_leases : HashMap<u64,(f32,u64)>= HashMap::new(); //{ri, (alpha, long_lease)}
        let mut num_hits : u64 = 0;
        let mut trace_length : u64 = 0;
        
        let phase_ids: Vec<&u64> = samples_per_phase.keys().collect();

        if verbose {
            println!("---------Dump RI Hists------------");
            super::io::debug::print_ri_hists(&ri_hists);
            println!("---------Dump Samples Per Phase---");
            println!("{:?}",&samples_per_phase);
        }

        //initialize ppucs
        let mut ppuc_tree = BinaryHeap::new();

        for (&ref_id, ri_hist) in ri_hists.ri_hists.iter(){
            let ppuc_vec = get_ppuc(ref_id,0,ri_hist);
            for ppuc in ppuc_vec.iter(){
                ppuc_tree.push(*ppuc);
            }
        }

        //initialize cost + budget
        for (&phase,&num) in samples_per_phase.iter(){
            cost_per_phase.insert(phase, 0);  
            budget_per_phase.insert(phase, num * cache_size * sample_rate);
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

        //initialize leases
        for (&ref_id, _ ) in ri_hists.ri_hists.iter(){
            leases.insert(ref_id,0);
        }

        loop {
            new_lease = match ppuc_tree.pop(){
                Some(i) => i,
                None => return Some((leases,dual_leases,trace_length - num_hits)),
            };

            //continue to pop until we have a ppuc with the right base_lease
            if new_lease.old_lease != *leases.get(&new_lease.ref_id).unwrap(){
                continue;
            }

            if dual_leases.contains_key(&new_lease.ref_id){
                continue;
            }

            let old_lease = *leases.get(&new_lease.ref_id).unwrap();
            //check for capacity
            let mut acceptable_lease = true;
            let mut new_phase_ref_cost = HashMap::new(); 
            for (&phase,&current_cost) in cost_per_phase.iter(){
                let new_cost = match cshel{
                    true => cshel_phase_ref_cost(sample_rate,
                                                      phase,
                                                      new_lease.ref_id,
                                                      old_lease,
                                                      new_lease.lease,
                                                      &ri_hists),
                    false => shel_phase_ref_cost(sample_rate,
                                                      phase,
                                                      new_lease.ref_id,
                                                      old_lease,
                                                      new_lease.lease,
                                                      &ri_hists),
                };

                new_phase_ref_cost.insert(phase,new_cost);
                if (new_cost + current_cost) > *budget_per_phase.get(&phase).unwrap() {
                    acceptable_lease = false;
                }
            }
            if verbose & debug {
                println!("Debug: NEW_PHASE_REF_COST {:?}",&new_phase_ref_cost);
            }

            if acceptable_lease {
                num_hits += new_lease.new_hits * sample_rate;
                //update cache use
                for phase in &phase_ids {
                    cost_per_phase.insert(**phase, 
                                          cost_per_phase.get(*phase).unwrap() + new_phase_ref_cost.get(phase).unwrap()); 
                }
                //update leases
                leases.insert(new_lease.ref_id,new_lease.lease);
                //push new ppucs
                let ppuc_vec = get_ppuc(new_lease.ref_id,
                                        new_lease.lease,
                                        ri_hists.ri_hists.get(&new_lease.ref_id).unwrap());

                for ppuc in ppuc_vec.iter(){
                    ppuc_tree.push(*ppuc);
                }

                if verbose {
                    println!("Assigned lease {:x} to reference ({},{:x}). Predicted miss count : {}", 
                             new_lease.lease, (new_lease.ref_id & 0xFF000000) >> 24, 
                             new_lease.ref_id & 0x00FFFFFF,
                             trace_length - num_hits);
                }
            }

            else {
                //unacceptable lease, must assign a dual lease
                let mut alpha = 1.0;
                for (&phase,&current_cost) in cost_per_phase.iter(){
                    let &phase_ref_cost   = new_phase_ref_cost.get(&phase).unwrap(); 
                    if phase_ref_cost > 0 {
                        if *budget_per_phase.get(&phase).unwrap() < current_cost{
                            println!("
                            ERROR: current cost exceeds budget
                            *budget_per_phase.get(&phase).unwrap():  {}
                            currenc_cost:                            {}
                            ",
                            *budget_per_phase.get(&phase).unwrap(),
                            current_cost
                            );
                            panic!();
                        }

                        let remaining_budget = *budget_per_phase.get(&phase).unwrap() - current_cost; 
                        alpha = super::helpers::float_min(alpha, 
                                                          remaining_budget as f32 / 
                                                            phase_ref_cost as f32);
                    }
                }

                if alpha > 0.0{
                    //update hits
                    num_hits += (new_lease.new_hits as f32 * sample_rate as f32 * alpha).round() as u64;
                    //update cache use
                    for phase in &phase_ids{
                        cost_per_phase.insert(**phase, 
                                          cost_per_phase.get(*phase).unwrap() + 
                                            (*new_phase_ref_cost.get(&phase).unwrap() as f32 * alpha).round()
                                              as u64);
                        //fix floating point precision error leading to "overallocation"
                        if cost_per_phase.get(*phase).unwrap()>budget_per_phase.get(*phase).unwrap(){
                            cost_per_phase.insert(**phase,*budget_per_phase.get(*phase).unwrap());
                        }
                    }
                }

                //update dual lease hashmap
                //inserting with alpha=0 is still valuable, since it tells us to 
                //ignore further lease increases of that reference. 
                dual_leases.insert(new_lease.ref_id,(alpha,new_lease.lease));

                if verbose {
                    println!("Assigned dual lease ({:x},{}) to reference ({},{:x}). Predicted miss count: {}", 
                              new_lease.lease, 
                              alpha, 
                              (new_lease.ref_id & 0xFF000000) >> 24,
                              new_lease.ref_id & 0x00FFFFFF,
                              trace_length - num_hits);
                }
            }
            if verbose & debug { 
                for phase in &phase_ids{
                    println!("Debug: phase: {}",phase);
                    println!("Debug:    cost_per_phase:   {:?}",
                             cost_per_phase.get(&phase).unwrap());
                    println!("Debug:    budget_per_phase: {:?}",
                             budget_per_phase.get(&phase).unwrap());
                }
                println!("=================OUTER LOOP ITER===============");
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;


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
        process_sample_head_cost(&mut ri_hists,1,12,20,(10,1)); //reference 1, phase 0
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

    #[test]
    fn tail_cost_cross_phase(){
        let mut ri_hists = HashMap::new();
        process_sample_head_cost(&mut ri_hists,1,100,170,(100,1));
        process_sample_head_cost(&mut ri_hists,1,50,120,(100,1));

        process_sample_tail_cost(&mut ri_hists,1,100,170,(100,1));
        process_sample_tail_cost(&mut ri_hists,1,50,120,(100,1));

        let hist_struct = RIHists::new(ri_hists);

        print_ri_hists(&hist_struct);
        assert_eq!(hist_struct.get_ref_ri_phase_cost(1,50,0).1,30);
        assert_eq!(hist_struct.get_ref_ri_phase_cost(1,50,1).1,20);
    }
}
