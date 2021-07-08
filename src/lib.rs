use std::collections::HashMap;
use std::collections::BinaryHeap;
use core::cmp::Ordering;
//use csv::{ReaderBuilder};
use serde::{Serialize, Deserialize};
//use std::fmt;
use std::u64;

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

pub struct RIHists{
    ri_hists: HashMap<u64,HashMap<u64,(u64,HashMap<u64,(u64,u64)>)>>,
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
}
impl PartialOrd for PPUC {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering>{
        other.ppuc.partial_cmp(&self.ppuc)
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

fn build_phase_transitions(input_file:&str) -> Vec<(u64,u64)>{
    println!("Reading input from: {}", input_file);

    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .from_path(input_file)
        .unwrap();

    let mut sample_hash = HashMap::new();

    for result in rdr.deserialize() {
        let sample: Sample = result.unwrap();
        let ri = u64::from_str_radix(&sample.ri,16).unwrap();
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

    //sorted_transitions.into_iter().cloned().collect()
}

fn binary_search(vector: &Vec<(u64,u64)>,value:u64) -> Option<(u64,u64)>{
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

fn process_sample_head_cost(ri_hists: &mut HashMap<u64,HashMap<u64,(u64,HashMap<u64,(u64,u64)>)>>,
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
    let next_phase_head_cost= std::cmp::max(time as i64 - next_phase_tuple.0 as i64,0) as u64;
    ri_tuple.1.entry(phase_id).or_insert_with(|| (0,0)).0 += this_phase_head_cost;

    if next_phase_head_cost>0 {
        ri_tuple.1.entry(next_phase_tuple.1).or_insert_with(|| (0,0)).0 += next_phase_head_cost;
    }
}

fn process_sample_tail_cost(ri_hists: &mut HashMap<u64,HashMap<u64,(u64,HashMap<u64,(u64,u64)>)>>,
                  phase_id_ref: u64,
                  ri: u64,
                  time: u64,
                  next_phase_tuple: (u64,u64)){

    let phase_id = (phase_id_ref & 0xFF000000)>>24;
    let start_time = time - ri;

    let ref_hist = ri_hists.entry(phase_id_ref).or_insert_with(|| HashMap::new());

    //this heinous code exists so we can iterate through a hashmap while modifying it
    let ris : Vec<&u64> = ref_hist.keys().collect();
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
        let count_phase_cost_tuple = ref_hist.entry(ri_other).or_insert_with(|| (0,HashMap::new()));

        let this_phase_tail_cost = std::cmp::min(next_phase_tuple.0 - start_time, ri_other);
        let next_phase_tail_cost = std::cmp::max(0,start_time as i64 + ri_other as i64 - next_phase_tuple.0 as i64) as u64;

        count_phase_cost_tuple.1.entry(phase_id).or_insert_with(|| (0,0)).1 += this_phase_tail_cost;

        if next_phase_tail_cost > 0 {
            count_phase_cost_tuple.1.entry(next_phase_tuple.1).or_insert_with(|| (0,0)).1 += next_phase_tail_cost;
        }
    }
}

//Build ri hists in the following form
//{ref_id,
//  {ri,
//    (count,{phase_id,cost})}}

//TODO: update this function to accumulate tail costs and return the form:
//{ref_id,
//  {ri,
//    (count,{phase_id,(head_cost,tail_cost)})}}
//
// Semantically, head cost refers to the accumulation of cost from reuses with length ri which span 
// phase boundaries.
//
// Tail cost refers to the accumulation of cost from reuses greater than ri, whose lease spans a
// phase boundary at l=ri
pub fn build_ri_hists(input_file:&str)-> (RIHists,HashMap<u64,u64>){
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
        let time = sample.time;
        let phase_id_ref = u64::from_str_radix(&sample.phase_id_ref,16).unwrap();
        let next_phase_tuple = match binary_search(&phase_transitions,time-ri){
            Some(v) => v,
            None => (time+1,0),
        };

        process_sample_head_cost(&mut ri_hists,phase_id_ref,ri,time,next_phase_tuple);

        let phase_id = (phase_id_ref & 0xFF000000)>>24;
        *samples_per_phase.entry(phase_id).or_insert_with(|| 0)+=1;
    }

    //second pass for tail costs
    for result in rdr.deserialize() {
        let sample: Sample = result.unwrap();
        let ri = u64::from_str_radix(&sample.ri,16).unwrap();
        let time = sample.time;
        let phase_id_ref = u64::from_str_radix(&sample.phase_id_ref,16).unwrap();
        let next_phase_tuple = match binary_search(&phase_transitions,time-ri){
            Some(v) => v,
            None => (time+1,0),
        };
        process_sample_tail_cost(&mut ri_hists,phase_id_ref,ri,time,next_phase_tuple);
    }

    (RIHists::new(ri_hists),samples_per_phase)
}

pub fn print_ri_hists(rihists: &RIHists){
    for (ref_id, ref_ri_hist) in &rihists.ri_hists{
        println!("{}:",ref_id);
        for (ri,tuple) in ref_ri_hist{
            println!(" | ri {}: count {}",ri, tuple.0);
            for (phase_id,cost) in &tuple.1{
                println!(" | | phase {} head_cost {} tail_cost {}",phase_id,cost.0,cost.1);
            }
        }
    }
}


fn get_phase_ref_cost(phase: u64,ref_id: u64,old_lease: u64,new_lease: u64,ri_hists: &RIHists) -> u64 {

    //putting this here for the compiler to not yell at me
    println!("{} {} {} {}", phase,ref_id,old_lease,new_lease);
    print_ri_hists(&ri_hists);
    0
/*
  //TODO rewrite this function to account for head + tail costs
    let mut lease_cost_table = HashMap::new();
    lease_cost_table.insert(0,0);
    let mut total_count = 0;

    let ri_hist = ri_hists.ri_hists.get(&ref_id).unwrap();
    let mut ri_hist_clone = ri_hist.clone();

    let head_cost = 0;
    ri_hist_clone.sort_by(|a,b| a.0.cmp(&b.0));
    for (&ri,&count) in ri_hist_clone.iter(){
        head_cost += count *  ri;
        let tail_cost = (total_count - hits) * ri;
        lease_cost_table.insert(ri, head_cost + tail_cost);
    }

    lease_cost_table.get(&new_lease).unwrap()-lease_cost_table.get(&old_lease).unwrap()
*/
}

pub fn c_shel(ri_hists : &RIHists, 
              cache_size : u64, 
              sample_rate : u64, 
              samples_per_phase : HashMap<u64,u64>) -> Option<(HashMap<u64,u64>,HashMap<u64,(f32,u64)>)>{

    let mut new_lease: PPUC;
    let mut cost_per_phase = HashMap::new();
    let mut budget_per_phase = HashMap::new();
    let mut leases = HashMap::new(); //{ri, lease}
    let mut dual_leases : HashMap<u64,(f32,u64)>= HashMap::new(); //{ri, (alpha, long_lease)}
    
    let phase_ids: Vec<&u64> = samples_per_phase.keys().collect();

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
    }

    //initialize leases
    for (&ref_id, _ ) in ri_hists.ri_hists.iter(){
        leases.insert(ref_id,0);
    }

    loop {
        new_lease = match ppuc_tree.pop(){
            Some(i) => i,
            None => return Some((leases,dual_leases)),
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
        for (&phase,&current_cost) in cost_per_phase.iter(){
            let new_cost = get_phase_ref_cost(phase,new_lease.ref_id,old_lease,new_lease.lease,&ri_hists);
            if (new_cost + current_cost) > *budget_per_phase.get(&phase).unwrap() {
                acceptable_lease = false;
                break;
            }
        }

        if acceptable_lease {
            //update cache use
            for phase in &phase_ids{
                cost_per_phase.insert(**phase, cost_per_phase.get(*phase).unwrap() + get_phase_ref_cost(**phase,new_lease.ref_id,old_lease,new_lease.lease,&ri_hists));
            }
            //update leases
            leases.insert(new_lease.ref_id,new_lease.lease);

            //push new ppucs
            let ppuc_vec = get_ppuc(new_lease.ref_id,new_lease.lease,ri_hists.ri_hists.get(&new_lease.ref_id).unwrap());
            for ppuc in ppuc_vec.iter(){
                ppuc_tree.push(*ppuc);
            }
        }

        else {
            //unacceptable lease, must assign a dual lease
            let mut alpha = 1.0;
            for (&phase,&current_cost) in cost_per_phase.iter(){
                alpha = float_min(alpha,
                                  get_phase_ref_cost(phase,new_lease.ref_id,old_lease,new_lease.lease,&ri_hists) as f32/
                                 (*budget_per_phase.get(&phase).unwrap() - current_cost) as f32);
            }

            //update cache use
            for phase in &phase_ids{
                cost_per_phase.insert(**phase, 
                                      cost_per_phase.get(*phase).unwrap() + 
                                         (get_phase_ref_cost(**phase,new_lease.ref_id,old_lease,new_lease.lease,&ri_hists) 
                                            as f32 * alpha).round() as u64);
            }

            //update dual lease hashmap
            dual_leases.insert(new_lease.ref_id,(alpha,new_lease.lease));
        }
    }
    //None
}

fn float_min(a: f32, b:f32) -> f32{
    if a.lt(&b){ 
        return a;
    }
    b
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
        }
    ).collect()
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
    fn test_cross_phase_cost(){

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
}
