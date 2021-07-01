use std::collections::HashMap;
use csv::{ReaderBuilder};
use serde::{Serialize, Deserialize};
use std::fmt;
use std::i64;

#[derive(Deserialize)]
#[derive(Debug)]
struct Sample{
    phase_id_ref: String,
    ri: String,
    tag: String,
    time: i64,
}

#[derive(Serialize)]
struct PhaseTime{
    phase_id: u16,
    start_time: u64,
}

pub struct RIHists{
    ri_hists: HashMap<i64,HashMap<i64,(u64,HashMap<i64,i64>)>>,
}

impl RIHists {
    pub fn new(ri_hists_input: HashMap<i64,HashMap<i64,(u64,HashMap<i64,i64>)>>) -> Self{
        RIHists{
            ri_hists: ri_hists_input, 
        }
    }

    pub fn get_ref_hist(&self,ref_id: i64) -> &HashMap<i64,(u64,HashMap<i64,i64>)>{
        return self.ri_hists.get(&ref_id).unwrap();
    }

    pub fn get_ref_ri_count(&self,ref_id: i64,ri:i64) -> u64 {
        return self.ri_hists.get(&ref_id).unwrap().get(&ri).unwrap().0;
    }

    pub fn get_ref_ri_cost(&self,ref_id: i64,ri:i64) -> &HashMap<i64,i64> {
        return &self.ri_hists.get(&ref_id).unwrap().get(&ri).unwrap().1;
    }

    pub fn get_ref_ri_phase_cost(&self,ref_id: i64, ri:i64, phase: i64) -> i64 {
        return *self.ri_hists.get(&ref_id).unwrap().get(&ri).unwrap().1.get(&phase).unwrap();
    }

}

fn build_phase_transitions(input_file:&str) -> Vec<(i64,i64)>{
    println!("Reading input from: {}", input_file);

    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .from_path(input_file)
        .unwrap();

    let mut sample_hash = HashMap::new();

    for result in rdr.deserialize() {
        let sample: Sample = result.unwrap();
        let ri = i64::from_str_radix(&sample.ri,16).unwrap();
        let mut time = sample.time;
        let phase_id_ref = i64::from_str_radix(&sample.phase_id_ref,16).unwrap();

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

fn binary_search(vector: &Vec<(i64,i64)>,value:i64) -> Option<(i64,i64)>{
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

fn process_sample(ri_hists: &mut HashMap<i64,HashMap<i64,(u64,HashMap<i64,i64>)>>,
                  phase_id_ref: i64,
                  ri: i64,
                  time: i64,
                  next_phase_tuple: (i64,i64)){

    let phase_id = (phase_id_ref & 0xFF000000)>>24;
    let start_time = time - ri;

    //increment count
    //let ref_hist = ri_hists.get(&phase_id_ref).or_else(|| Some(&Box::new(HashMap::new()))).unwrap();
    //let ri_tuple = ref_hist.get(&ri).or_else(|| Some(&(0,HashMap::new()))).unwrap();
    let ref_hist = ri_hists.entry(phase_id_ref).or_insert_with(|| HashMap::new());
    let ri_tuple = ref_hist.entry(ri).or_insert_with(|| (0,HashMap::new()));
    ri_tuple.0 +=1;

    //increment costs
    let this_phase_cost= std::cmp::min(next_phase_tuple.0 - start_time,ri);
    let next_phase_cost= std::cmp::max(time - next_phase_tuple.0,0);
    *ri_tuple.1.entry(phase_id).or_insert_with(|| 0)+=this_phase_cost;

    if next_phase_cost>0 {
        *ri_tuple.1.entry(next_phase_tuple.1).or_insert_with(|| 0)+=next_phase_cost;
    }

    //println!("{:?}",ri_hists);
}
//Build ri hists in the following form
//{ref_id,
//  {ri,
//    (count,{phase_id,cost})}}
pub fn build_ri_hists(input_file:&str)-> RIHists{
    let phase_transitions = build_phase_transitions(input_file);
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .from_path(input_file)
        .unwrap();


    let mut ri_hists = HashMap::new();

    for result in rdr.deserialize() {
        
        let sample: Sample = result.unwrap();
        let ri = i64::from_str_radix(&sample.ri,16).unwrap();
        let mut time = sample.time;
        let phase_id_ref = i64::from_str_radix(&sample.phase_id_ref,16).unwrap();
        let next_phase_tuple = match binary_search(&phase_transitions,time-ri){
            Some(v) => v,
            None => (time+1,0),
        };
        process_sample(&mut ri_hists,phase_id_ref,ri,time,next_phase_tuple);
    }
    RIHists::new(ri_hists)
}

pub fn print_ri_hists(rihists: &RIHists){
    for (ref_id, ref_ri_hist) in &rihists.ri_hists{
        println!("{}:",ref_id);
        for (ri,tuple) in ref_ri_hist{
            println!(" | ri {}: count {}",ri, tuple.0);
            for (phase_id,cost) in &tuple.1{
                println!(" | | phase {} cost {}",phase_id,cost);
            }
        }
        //println!("{}: {:?}",ref_id,hashmap);
    }
}

fn c_shel(ri_hists:i32)->i32{
    0
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
        assert_eq!(binary_search(&a,-10),Some((0,0)));
        assert_eq!(binary_search(&a,1025),None);
    }

    #[test]
    fn test_process_sample(){
        let mut ri_hists = HashMap::new();
        process_sample(&mut ri_hists,1,12,20,(100000,1));
        let hist_struct = RIHists::new(ri_hists);

        assert_eq!(hist_struct.get_ref_ri_count(1,12),1);
        assert_eq!(hist_struct.get_ref_ri_phase_cost(1,12,0),12);

    }
    #[test]
    fn test_cross_phase_cost(){
        let mut ri_hists = HashMap::new();
        process_sample(&mut ri_hists,1,12,20,(10,1));//reference 1, phase 0
        let hist_struct = RIHists::new(ri_hists);

        assert_eq!(hist_struct.get_ref_ri_count(1,12),1);
        assert_eq!(hist_struct.get_ref_ri_phase_cost(1,12,0),2);
        assert_eq!(hist_struct.get_ref_ri_phase_cost(1,12,1),10);
    }
}
