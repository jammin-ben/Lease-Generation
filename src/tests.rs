use std::collections::HashMap;

use super::helpers::*;
use super::io::debug::*;
use super::lease_gen::*;

#[test]
fn it_works() {
    assert_eq!(2 + 2, 4);
}

#[test]
fn test_bin_search() {
    let a = vec![(0, 0), (5, 1), (120, 2), (288, 0), (1025, 1)];
    assert_eq!(binary_search(&a, 2), Some((5, 1)));
    assert_eq!(binary_search(&a, 120), Some((288, 0)));
    assert_eq!(binary_search(&a, 1024), Some((1025, 1)));
    assert_eq!(binary_search(&a, 2048), None);
    assert_eq!(binary_search(&a, 1025), None);
}

#[test]
fn test_process_sample_head_cost() {
    let mut ri_hists = HashMap::new();
    process_sample_head_cost(&mut ri_hists, 1, 12, 20, (100000, 1));
    let hist_struct = RIHists::new(ri_hists);

    assert_eq!(hist_struct.get_ref_ri_count(1, 12), 1);
    assert_eq!(hist_struct.get_ref_ri_phase_cost(1, 12, 0), (12, 0));
}

#[test]
fn test_cross_phase_head_cost() {
    let mut ri_hists = HashMap::new();
    process_sample_head_cost(&mut ri_hists, 1, 12, 8, (10, 1)); //reference 1, phase 0
    let hist_struct = RIHists::new(ri_hists);

    assert_eq!(hist_struct.get_ref_ri_count(1, 12), 1);
    assert_eq!(hist_struct.get_ref_ri_phase_cost(1, 12, 0), (2, 0));
    assert_eq!(hist_struct.get_ref_ri_phase_cost(1, 12, 1), (10, 0));
}

#[test]
fn test_get_ppuc() {
    let mut ri_hist = HashMap::new();

    ri_hist.insert(3, (1, HashMap::new()));
    ri_hist.insert(5, (7, HashMap::new()));
    ri_hist.insert(17, (4, HashMap::new()));
    ri_hist.insert(19, (3, HashMap::new()));
    let ppucs = get_ppuc(1, 0, &ri_hist);
    println!("{:?}", ppucs); //According to Asplos19, the numbers are [.02,.11,.08,.09];
}

#[test]
fn test_process_sample_tail_cost() {
    let mut ri_hists = HashMap::new();
    let ri_short = 10;
    let ri_long = 100;
    let ri_very_long = 1000;
    process_sample_head_cost(&mut ri_hists, 1, ri_short, 10, (10000, 1));
    process_sample_head_cost(&mut ri_hists, 1, ri_short, 20, (10000, 1));
    process_sample_head_cost(&mut ri_hists, 1, ri_long, 200, (10000, 1));
    process_sample_head_cost(&mut ri_hists, 1, ri_very_long, 2200, (10000, 1));

    process_sample_tail_cost(&mut ri_hists, 1, ri_short, 10, (10000, 1));
    process_sample_tail_cost(&mut ri_hists, 1, ri_short, 20, (10000, 1));
    process_sample_tail_cost(&mut ri_hists, 1, ri_long, 200, (10000, 1));
    process_sample_tail_cost(&mut ri_hists, 1, ri_very_long, 2200, (10000, 1));

    let hist_struct = RIHists::new(ri_hists);

    // print_ri_hists(&hist_struct);
    assert_eq!(hist_struct.get_ref_ri_phase_cost(1, ri_short, 0).1, 20);
    assert_eq!(hist_struct.get_ref_ri_phase_cost(1, ri_long, 0).1, 100);
    assert_eq!(hist_struct.get_ref_ri_phase_cost(1, ri_very_long, 0).1, 0);
}

//pub fn process_sample_head_cost(ri_hists: &mut HashMap<u64,HashMap<u64,(u64,HashMap<u64,(u64,u64)>)>>,
//p                 phase_id_ref: u64,
//p              ri: u64,
//p           use_time: u64,
//p        next_phase_tuple: (u64,u64)){
#[test]
fn tail_cost_cross_phase() {
    let mut ri_hists = HashMap::new();
    process_sample_head_cost(&mut ri_hists, 1, 100, 70, (100, 1));
    process_sample_head_cost(&mut ri_hists, 1, 50, 70, (100, 1));

    process_sample_tail_cost(&mut ri_hists, 1, 100, 70, (100, 1));
    process_sample_tail_cost(&mut ri_hists, 1, 50, 70, (100, 1));

    let hist_struct = RIHists::new(ri_hists);

    print_ri_hists(&hist_struct);
    assert_eq!(hist_struct.get_ref_ri_phase_cost(1, 50, 0).1, 30);
    assert_eq!(hist_struct.get_ref_ri_phase_cost(1, 50, 1).1, 20);
}

#[test]
fn negative_ri() {
    let mut ri_hists = HashMap::new();
    process_sample_head_cost(&mut ri_hists, 1, 50, 80, (100, 1));
    process_sample_head_cost(&mut ri_hists, 1, i32::max as u64, 90, (100, 1));

    process_sample_tail_cost(&mut ri_hists, 1, 50, 80, (100, 1));
    process_sample_tail_cost(&mut ri_hists, 1, i32::max as u64, 90, (100, 1));
    let hist_struct = RIHists::new(ri_hists);

    assert_eq!(hist_struct.get_ref_ri_phase_cost(1, 50, 0).0, 20);
    assert_eq!(hist_struct.get_ref_ri_phase_cost(1, 50, 1).0, 30);
    assert_eq!(hist_struct.get_ref_ri_phase_cost(1, 50, 0).1, 10);
    assert_eq!(hist_struct.get_ref_ri_phase_cost(1, 50, 1).1, 40);
}
