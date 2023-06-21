use serde::{Deserialize, Serialize};
use std::collections::BinaryHeap;
use std::collections::HashMap;
use std::convert::TryInto;
use std::fs::File;
use std::io::Write;

#[derive(Deserialize, Debug)]
struct Sample {
    phase_id_ref: String,
    ri: String,
    tag: String,
    time: u64,
}

#[derive(Serialize)]
struct PhaseTime {
    phase_id: u16,
    start_time: u64,
}

pub fn get_prl_hists(
    input_file: &str,
    num_bins: u64,
    set_mask: u32,
) -> (super::lease_gen::BinnedRIs, super::lease_gen::BinFreqs, u64) {
    let mut curr_bin: u64 = 0;
    let mut curr_bin_dict = HashMap::<u64, u64>::new();
    let mut bin_freqs = HashMap::<u64, HashMap<u64, u64>>::new();
    let mut bin_ri_distributions = HashMap::<u64, HashMap<u64, HashMap<u64, u64>>>::new();
    let mut curr_ri_distribution_dict = HashMap::<u64, HashMap<u64, u64>>::new();
    let mut last_address: u64 = 0;
    let mut all_keys: Vec<u64> = Vec::new();

    bin_freqs.insert(0, curr_bin_dict.clone());
    bin_ri_distributions.insert(0, curr_ri_distribution_dict.clone());

    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .from_path(input_file)
        .unwrap();
    for result in rdr.deserialize() {
        let sample: Sample = result.unwrap();
        last_address = sample.time;
    }
    let bin_width = ((last_address as f64 / (num_bins as f64)) as f64).ceil() as u64;
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .from_path(input_file)
        .unwrap();
    for result in rdr.deserialize() {
        let sample: Sample = result.unwrap();

        //if outside of current bin, moved to the next
        if sample.time > curr_bin + bin_width {
            //store the RI and frequency data for the old bin
            bin_freqs.insert(curr_bin, curr_bin_dict.clone());
            bin_ri_distributions.insert(curr_bin, curr_ri_distribution_dict.clone());
            //initalize storage for the new bin
            curr_bin_dict.clear();
            curr_ri_distribution_dict.clear();
            curr_bin += bin_width;
        }
        let tag = u32::from_str_radix(&sample.tag, 16).unwrap();
        let set = (tag & set_mask) as u64;

        let addr = u64::from_str_radix(&sample.phase_id_ref, 16).unwrap() | (set << 32);
        let ri = u64::from_str_radix(&sample.ri, 16).unwrap();

        if curr_bin_dict.contains_key(&addr) {
            curr_bin_dict.insert(addr, curr_bin_dict.get(&addr).unwrap() + 1);
        } else {
            curr_bin_dict.insert(addr, 1);
        }
        //if the address isn't a key in the outer level hashmap,
        //add a hashmap corresponding to that key; then add the
        //key value pair (ri,1) to that hashmap
        //
        //if the address is a key in the outer level hashmap,
        //check to see if the current ri is a key in the corresponding inner hashmap.
        //if yes, then increment the count for that RI, if not, insert a new key value pair (ri,1).
        *curr_ri_distribution_dict
            .entry(addr)
            .or_insert(HashMap::new())
            .entry(ri)
            .or_insert(0) += 1;

        //create an array of all found references
        if !all_keys.iter().any(|&i| i == addr) {
            all_keys.push(addr);
        }
    }

    //store the frequency and RI data for the last bin
    bin_freqs.insert(curr_bin, curr_bin_dict.clone());
    bin_ri_distributions.insert(curr_bin, curr_ri_distribution_dict.clone());
    //if a reference is not in a bin, add it with a frequency count of 0
    let temp = bin_freqs.clone();
    for (bin, _addrs) in &temp {
        let bin_freqs_temp = bin_freqs.entry(*bin).or_insert(HashMap::new());
        for key in &all_keys {
            bin_freqs_temp.entry(*key).or_insert(0);
        }
    }
    (
        super::lease_gen::BinnedRIs::new(bin_ri_distributions),
        super::lease_gen::BinFreqs::new(bin_freqs),
        bin_width,
    )
}

pub fn build_phase_transitions(input_file: &str) -> (Vec<(u64, u64)>, usize, u64) {
    println!("Reading input from: {}", input_file);
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .from_path(input_file)
        .unwrap();
    let mut u_tags = HashMap::<u64, bool>::new();
    let mut sample_hash = HashMap::new();
    let mut last_sample_time: u64 = 0;
    let mut sample_num: u64 = 0;
    for result in rdr.deserialize() {
        let sample: Sample = result.unwrap();
        let ri = u64::from_str_radix(&sample.ri, 16).unwrap();
        //don't use end of benchmark infinite RIs

        //store unique tags
        u_tags.insert(u64::from_str_radix(&sample.tag, 16).unwrap(), false);
        let phase_id_ref = u64::from_str_radix(&sample.phase_id_ref, 16).unwrap();
        let phase_id = (phase_id_ref & 0xFF000000) >> 24;
        let reuse_time = sample.time;
        let use_time;
        if (ri as i32) < 0 {
            use_time = reuse_time - (ri as i32).abs() as u64;
        } else {
            use_time = reuse_time - ri as u64;
        }
        sample_hash.insert(use_time, phase_id);
        //get empircal sampling rate
        last_sample_time = sample.time;
        sample_num += 1;
    }
    //empircally calculate sampling rate
    let sampling_rate = (last_sample_time as f64 / sample_num as f64).round() as u64;
    println!("empirical sampling_rate:{}", sampling_rate);
    //every data block is associated with at least one miss in the absense of hardware prefetching.
    let first_misses = u_tags.len();

    let mut sorted_samples: Vec<_> = sample_hash.iter().collect();
    sorted_samples.sort_by_key(|a| a.0);

    //get phase transitions
    let mut phase_transitions = HashMap::new(); //(time,phase start)
    let mut phase = 0;
    phase_transitions.insert(0, phase);

    for s in sorted_samples.iter() {
        if *s.1 != phase {
            phase_transitions.insert(*s.0, *s.1);
            phase = *s.1;
        }
    }

    let mut sorted_transitions: Vec<_> = phase_transitions.iter().collect();
    sorted_transitions.sort_by_key(|a| a.0);

    return (
        sorted_transitions
            .iter()
            .map(|&x| (*(x.0), *(x.1)))
            .collect(),
        first_misses,
        sampling_rate,
    );
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
pub fn build_ri_hists(
    input_file: &str,
    cshel: bool,
    set_mask: u32,
) -> (super::lease_gen::RIHists, HashMap<u64, u64>, usize, u64) {
    let (phase_transitions, first_misses, sampling_rate) = build_phase_transitions(input_file);
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
            let ri = u32::from_str_radix(&sample.ri, 16).unwrap();
            let reuse_time = sample.time;
            let use_time;

            let mut ri_signed = ri as i32;

            if ri_signed < 0 {
                use_time = reuse_time - (!ri_signed + 1) as u64;
                ri_signed = 16777215; //canonical value for negatives
            } else {
                use_time = reuse_time - ri_signed as u64;
            }

            let phase_id_ref = u64::from_str_radix(&sample.phase_id_ref, 16).unwrap();

            let tag = u32::from_str_radix(&sample.tag, 16).unwrap();
            let set = (tag & set_mask) as u64;
            let set_phase_id_ref = phase_id_ref | set << 32;

            let next_phase_tuple = match super::helpers::binary_search(&phase_transitions, use_time)
            {
                Some(v) => v,
                None => (reuse_time + 1, 0),
            };
            super::lease_gen::process_sample_head_cost(
                &mut ri_hists,
                set_phase_id_ref,
                ri_signed as u64,
                use_time,
                next_phase_tuple,
            );

            let phase_id = (phase_id_ref & 0xFF000000) >> 24;
            *samples_per_phase.entry(phase_id).or_insert_with(|| 0) += 1;
        }

        let mut rdr = csv::ReaderBuilder::new()
            .has_headers(false)
            .from_path(input_file)
            .unwrap();
        println!("before second pass!");
        //second pass for tail costs
        for result in rdr.deserialize() {
            let sample: Sample = result.unwrap();
            let ri = u32::from_str_radix(&sample.ri, 16).unwrap();

            let reuse_time = sample.time;
            let use_time;

            let mut ri_signed = ri as i32;

            if ri_signed < 0 {
                use_time = reuse_time - (!ri_signed + 1) as u64;
                ri_signed = 16777215; //canonical value for negatives
            } else {
                use_time = reuse_time - ri_signed as u64;
            }

            let phase_id_ref = u64::from_str_radix(&sample.phase_id_ref, 16).unwrap();
            let tag = u32::from_str_radix(&sample.tag, 16).unwrap();
            let set = (tag & set_mask) as u64;
            let set_phase_id_ref = phase_id_ref | set << 32;

            let next_phase_tuple = match super::helpers::binary_search(&phase_transitions, use_time)
            {
                Some(v) => v,
                None => (reuse_time + 1, 0),
            };

            super::lease_gen::process_sample_tail_cost(
                &mut ri_hists,
                set_phase_id_ref,
                ri_signed as u64,
                use_time,
                next_phase_tuple,
            );
        }
    }
    //if not doing C-SHEL generates RI distribution with head
    //and tail costs set to 0 (significantly faster)
    else {
        for result in rdr.deserialize() {
            let sample: Sample = result.unwrap();
            let mut ri = u64::from_str_radix(&sample.ri, 16).unwrap();
            let _reuse_time = sample.time;
            //if sample is negative, there is no reuse
            let ri_signed = ri as i32;
            if ri_signed < 0 {
                ri = 16777215;
            }

            let phase_id_ref = u64::from_str_radix(&sample.phase_id_ref, 16).unwrap();
            let tag = u32::from_str_radix(&sample.tag, 16).unwrap();
            let set = (tag & set_mask) as u64;
            let set_phase_id_ref = phase_id_ref | set << 32;
            let phase_id = (phase_id_ref & 0xFF000000) >> 24;
            *samples_per_phase.entry(phase_id).or_insert(0) += 1;
            //if the reference isn't in the distrubtion add it
            //if an ri for that reference isn't in the distrubtion
            //add it as key with a value of (ri_count,{phaseID,(0,0)))
            //if an ri for that reference is in the distribution
            //increment the ri count by 1
            ri_hists
                .entry(set_phase_id_ref)
                .or_insert(HashMap::new())
                .entry(ri)
                .and_modify(|e| e.0 += 1)
                .or_insert((1, HashMap::new()))
                .1
                .entry(phase_id)
                .or_insert((0, 0));
        }
    }
    (
        super::lease_gen::RIHists::new(ri_hists),
        samples_per_phase,
        first_misses,
        sampling_rate,
    )
}

pub fn dump_leases(
    leases: HashMap<u64, u64>,
    dual_leases: HashMap<u64, (f64, u64)>,
    lease_hits: HashMap<u64, HashMap<u64, u64>>,
    trace_length: u64,
    output_file: &str,
    sampling_rate: u64,
    first_misses: usize,
) -> Vec<(u64, u64, u64, u64, f64)> {
    let mut num_hits = 0;
    //create lease output vector
    let mut lease_vector: Vec<(u64, u64, u64, u64, f64)> = Vec::new();
    for (&phase_address, &lease) in leases.iter() {
        let lease = if lease > 0 { lease } else { 1 };
        let phase = (phase_address & 0xFF000000) >> 24;
        let address = phase_address & 0x00FFFFFF;
        if dual_leases.contains_key(&phase_address) {
            lease_vector.push((
                phase,
                address,
                lease,
                dual_leases.get(&phase_address).unwrap().1,
                1.0 - dual_leases.get(&phase_address).unwrap().0,
            ));
        } else {
            lease_vector.push((phase, address, lease, 0, 1.0));
        }
    }
    lease_vector.sort_by_key(|a| (a.0, a.1)); //sort by phase and then by reference
                                              //get number of predicted misses
    for (phase, address, lease_short, lease_long, percentage) in lease_vector.iter() {
        //reassemble phase address
        let phase_address = address | phase << 24;
        //we are assuming that our sampling captures all RIS
        //by assuming the distribution is normal
        //thus if an RI for a reference didn't occur during runtime
        //(i.e., the base lease of 1 that all references get)
        //we can assume the number of hits it gets is zero.
        if lease_hits.get(&phase_address).unwrap().get(lease_short) != None {
            num_hits += (*lease_hits
                .get(&phase_address)
                .unwrap()
                .get(lease_short)
                .unwrap() as f64
                * (percentage))
                .round() as u64;
        }
        if lease_hits.get(&phase_address).unwrap().get(lease_long) != None {
            num_hits += (*lease_hits
                .get(&phase_address)
                .unwrap()
                .get(lease_long)
                .unwrap() as f64
                * (1.0 - percentage))
                .round() as u64;
        }
    }
    println!("Writing output to: {}", output_file);
    let mut file = File::create(output_file).expect("create failed");

    file.write_all(
        &format!(
            "Dump predicted miss count (no contention misses): {}\n",
            trace_length - num_hits * sampling_rate + first_misses as u64
        )[..]
            .as_bytes(),
    )
    .expect("write failed");

    file.write_all("Dump formated leases\n".as_bytes())
        .expect("write failed");

    for (phase, address, lease_short, lease_long, percentage) in lease_vector.iter() {
        file.write_all(
            &format!(
                "{:x}, {:x}, {:x}, {:x}, {}\n",
                phase, address, lease_short, lease_long, percentage
            )[..]
                .as_bytes(),
        )
        .expect("write failed");
    }
    return lease_vector;
}
// function for generating c-files
pub fn gen_lease_c_file(
    mut lease_vector: Vec<(u64, u64, u64, u64, f64)>,
    llt_size: u64,
    max_num_scopes: u64,
    mem_size: u64,
    output_file: String,
    discretize_width: u64,
) {
    let mut phase_lease_arr: HashMap<u64, HashMap<u64, (u64, u64, f64, bool)>> = HashMap::new();
    let mut phases: Vec<u64> = Vec::new();
    for lease in lease_vector.iter() {
        if !phases.contains(&lease.0) {
            phases.push(lease.0);
        }
    }
    //due to the way the lease cache operates, phases skipped
    //if there are phases with no leases, assign a dummy lease to the skipped phase
    //Since we have no way of knowing without adding
    //another dependency how many phases a program has
    //create dummy phases for all phases that can fit in memory that aren't represented
    for phase in 0..max_num_scopes {
        if !phases.contains(&phase) {
            lease_vector.push((phase, 0, 0, 0, 1.0));
        }
    }

    //convert lease vector to hashmap of leases per phase
    for (phase, address, lease_short, lease_long, percentage) in lease_vector.iter() {
        phase_lease_arr
            .entry(*phase)
            .or_insert(HashMap::new())
            .entry(*address)
            .or_insert((
                *lease_short,
                *lease_long,
                *percentage,
                if lease_long > &0 { true } else { false },
            ));
    }
    let default_lease = 1;

    //make sure each phase can fit in the specified LLT
    for (phase, phase_leases) in phase_lease_arr.iter() {
        if phase_leases.len() > llt_size as usize {
            println!(
                "Leases for Phase {} don't fit in lease lookup table!",
                phase
            );
            panic!();
        }
    }

    //make sure that all phases can fit in the memory allocated
    if *phases.iter().max().unwrap() > max_num_scopes {
        println!(
            "Error: phases cannot fit in specified {} byte memory",
            mem_size
        );
        panic!();
    }

    //write header
    let mut file = std::fs::File::create(output_file).expect("create failed");
    file.write_all("#include \"stdint.h\"\n\n".as_bytes())
        .expect("write failed");
    file.write_all(
            format!(
                "static uint32_t lease[{}] __attribute__((section (\".lease\"))) __attribute__ ((__used__)) = {{\n",
                mem_size/4)
            .as_bytes())
            .expect("write failed");
    file.write_all("// lease header\n".as_bytes())
        .expect("write failed");
    let mut phase_index: u64 = 0; //len returns usize which can't directly substituted as u64
    for i in 0..phase_lease_arr.len() {
        let phase_leases = phase_lease_arr.get(&phase_index).unwrap();
        phase_index = phase_index + 1; //increment to next phase
        file.write_all(format!("// phase {}\n", i).as_bytes())
            .expect("write failed");

        let mut dual_lease_ref = (0, 0, 1.0);
        let mut lease_phase: Vec<(u64, u64)> = Vec::new();
        let dual_lease_found = false;
        for (lease_ref, lease_data) in phase_leases.iter() {
            //convert hashmap of leases for phase to vector
            lease_phase.push((*lease_ref, lease_data.0));
            //get dual lease if it exists;
            if lease_data.3 == true && dual_lease_found == false {
                dual_lease_ref = (*lease_ref, lease_data.1, lease_data.2);
            }
        }
        lease_phase.sort_by_key(|a| a.0);
        //output config
        for j in 0..16 {
            if j == 0 {
                file.write_all(
                    format!("\t0x{:08x},\t// default lease\n", default_lease).as_bytes(),
                )
                .expect("write failed");
            } else if j == 1 {
                file.write_all(
                    format!("\t0x{:08x},\t// long lease value\n", dual_lease_ref.1).as_bytes(),
                )
                .expect("write failed");
            } else if j == 2 {
                file.write_all(
                    format!(
                        "\t0x{:08x},\t// short lease probability\n",
                        discretize(dual_lease_ref.2, discretize_width)
                    )
                    .as_bytes(),
                )
                .expect("write failed");
            } else if j == 3 {
                file.write_all(
                    format!(
                        "\t0x{:08x},\t// num of references in phase\n",
                        phase_leases.len()
                    )
                    .as_bytes(),
                )
                .expect("write failed");
            } else if j == 4 {
                file.write_all(
                    format!(
                        "\t0x{:08x},\t// dual lease ref (word address)\n",
                        dual_lease_ref.0 >> 2
                    )
                    .as_bytes(),
                )
                .expect("write failed");
            } else {
                file.write_all(format!("\t0x{:08x},\t // unused\n", 0).as_bytes())
                    .expect("write failed");
            }
        }
        let field_list = ["reference address", "lease0 value"];

        // loop through lease fields
        for k in 0..2 {
            file.write_all(format!("\t//{}\n\t", field_list[k]).as_bytes())
                .expect("write failed");

            for j in 0..llt_size {
                if j < phase_leases.len().try_into().unwrap() {
                    if k == 0 {
                        file.write_all(format!("0x{:08x}", lease_phase[j as usize].0).as_bytes())
                            .expect("write failed");
                    } else {
                        file.write_all(format!("0x{:08x}", lease_phase[j as usize].1).as_bytes())
                            .expect("write failed");
                    }
                } else {
                    file.write_all(format!("0x{:08x}", 0).as_bytes())
                        .expect("write failed");
                }
                //print delimiter
                if j + 1 == llt_size && k == 1 && i + 1 == phase_lease_arr.len() {
                    file.write_all(format!("\n").as_bytes())
                        .expect("write failed");
                } else if j + 1 == llt_size {
                    file.write_all(format!(",\n").as_bytes())
                        .expect("write failed");
                } else if ((j + 1) % 10) == 0 {
                    file.write_all(format!(",\n\t").as_bytes())
                        .expect("write failed");
                } else {
                    file.write_all(format!(", ").as_bytes())
                        .expect("write failed");
                }
            }
        }
    }
    file.write_all(format!("}};").as_bytes())
        .expect("write failed");
}

pub fn discretize(percentage: f64, discretization: u64) -> u64 {
    let percentage_binary =
        (percentage * ((2 << (discretization - 1)) as f64) - 1.0).round() as u64;
    return percentage_binary;
}

pub mod debug {
    pub fn print_ri_hists(rihists: &super::super::lease_gen::RIHists) {
        for (ref_id, ref_ri_hist) in &rihists.ri_hists {
            println!(
                "({},0x{:x}):",
                (ref_id & 0xFF000000) >> 24,
                ref_id & 0x00FFFFFF
            );
            for (ri, tuple) in ref_ri_hist {
                println!(" | ri 0x{:x}: count {}", ri, tuple.0);
                for (phase_id, cost) in &tuple.1 {
                    println!(
                        " | | phase {} head_cost {} tail_cost {}",
                        phase_id, cost.0, cost.1
                    );
                }
            }
        }
    }
    //<u64,HashMap<u64,HashMap<u64,u64>>>
    pub fn print_binned_hists(binned_ris: &super::super::lease_gen::BinnedRIs) {
        for (bin, ref_ri_hist) in &binned_ris.bin_ri_distribution {
            println!("Bin:{}", bin);
            for (ref_id, ri_hist) in ref_ri_hist {
                println!(" | ref 0x{:x}:", ref_id);
                for (ri, count) in ri_hist {
                    println!(" | | ri 0x{:x}: count {}", ri, count);
                }
            }
        }
    }

    pub fn destructive_print_ppuc_tree(
        ppuc_tree: &mut super::BinaryHeap<super::super::lease_gen::PPUC>,
    ) {
        while ppuc_tree.peek() != None {
            println!("ppuc: {:?}", ppuc_tree.pop().unwrap());
        }
    }
}
