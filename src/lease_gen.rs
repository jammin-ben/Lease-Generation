use core::cmp::Ordering;
use std::collections::BinaryHeap;
use std::collections::HashMap;
#[derive(Debug, Clone)]
pub struct BinFreqs {
    pub bin_freqs: HashMap<u64, HashMap<u64, u64>>,
}
#[derive(Debug, Clone)]
pub struct BinnedRIs {
    pub bin_ri_distribution: HashMap<u64, HashMap<u64, HashMap<u64, u64>>>,
}

impl BinFreqs {
    pub fn new(bin_freqs_input: HashMap<u64, HashMap<u64, u64>>) -> Self {
        BinFreqs {
            bin_freqs: bin_freqs_input,
        }
    }
}

impl BinnedRIs {
    pub fn new(bin_ri_input: HashMap<u64, HashMap<u64, HashMap<u64, u64>>>) -> Self {
        BinnedRIs {
            bin_ri_distribution: bin_ri_input,
        }
    }
}

pub struct RIHists {
    pub ri_hists: HashMap<u64, HashMap<u64, (u64, HashMap<u64, (u64, u64)>)>>,
}

impl RIHists {
    pub fn new(
        ri_hists_input: HashMap<u64, HashMap<u64, (u64, HashMap<u64, (u64, u64)>)>>,
    ) -> Self {
        RIHists {
            ri_hists: ri_hists_input,
        }
    }

    pub fn get_ref_hist(&self, ref_id: u64) -> &HashMap<u64, (u64, HashMap<u64, (u64, u64)>)> {
        return self.ri_hists.get(&ref_id).unwrap();
    }

    pub fn get_ref_ri_count(&self, ref_id: u64, ri: u64) -> u64 {
        return self.ri_hists.get(&ref_id).unwrap().get(&ri).unwrap().0;
    }

    pub fn get_ref_ri_cost(&self, ref_id: u64, ri: u64) -> &HashMap<u64, (u64, u64)> {
        return &self.ri_hists.get(&ref_id).unwrap().get(&ri).unwrap().1;
    }

    pub fn get_ref_ri_phase_cost(&self, ref_id: u64, ri: u64, phase: u64) -> (u64, u64) {
        return *self
            .ri_hists
            .get(&ref_id)
            .unwrap()
            .get(&ri)
            .unwrap()
            .1
            .get(&phase)
            .unwrap();
    }
}

#[derive(Debug, Copy, Clone)]
pub struct PPUC {
    ppuc: f64,
    lease: u64,
    old_lease: u64,
    ref_id: u64,
    new_hits: u64,
}
impl PartialOrd for PPUC {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.ppuc.partial_cmp(&other.ppuc)
    }
}
impl Ord for PPUC {
    fn cmp(&self, other: &Self) -> Ordering {
        other.ppuc.partial_cmp(&self.ppuc).unwrap()
    }
}
impl PartialEq for PPUC {
    fn eq(&self, other: &Self) -> bool {
        other.ppuc.eq(&self.ppuc)
    }
}
impl Eq for PPUC {}

pub fn process_sample_head_cost(
    ri_hists: &mut HashMap<u64, HashMap<u64, (u64, HashMap<u64, (u64, u64)>)>>,
    phase_id_ref: u64,
    ri: u64,
    use_time: u64,
    next_phase_tuple: (u64, u64),
) {
    let phase_id = (phase_id_ref & 0xFF000000) >> 24;

    //increment count
    let ref_hist = ri_hists
        .entry(phase_id_ref)
        .or_insert_with(|| HashMap::new());
    let ri_tuple = ref_hist.entry(ri).or_insert_with(|| (0, HashMap::new()));
    ri_tuple.0 += 1;

    //increment head costs
    let this_phase_head_cost = std::cmp::min(next_phase_tuple.0 - use_time, ri);
    let next_phase_head_cost =
        std::cmp::max(use_time as i64 + ri as i64 - next_phase_tuple.0 as i64, 0) as u64;

    ri_tuple.1.entry(phase_id).or_insert_with(|| (0, 0)).0 += this_phase_head_cost;

    if next_phase_head_cost > 0 {
        ri_tuple
            .1
            .entry(next_phase_tuple.1)
            .or_insert_with(|| (0, 0))
            .0 += next_phase_head_cost;
    }
}

pub fn process_sample_tail_cost(
    ri_hists: &mut HashMap<u64, HashMap<u64, (u64, HashMap<u64, (u64, u64)>)>>,
    phase_id_ref: u64,
    ri: u64,
    use_time: u64,
    next_phase_tuple: (u64, u64),
) {
    let phase_id = (phase_id_ref & 0xFF000000) >> 24;

    let ref_hist = ri_hists
        .entry(phase_id_ref)
        .or_insert_with(|| HashMap::new());

    //this heinous code exists so we can iterate through a HashMap while modifying it
    let ris: Vec<&u64> = ref_hist.keys().collect();
    let mut ris_keys: Vec<u64> = Vec::new();
    for ri_other in ris {
        ris_keys.push(*ri_other);
    }

    //increment tail costs
    for ri_other in ris_keys {
        //no tail cost if the other ri is greater
        if ri_other >= ri {
            continue;
        }
        let count_phase_cost_tuple = ref_hist
            .entry(ri_other)
            .or_insert_with(|| (0, HashMap::new()));

        let this_phase_tail_cost = std::cmp::min(next_phase_tuple.0 - use_time, ri_other);
        let next_phase_tail_cost = std::cmp::max(
            0,
            use_time as i64 + ri_other as i64 - next_phase_tuple.0 as i64,
        ) as u64;

        count_phase_cost_tuple
            .1
            .entry(phase_id)
            .or_insert_with(|| (0, 0))
            .1 += this_phase_tail_cost;

        if next_phase_tail_cost > 0 {
            count_phase_cost_tuple
                .1
                .entry(next_phase_tuple.1)
                .or_insert_with(|| (0, 0))
                .1 += next_phase_tail_cost;
        }
    }
}

fn cshel_phase_ref_cost(
    sample_rate: u64,
    phase: u64,
    ref_id: u64,
    old_lease: u64,
    new_lease: u64,
    ri_hists: &RIHists,
) -> u64 {
    let mut old_cost = 0;
    let mut new_cost = 0;
    //if a set in a phase does not have this reference
    if ri_hists.ri_hists.get(&ref_id) == None {
        return 0;
    }
    let ri_hist = ri_hists.ri_hists.get(&ref_id).unwrap();

    for (&ri, (_, phase_cost_hashmap)) in ri_hist.iter() {
        let (phase_head_cost, phase_tail_cost) = match phase_cost_hashmap.get(&phase) {
            Some((a, b)) => (*a, *b),
            None => (0, 0),
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

fn shel_phase_ref_cost(
    sample_rate: u64,
    phase: u64,
    ref_id: u64,
    old_lease: u64,
    new_lease: u64,
    ri_hists: &RIHists,
) -> u64 {
    //if a set in a phase does not have this reference
    if ri_hists.ri_hists.get(&ref_id) == None {
        return 0;
    }
    let ref_ri_hist: &HashMap<u64, (u64, HashMap<u64, (u64, u64)>)> =
        ri_hists.ri_hists.get(&ref_id).unwrap();
    let ri_hist: Vec<(u64, u64)> = ref_ri_hist.iter().map(|(k, v)| (*k, v.0)).collect();
    let mut old_cost = 0;
    let mut new_cost = 0;
    if phase != (ref_id & 0xFF000000) >> 24 {
        return 0;
    }
    for (ri, count) in ri_hist.iter() {
        if *ri <= old_lease {
            old_cost += *count * *ri;
        } else {
            old_cost += *count * old_lease;
        }
        if *ri <= new_lease {
            new_cost += *count * *ri;
        } else {
            new_cost += *count * new_lease;
        }
    }
    (new_cost - old_cost) * sample_rate
}

pub fn get_ppuc(
    ref_id: u64,
    base_lease: u64,
    ref_ri_hist: &HashMap<u64, (u64, HashMap<u64, (u64, u64)>)>,
) -> Vec<PPUC> {
    let ri_hist: Vec<(u64, u64)> = ref_ri_hist.iter().map(|(k, v)| (*k, v.0)).collect();
    let total_count = ri_hist.iter().fold(0, |acc, (_k, v)| acc + v);
    let mut hits = 0;
    let mut head_cost = 0;

    let mut lease_hit_table = HashMap::new();
    let mut lease_cost_table = HashMap::new();

    //prevent kernel panic for a base lease that doesn't correspond to sampled ri for a reference
    lease_hit_table.insert(base_lease, 0);
    lease_cost_table.insert(base_lease, 0);

    let mut ri_hist_clone = ri_hist.clone();
    ri_hist_clone.sort_by(|a, b| a.0.cmp(&b.0));

    for (ri, count) in ri_hist_clone.iter() {
        hits += *count;
        head_cost += *count * *ri;
        let tail_cost = (total_count - hits) * *ri;

        lease_hit_table.insert(*ri, hits);
        lease_cost_table.insert(*ri, head_cost + tail_cost);
    }

    ri_hist_clone
        .iter()
        .map(|(k, _v)| k)
        .filter(|k| **k > base_lease)
        .map(|k| PPUC {
            ppuc: ((*lease_hit_table.get(k).unwrap() - *lease_hit_table.get(&base_lease).unwrap())
                as f64
                / (*lease_cost_table.get(k).unwrap() - *lease_cost_table.get(&base_lease).unwrap())
                    as f64),
            lease: *k,
            old_lease: base_lease,
            ref_id: ref_id,
            new_hits: *lease_hit_table.get(k).unwrap() - *lease_hit_table.get(&base_lease).unwrap(),
        })
        .collect()
}

pub fn get_avg_lease(distribution: &BinnedRIs, addr: &u64, bin: u64, lease: u64) -> u64 {
    let mut total = 0;
    for (ri, freq) in distribution
        .bin_ri_distribution
        .get(&bin)
        .unwrap()
        .get(&addr)
        .unwrap()
    {
        if *ri <= lease && *ri > 0 {
            total += ri * freq;
        } else {
            total += lease * freq;
        }
    }
    return total;
}

pub fn prl(
    bin_width: u64,
    ri_hists: &RIHists,
    binned_ris: &BinnedRIs,
    binned_freqs: &BinFreqs,
    sample_rate: u64,
    cache_size: u64,
    discretize: u64,
    samples_per_phase: &HashMap<u64, u64>,
    verbose: bool,
    debug: bool,
    set_mask: u32,
) -> Option<(
    HashMap<u64, u64>,
    HashMap<u64, (f64, u64)>,
    HashMap<u64, HashMap<u64, u64>>,
    u64,
)> {
    let mut new_lease: PPUC;
    let mut dual_leases: HashMap<u64, (f64, u64)> = HashMap::new(); //{ref_id, (alpha, long_lease)}
    let mut trace_length: u64 = 0;
    let mut bin_endpoints: Vec<u64> = Vec::new();
    let mut num_unsuitable: u64;
    let mut ppuc_tree = BinaryHeap::new();
    let mut impact_dict: HashMap<u64, HashMap<u64, f64>> = HashMap::new();
    let mut bin_ranks: HashMap<u64, f64> = HashMap::new();
    let mut sorted_bins: Vec<(u64, f64)> = Vec::new();
    let mut acceptable_ratio: f64;
    let mut neg_impact;
    let mut lease_hits = HashMap::new();
    let mut num_full_bins;
    let mut bin_saturation: HashMap<u64, HashMap<u64, f64>> = HashMap::new();
    let mut leases: HashMap<u64, u64> = HashMap::new();

    let num_sets = set_mask as u64 + 1;
    let bin_target: u64 = bin_width * cache_size / num_sets;
    //threshold for meaningful dual lease
    let min_alpha = 1.0
        - (((2 << (discretize - 1)) as f64) - 1.5 as f64)
            / (((2 << (discretize - 1)) as f64) - 1.0 as f64);

    if verbose {
        println!("---------Dump Binned RI Hists------------");
        super::io::debug::print_binned_hists(&binned_ris);
        println!("---------Dump Reference Frequency per bin---");
        println!("{:?}", &binned_freqs);
    }

    if verbose {
        println!("bin_width:  {}", bin_width);
    }
    for key in binned_freqs.bin_freqs.keys() {
        bin_endpoints.push(*key);
    }
    //each bin will have all addresses although freq may be 0
    let mut addrs: Vec<u64> = Vec::new();
    for key in binned_freqs.bin_freqs.get(&0).unwrap().keys() {
        addrs.push(*key);
    }

    for endpoint in bin_endpoints.iter() {
        for set in 0..num_sets {
            bin_saturation
                .entry(*endpoint)
                .or_insert(HashMap::new())
                .entry(set)
                .or_insert(0.0);
        }
    }
    //make all references have lease of 1
    for addr in addrs {
        leases.insert(addr & 0x00FFFFFF, 1);
        // update saturation to take into account each reference having a lease of 1
        for (bin, _sat) in bin_saturation.clone() {
            for set in 0..num_sets {
                if binned_ris
                    .bin_ri_distribution
                    .get(&bin)
                    .unwrap()
                    .contains_key(&addr)
                {
                    let old_avg_lease = get_avg_lease(binned_ris, &addr, bin, 0);
                    let avg_lease = get_avg_lease(binned_ris, &addr, bin, 1);
                    let impact = (avg_lease as f64 - old_avg_lease as f64) * &(sample_rate as f64);
                    let bin_saturation_set = bin_saturation.get_mut(&bin).unwrap();
                    bin_saturation_set.insert(set, bin_saturation_set.get(&set).unwrap() + impact);
                }
                //init impact dict for later
                impact_dict
                    .entry(bin)
                    .or_insert(HashMap::new())
                    .entry(set)
                    .or_insert(0.0);
            }
        }
    }

    for (_phase, &num) in samples_per_phase.iter() {
        trace_length += num * sample_rate;
    }

    for (&ref_id, ri_hist) in ri_hists.ri_hists.iter() {
        let ppuc_vec = get_ppuc(ref_id, 0, ri_hist);
        for ppuc in ppuc_vec.iter() {
            ppuc_tree.push(*ppuc);
        }
    }
    // get lease hits assuming a base lease of 0
    for _r in ppuc_tree.clone() {
        let lease = ppuc_tree.pop().unwrap();
        lease_hits
            .entry(lease.ref_id)
            .or_insert(HashMap::new())
            .entry(lease.lease)
            .or_insert(lease.new_hits);
    }
    //reinitallize ppuc tree, assuming a base lease of 1
    for (&ref_id, ri_hist) in ri_hists.ri_hists.iter() {
        let ppuc_vec = get_ppuc(ref_id, 1, ri_hist);
        for ppuc in ppuc_vec.iter() {
            ppuc_tree.push(*ppuc);
        }
    }

    loop {
        new_lease = match ppuc_tree.pop() {
            Some(i) => i,
            None => return Some((leases, dual_leases, lease_hits, trace_length)),
        };

        //continue to pop until we have a ppuc with the right base_lease
        if new_lease.old_lease != *leases.get(&(new_lease.ref_id & 0xFFFFFFFF)).unwrap() {
            continue;
        }
        //won't assign a reference to a reference that has recieved a dual lease
        if dual_leases.contains_key(&(new_lease.ref_id & 0xFFFFFFFF)) {
            continue;
        }
        neg_impact = false;
        num_unsuitable = 0;
        let addr = new_lease.ref_id;
        for (bin, bin_sat_set) in &bin_saturation {
            for (set, _sat) in bin_sat_set {
                let mut impact: f64 = 0.0;
                let set_addr = (addr & 0xFFFFFFFF) | (set << 32);
                if binned_ris
                    .bin_ri_distribution
                    .get(&bin)
                    .unwrap()
                    .contains_key(&set_addr)
                {
                    let old_avg_lease = get_avg_lease(
                        binned_ris,
                        &set_addr,
                        *bin,
                        *leases.get(&(addr & 0xFFFFFFFF)).unwrap(),
                    );
                    let avg_lease = get_avg_lease(binned_ris, &set_addr, *bin, new_lease.lease);
                    impact = (avg_lease as f64 - old_avg_lease as f64) * &(sample_rate as f64);
                    //don't assign leases that decrease bin saturation
                    neg_impact = if impact >= 0.0 { false } else { true };
                    impact_dict
                        .get_mut(&bin)
                        .unwrap()
                        .insert(*set, impact as f64);
                } else {
                    impact_dict.get_mut(&bin).unwrap().insert(*set, 0 as f64);
                }
                if (bin_saturation.get(&bin).unwrap().get(&set).unwrap() + impact)
                    > bin_target as f64
                {
                    num_unsuitable += 1;
                }
            }
        }
        for (bin, sat_set) in &bin_saturation.clone() {
            for (set, sat) in sat_set {
                if verbose && debug {
                    println!(
                        "bin: {} set: {} current capacity: {:.7} impact: {:.7}",
                        bin / bin_width,
                        set,
                        sat / bin_width as f64,
                        impact_dict.get(bin).unwrap().get(set).unwrap() / bin_width as f64
                    );
                }
            }
        }
        if verbose {
            println!("addr:{:x} ri:{:x}", addr & 0xFFFFFFFF, new_lease.lease);
        }
        //skip lease, if it makes it worse
        if neg_impact {
            if verbose {
                println!("lease value would be negative");
            }
            continue;
        }
        if num_unsuitable < 1 {
            leases.insert(addr & 0xFFFFFFFF, new_lease.lease);
            //push new ppucs
            let ppuc_vec = get_ppuc(
                new_lease.ref_id,
                new_lease.lease,
                ri_hists.ri_hists.get(&new_lease.ref_id).unwrap(),
            );

            for ppuc in ppuc_vec.iter() {
                ppuc_tree.push(*ppuc);
            }
            let mut print_string: String = String::new();
            for (bin, sat_set) in &bin_saturation.clone() {
                for (set, sat) in sat_set {
                    let set_addr = (addr & 0xFFFFFFFF) | (set << 32);
                    if binned_ris
                        .bin_ri_distribution
                        .get(bin)
                        .unwrap()
                        .contains_key(&set_addr)
                    {
                        bin_saturation
                            .get_mut(bin)
                            .unwrap()
                            .insert(*set, sat + impact_dict.get(bin).unwrap().get(set).unwrap());
                    }
                    print_string = format!(
                        "{:} {1:.7}",
                        print_string,
                        &(bin_saturation.get(bin).unwrap().get(set).unwrap() / bin_width as f64)
                    );
                }
            }
            if verbose {
                println!(
                    "assigning lease: {:x} to reference {:x}",
                    new_lease.lease, addr
                );
                println!("Average cache occupancy per bin: [{}]", print_string);
            }
        } else {
            num_full_bins = 0;
            let mut alpha = 1.0;
            for (bin, sat_set) in &bin_saturation {
                for (set, sat) in sat_set {
                    if sat >= &(bin_target as f64) {
                        num_full_bins += 1
                    }
                    let new_capacity = sat + impact_dict.get(bin).unwrap().get(set).unwrap();

                    if &new_capacity >= &(bin_target as f64) {
                        if *impact_dict.get(bin).unwrap().get(set).unwrap() != 0.0 {
                            //get minimum set alpha for bin
                            let set_alpha = ((bin_target as f64) - sat)
                                / *impact_dict.get(bin).unwrap().get(set).unwrap();
                            if set_alpha < alpha {
                                alpha = set_alpha;
                            }
                        }
                    }
                }
                bin_ranks.insert(*bin, alpha);
            }
            for (num, key_val_pair) in bin_ranks.iter().enumerate() {
                if sorted_bins.get_mut(num) == None {
                    sorted_bins.push((*key_val_pair.0, *key_val_pair.1));
                } else {
                    sorted_bins[num] = (*key_val_pair.0, *key_val_pair.1);
                }
            }
            sorted_bins.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
            acceptable_ratio = if num_full_bins == 0 {
                sorted_bins[0].1
            } else {
                0.0
            };

            if acceptable_ratio > min_alpha {
                dual_leases.insert(addr, (acceptable_ratio as f64, new_lease.lease));
                let mut print_string: String = String::new();

                for (bin, sat_set) in &bin_saturation.clone() {
                    for (set, sat) in sat_set {
                        let set_addr = (addr & 0xFFFFFFFF) | (set << 32);
                        if binned_ris
                            .bin_ri_distribution
                            .get(bin)
                            .unwrap()
                            .contains_key(&set_addr)
                        {
                            bin_saturation.get_mut(bin).unwrap().insert(
                                *set,
                                sat + impact_dict.get(bin).unwrap().get(set).unwrap()
                                    * acceptable_ratio,
                            );
                        }
                        print_string = format!(
                            "{:} {1:.7}",
                            print_string,
                            &(bin_saturation.get(bin).unwrap().get(set).unwrap()
                                / bin_width as f64)
                        );
                    }
                }
                if verbose {
                    println!(
                        "Assigning dual lease {} to address {} with percentage: {}",
                        new_lease.lease, addr, acceptable_ratio
                    );

                    println!("Average cache occupancy per bin: [{}]", print_string);
                }
            }
        }
    }
}

pub fn get_num_leases_per_phase(leases: &HashMap<u64, u64>) -> HashMap<u64, u64> {
    let mut references_per_phase: HashMap<u64, u64> = HashMap::new();
    for (phase_id_x_reference, _lease) in leases.iter() {
        let phase_id = (phase_id_x_reference & 0xFF000000) >> 24;
        references_per_phase
            .entry(phase_id)
            .and_modify(|e| *e += 1)
            .or_insert(1);
    }
    return references_per_phase;
}

pub fn prune_leases_to_fit_llt(
    leases: HashMap<u64, u64>,
    dual_leases: HashMap<u64, (f64, u64)>,
    ri_hists: &RIHists,
    llt_size: u64,
) -> (HashMap<u64, u64>, HashMap<u64, (f64, u64)>) {
    let mut pruned_leases: HashMap<u64, u64> = HashMap::new();
    let mut pruned_dual_leases: HashMap<u64, (f64, u64)> = HashMap::new();
    let references_per_phase: HashMap<u64, u64> = get_num_leases_per_phase(&leases);

    for (phase_id, _lease_count) in references_per_phase.iter() {
        //loop through phases
        let mut importance_per_reference: HashMap<u64, u64> = HashMap::new();

        //this is globally sorting leases by importance
        //need to be locally sorting them per phase
        for (reference, _lease) in leases.iter() {
            let reference_phase_id = (reference & 0xFF000000) >> 24;
            //if this reference is not in the current phase, pass instead of inserting
            if reference_phase_id != *phase_id {
                continue;
            }
            let ri_hist = ri_hists.get_ref_hist(*reference);
            let mut count = 0;
            //need to sum over this
            for (_ri, count_cost_tuple) in ri_hist {
                count += count_cost_tuple.0;
            }
            importance_per_reference.entry(*reference).or_insert(count);
        }

        let mut importance_vec: Vec<_> = importance_per_reference.iter().collect();
        importance_vec.sort_by(|a, b| a.1.cmp(b.1).reverse());

        for i in 0..llt_size {
            if i == importance_vec.len() as u64 {
                break;
            }

            //add the top llt_size most important leases to the pruned vector
            let reference_id = importance_vec[i as usize].0;
            pruned_leases
                .entry(*reference_id)
                .or_insert(*leases.get(reference_id).unwrap());

            if dual_leases.contains_key(reference_id) {
                pruned_dual_leases
                    .entry(*reference_id)
                    .or_insert(*dual_leases.get(reference_id).unwrap());
            }
            //println!("Inserted successfully");
        }
    }
    return (pruned_leases, pruned_dual_leases);
}

//Output:
//leases: Hashmap<u64,u64>
//dual_leases: HashMap<u64, (f64, u64)>
//lease_hits: HashMap<u64, HashMap<u64,u64>>
//trace_length: u64
pub fn shel_cshel(
    cshel: bool,
    ri_hists: &RIHists,
    cache_size: u64,
    sample_rate: u64,
    samples_per_phase: &HashMap<u64, u64>,
    discretize: u64,
    verbose: bool,
    debug: bool,
    set_mask: u32,
) -> Option<(
    HashMap<u64, u64>,
    HashMap<u64, (f64, u64)>,
    HashMap<u64, HashMap<u64, u64>>,
    u64,
)> {
    let mut new_lease: PPUC;
    let mut cost_per_phase: HashMap<u64, HashMap<u64, u64>> = HashMap::new();
    let mut budget_per_phase: HashMap<u64, u64> = HashMap::new();
    let mut leases = HashMap::new(); //{ri, lease}
    let mut dual_leases: HashMap<u64, (f64, u64)> = HashMap::new(); //{ref_id, (alpha, long_lease)}
    let mut trace_length: u64 = 0;
    let mut lease_hits = HashMap::new();
    let mut dual_lease_phases: Vec<u64> = Vec::new();
    //{phase,(cost with alpha, cost if alpha was 1, ref ID)}
    let mut past_lease_values: HashMap<u64, (u64, u64)> = HashMap::new();
    let mut last_lease_cost: HashMap<u64, HashMap<u64, (u64, u64, u64)>> = HashMap::new();

    let num_sets = set_mask as u64 + 1;
    let phase_ids: Vec<&u64> = samples_per_phase.keys().collect();

    //since we can't run CSHEL without also running SHEL, don't output RI history twice
    if !cshel {
        if verbose {
            println!("---------Dump RI Hists------------");
            super::io::debug::print_ri_hists(&ri_hists);
            println!("---------Dump Samples Per Phase---");
            println!("{:?}", &samples_per_phase);
        }
    }

    //threshold for meaningful dual lease
    let min_alpha: f64 = 1.0
        - (((2 << (discretize - 1)) as f64) - 1.5 as f64)
            / (((2 << (discretize - 1)) as f64) - 1.0 as f64);
    //initialize ppucs
    let mut ppuc_tree = BinaryHeap::new();

    for (&ref_id, ri_hist) in ri_hists.ri_hists.iter() {
        let ppuc_vec = get_ppuc(ref_id, 0, ri_hist);
        for ppuc in ppuc_vec.iter() {
            ppuc_tree.push(*ppuc);
        }
    }

    // get lease hits assuming a base lease of 0
    for _r in ppuc_tree.clone() {
        let lease = ppuc_tree.pop().unwrap();
        //sum hits for reference over all sets
        *lease_hits
            .entry(lease.ref_id)
            .or_insert(HashMap::new())
            .entry(lease.lease)
            .or_insert(0) += lease.new_hits;
    }

    //reinitalize ppuc tree, assuming a base lease of 1
    for (&ref_id, ri_hist) in ri_hists.ri_hists.iter() {
        let ppuc_vec = get_ppuc(ref_id, 1, ri_hist);
        for ppuc in ppuc_vec.iter() {
            ppuc_tree.push(*ppuc);
        }
    }

    //initialize cost + budget
    for (&phase, &num) in samples_per_phase.iter() {
        budget_per_phase
            .entry(phase)
            .or_insert(num * cache_size / num_sets * sample_rate);
        trace_length += num * sample_rate;
    }

    if verbose {
        println!(
            "
        ---------------------
        Initial budget per phase: 
        {:?} 
        ---------------------",
            budget_per_phase
        );
    }
    //initialize leases to a default value of 1
    for (&ref_id, _) in ri_hists.ri_hists.iter() {
        leases.insert(ref_id & 0xFFFFFFFF, 1);
        let phase = (ref_id & 0xFF000000) >> 24;
        let phase_id_ref = ref_id & 0xFFFFFFFF;
        // get cost of assigning a lease of 1 for each set
        for set in 0..num_sets {
            let set_phase_id_ref = phase_id_ref | (set << 32);
            let new_cost = match cshel {
                true => cshel_phase_ref_cost(sample_rate, phase, set_phase_id_ref, 0, 1, &ri_hists),
                false => shel_phase_ref_cost(sample_rate, phase, set_phase_id_ref, 0, 1, &ri_hists),
            };
            *cost_per_phase
                .entry(phase)
                .or_insert(HashMap::new())
                .entry(set)
                .or_insert(0) += new_cost;
        }
    }
    if verbose {
        println!("costs per phase{:?}", cost_per_phase);
    }

    loop {
        new_lease = match ppuc_tree.pop() {
            //TERMINATION CONDITION 1
            Some(i) => i,
            None => return Some((leases, dual_leases, lease_hits, trace_length)),
        };
        let phase = (new_lease.ref_id & 0xFFFFFFFF) >> 24;
        let ref_id = new_lease.ref_id & 0xFFFFFFFF;

        //continue to pop until we have a ppuc with the right base_lease
        if new_lease.old_lease != *leases.get(&ref_id).unwrap() {
            continue;
        }

        let mut set_full = false;
        for set in 0..num_sets {
            if cost_per_phase.get(&phase).unwrap().get(&set).unwrap()
                == budget_per_phase.get(&phase).unwrap()
            {
                set_full = true;
                break;
            }
        }

        //if any set in phase is full, skip
        if set_full {
            continue;
        }
        //if we've already assigned dual leases to all phases, end
        if dual_lease_phases.len() == cost_per_phase.len() {
            //TERMINATION CONDITION 2
            return Some((leases, dual_leases, lease_hits, trace_length));
        }
        //if we've already assigned a dual lease for the phase
        if dual_lease_phases.contains(&phase) {
            continue;
        }

        let old_lease = *leases.get(&ref_id).unwrap();
        //check for capacity
        let mut acceptable_lease = true;
        let mut new_phase_ref_cost: HashMap<u64, HashMap<u64, u64>> = HashMap::new();
        for (&phase, current_cost) in cost_per_phase.iter() {
            //get cost of assigning a lease of 1 for each set
            for set in 0..num_sets {
                let set_phase_id_ref = ref_id | (set << 32);
                let new_cost = match cshel {
                    true => cshel_phase_ref_cost(
                        sample_rate,
                        phase,
                        set_phase_id_ref,
                        old_lease,
                        new_lease.lease,
                        &ri_hists,
                    ),
                    false => shel_phase_ref_cost(
                        sample_rate,
                        phase,
                        set_phase_id_ref,
                        old_lease,
                        new_lease.lease,
                        &ri_hists,
                    ),
                };

                new_phase_ref_cost
                    .entry(phase)
                    .or_insert(HashMap::new())
                    .entry(set)
                    .or_insert(new_cost);
                if (new_cost + current_cost.get(&set).unwrap())
                    > *budget_per_phase.get(&phase).unwrap()
                {
                    acceptable_lease = false;
                }
            }
        }
        if verbose & debug {
            println!("\nDebug: budgets per phase {:?}", &budget_per_phase);
            println!("Debug: Current cost budgets {:?}", &cost_per_phase);
            println!("Debug: NEW_PHASE_REF_COST {:?}", &new_phase_ref_cost);
        }
        if acceptable_lease {
            //update cache use
            for (phase, phase_set_costs) in cost_per_phase.iter_mut() {
                for (set, set_costs) in phase_set_costs.iter_mut() {
                    *set_costs += new_phase_ref_cost.get(phase).unwrap().get(set).unwrap();
                }
            }
            let phase = (new_lease.ref_id & 0xFF000000) >> 24;
            //store lease value we assign to the reference and
            //the value of the previously assigned lease for that reference
            past_lease_values.insert(
                new_lease.ref_id & 0xFFFFFFFF,
                (
                    new_lease.lease,
                    *leases.get(&(&new_lease.ref_id & 0xFFFFFFFF)).unwrap(),
                ),
            );

            if last_lease_cost.get_mut(&phase) == None {
                last_lease_cost.insert(phase, HashMap::new());
            }
            for set in 0..num_sets {
                last_lease_cost.get_mut(&phase).unwrap().insert(
                    set,
                    (
                        *new_phase_ref_cost.get(&phase).unwrap().get(&set).unwrap(),
                        *new_phase_ref_cost.get(&phase).unwrap().get(&set).unwrap(),
                        ref_id & 0xFFFFFFFF,
                    ),
                );
            }
            //update leases
            leases.insert(new_lease.ref_id & 0xFFFFFFFF, new_lease.lease);
            //push new ppucs
            let ppuc_vec = get_ppuc(
                new_lease.ref_id,
                new_lease.lease,
                ri_hists.ri_hists.get(&new_lease.ref_id).unwrap(),
            );

            for ppuc in ppuc_vec.iter() {
                ppuc_tree.push(*ppuc);
            }
            if verbose {
                println!(
                    "Assigned lease {:x} to reference ({},{:x}).",
                    new_lease.lease,
                    (new_lease.ref_id & 0xFF000000) >> 24,
                    new_lease.ref_id & 0x00FFFFFF
                );
            }
        } else {
            //unacceptable lease, must assign a dual lease
            let mut alpha = 1.0;
            let mut current_phase_alpha = 1.0;
            for (&phase, phase_set_current_cost) in cost_per_phase.iter() {
                let set_budget = *budget_per_phase.get(&phase).unwrap();
                for (&set, &current_set_cost) in phase_set_current_cost.iter() {
                    let &set_phase_ref_cost =
                        new_phase_ref_cost.get(&phase).unwrap().get(&set).unwrap();
                    if set_phase_ref_cost > 0 {
                        if set_budget < current_set_cost {
                            println!(
                                "
                            ERROR: current cost exceeds budget
                            *budget_per_phase.get(&phase)=.unwrap():  {}
                            currenc_cost:                            {}
                            ",
                                set_budget, current_set_cost
                            );
                            panic!();
                        }

                        let remaining_budget = set_budget - current_set_cost;
                        //get the best alpha for any set  (ignoring other phases) that we want for the current reference
                        if phase == (new_lease.ref_id & 0xFF000000) >> 24 {
                            current_phase_alpha = super::helpers::float_min(
                                current_phase_alpha,
                                remaining_budget as f64 / set_phase_ref_cost as f64,
                            );
                        }
                        alpha = super::helpers::float_min(
                            alpha,
                            remaining_budget as f64 / set_phase_ref_cost as f64,
                        );
                    }
                }
            }
            //if the alpha we wish to assign would result in
            //a long lease that is never used because the short lease
            //probabiliy will be 1 after descretizing, don't assign dual lease.
            if current_phase_alpha < min_alpha {
                println!("Assigning lease {:x} with percentage {} to reference ({},{:x}) would not be meaningful.", 
                             new_lease.lease,current_phase_alpha,(new_lease.ref_id & 0xFF000000) >> 24,
                             new_lease.ref_id & 0x00FFFFFF);
                continue;
            }

            if alpha > min_alpha {
                //update cache use
                for (phase, phase_set_costs) in cost_per_phase.iter_mut() {
                    let mut set_budget = *budget_per_phase.get(phase).unwrap();
                    for (set, set_costs) in phase_set_costs.iter_mut() {
                        *set_costs += (*new_phase_ref_cost.get(phase).unwrap().get(set).unwrap()
                            as f64
                            * alpha)
                            .round() as u64;
                        //fix floating point precision error leading
                        //to "overallocation" or underallocation
                        if set_costs > &mut set_budget {
                            *set_costs = set_budget;
                        }
                    }
                }
            }

            if cshel {
                //if there's no alpha that would assign a meaningful dual lease
                //that wouldn't put other phases over budget
                if alpha <= min_alpha {
                    let mut new_costs = HashMap::new();
                    let mut new_alpha = HashMap::new();
                    let mut adjust_lease = true;
                    let mut phase_alpha = 1.0;
                    for phase in &phase_ids {
                        for set in 0..num_sets {
                            let set_phase_ref_cost =
                                new_phase_ref_cost.get(&phase).unwrap().get(&set).unwrap();
                            //if the phase would be effected by the lease assignment
                            if set_phase_ref_cost > &0 {
                                //get phases that would be over budgeted by assigning the current lease.
                                //then subtract the cost of their prior dual lease (which may be, due to the default, a non-dual lease)
                                //and then add the spillover cost from the new leases
                                let past_cost_actual = if last_lease_cost.get(&phase) == None {
                                    0
                                } else if last_lease_cost.get(&phase).unwrap().get(&set) == None {
                                    0
                                } else {
                                    last_lease_cost.get(&phase).unwrap().get(&set).unwrap().0
                                };
                                if new_costs.get(&phase) == None {
                                    new_costs.insert(phase, HashMap::new());
                                }
                                let new_cost =
                                    cost_per_phase.get(&phase).unwrap().get(&set).unwrap()
                                        - past_cost_actual
                                        + (*set_phase_ref_cost as f64 * current_phase_alpha).round()
                                            as u64;
                                new_costs.get_mut(&phase).unwrap().insert(set, new_cost);
                                //if no lease adjustment can be made to keep the phase from being over budget
                                if new_costs.get(&phase).unwrap().get(&set).unwrap()
                                    > budget_per_phase.get(&phase).unwrap()
                                {
                                    adjust_lease = false;
                                    break;
                                }
                                let remaining_budget = *budget_per_phase.get(&phase).unwrap()
                                    - new_costs.get(&phase).unwrap().get(&set).unwrap();
                                //if cost of last lease was zero i.e., no prior lease for phase, then alpha will be 1 and will not be adjusted
                                let past_cost_max = if past_cost_actual != 0 {
                                    last_lease_cost.get(phase).unwrap().get(&set).unwrap().1
                                } else {
                                    0
                                };
                                if past_cost_max != 0 {
                                    //if previous long lease didn't fill phase, could be greater than one
                                    let set_phase_alpha = super::helpers::float_min(
                                        1.0,
                                        remaining_budget as f64 / past_cost_max as f64,
                                    );
                                    if set_phase_alpha <= min_alpha {
                                        let old_phase_ref = last_lease_cost
                                            .get(phase)
                                            .unwrap()
                                            .get(&set)
                                            .unwrap()
                                            .2;
                                        dual_leases.get(&old_phase_ref).unwrap().1;
                                        println!("Assigning adjusted dual lease {:x} with percentage {} to reference ({},{:x}) would not be meaningful.", 
                                    new_lease.lease,set_phase_alpha,phase,old_phase_ref);
                                        adjust_lease = false;
                                        break;
                                    }
                                    //need the minimum alpha of any set in the phase
                                    else if set_phase_alpha < phase_alpha {
                                        phase_alpha = set_phase_alpha;
                                    }

                                    new_alpha.insert(phase, phase_alpha);
                                }
                            }
                            //new costs is equal to old cost
                            else {
                                if new_costs.get(&phase) == None {
                                    new_costs.insert(&phase, HashMap::new());
                                }
                                new_costs.get_mut(&phase).unwrap().insert(
                                    set,
                                    *cost_per_phase.get(&phase).unwrap().get(&set).unwrap(),
                                );
                            }
                        }
                    }
                    if adjust_lease == true {
                        for phase in &phase_ids {
                            //if adjusting lease
                            for set in 0..num_sets {
                                if new_alpha.get(&phase) != None {
                                    let old_phase_cost_max =
                                        last_lease_cost.get(phase).unwrap().get(&set).unwrap().1;
                                    let old_phase_ref =
                                        last_lease_cost.get(phase).unwrap().get(&set).unwrap().2;
                                    let new_phase_cost = (old_phase_cost_max as f64
                                        * new_alpha.get(phase).unwrap())
                                        as u64;

                                    //if phase had a dual lease
                                    if dual_lease_phases.contains(phase) {
                                        dual_leases.insert(
                                            old_phase_ref,
                                            (
                                                *new_alpha.get(&phase).unwrap(),
                                                dual_leases.get(&old_phase_ref).unwrap().1,
                                            ),
                                        );
                                    }
                                    //if we are not currently assigning a dual lease to this phase
                                    //generate dual lease from the past two lease values of the last reference assigned in this phase
                                    else if **phase != new_lease.ref_id >> 24 {
                                        //set prior single lease as long lease value with new alpha
                                        dual_leases.insert(
                                            old_phase_ref,
                                            (
                                                *new_alpha.get(&phase).unwrap(),
                                                past_lease_values.get(&old_phase_ref).unwrap().0,
                                            ),
                                        );
                                        //set the lease two references back as the short lease value
                                        leases.insert(
                                            old_phase_ref,
                                            past_lease_values.get(&old_phase_ref).unwrap().1,
                                        );
                                        dual_lease_phases.push(**phase);
                                    }

                                    last_lease_cost.get_mut(&phase).unwrap().insert(
                                        set,
                                        (new_phase_cost, old_phase_cost_max, old_phase_ref),
                                    );
                                    //update phase costs
                                    cost_per_phase.get_mut(&phase).unwrap().insert(
                                        set,
                                        *new_costs.get(&phase).unwrap().get(&set).unwrap()
                                            + new_phase_cost,
                                    );
                                }
                                //if not adjusting the lease
                                else {
                                    //update phase costs
                                    cost_per_phase.get_mut(&phase).unwrap().insert(
                                        set,
                                        *new_costs.get(&phase).unwrap().get(&set).unwrap(),
                                    );
                                    //fix floating point precision error leading to "overallocation"
                                    if cost_per_phase.get(&phase).unwrap().get(&set).unwrap()
                                        > budget_per_phase.get(&phase).unwrap()
                                    {
                                        cost_per_phase
                                            .get_mut(&phase)
                                            .unwrap()
                                            .insert(set, *budget_per_phase.get(&phase).unwrap());
                                    }
                                }
                            }
                        }
                        alpha = current_phase_alpha;
                    } else {
                        //if we can't assign a dual lease without overflowing a phase
                        //without adjustment of past dual leases, with adjustment of past dual leases,
                        //or in the the unlikely case a phase is full with no dual lease

                        println!(
                            "Unable to assign lease {:x} with percentage {} to reference ({},{:x})",
                            new_lease.lease,
                            current_phase_alpha,
                            (new_lease.ref_id & 0xFF000000) >> 24,
                            new_lease.ref_id & 0x00FFFFFF
                        );
                        continue;
                    }
                }
            }

            let phase = (new_lease.ref_id & 0xFF000000) >> 24;

            //detect if set full
            let mut set_full = false;
            for set in 0..num_sets {
                if cost_per_phase.get(&phase).unwrap().get(&set).unwrap()
                    == budget_per_phase.get(&phase).unwrap()
                {
                    set_full = true;
                    break;
                }
            }
            //if last lease was a dual lease with alpha of 1 that didn't fill the budget, then it is actually a short lease and adjustments can be made to ensure
            //there is only 1 dual lease per phase.
            if alpha == 1.0 && set_full == false {
                //update leases
                leases.insert(new_lease.ref_id & 0xFFFFFFFF, new_lease.lease);
                //push new ppucs
                let ppuc_vec = get_ppuc(
                    new_lease.ref_id,
                    new_lease.lease,
                    ri_hists.ri_hists.get(&new_lease.ref_id).unwrap(),
                );

                for ppuc in ppuc_vec.iter() {
                    ppuc_tree.push(*ppuc);
                }
                if verbose {
                    println!(
                        "Assigned lease {:x} to reference ({},{:x}).",
                        new_lease.lease,
                        (new_lease.ref_id & 0xFF000000) >> 24,
                        new_lease.ref_id & 0x00FFFFFF
                    );
                }
            }
            //add dual lease
            else {
                //store cost of dual lease and store cost of lease with no dual lease and the reference for that lease
                for set in 0..num_sets {
                    if last_lease_cost.get_mut(&phase) == None {
                        last_lease_cost.entry(phase).or_insert(HashMap::new());
                    }

                    last_lease_cost.get_mut(&phase).unwrap().insert(
                        set,
                        (
                            (*new_phase_ref_cost.get(&phase).unwrap().get(&set).unwrap() as f64
                                * alpha)
                                .round() as u64,
                            *new_phase_ref_cost.get(&phase).unwrap().get(&set).unwrap(),
                            ref_id & 0xFFFFFFFF,
                        ),
                    );
                }

                dual_lease_phases.push(phase);
                //update dual lease HashMap
                dual_leases.insert(new_lease.ref_id & 0xFFFFFFFF, (alpha, new_lease.lease));

                if verbose {
                    println!(
                        "Assigned dual lease ({:x},{}) to reference ({},{:x}).",
                        new_lease.lease,
                        alpha,
                        (new_lease.ref_id & 0xFF000000) >> 24,
                        new_lease.ref_id & 0x00FFFFFF
                    );
                }
            }
        } //unacceptable lease

        if verbose & debug {
            for (phase, num) in samples_per_phase.iter() {
                for set in 0..num_sets {
                    println!(
                        "Debug: phase: {}. set: {} Allocation: {}",
                        phase,
                        set,
                        cost_per_phase.get(&phase).unwrap().get(&set).unwrap()
                            / (num * sample_rate)
                    );
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

        if verbose {
            let mut hits_from_old_lease = 0;

            if lease_hits.get(&new_lease.ref_id).unwrap().get(&old_lease) != None {
                hits_from_old_lease = *lease_hits
                    .get(&new_lease.ref_id)
                    .unwrap()
                    .get(&old_lease)
                    .unwrap();
            }
            let mut hits_from_new_lease = *lease_hits
                .get(&new_lease.ref_id)
                .unwrap()
                .get(&new_lease.lease)
                .unwrap();
            let long_lease_percentage: f64;
            if dual_leases.get(&new_lease.ref_id) != None {
                long_lease_percentage = dual_leases.get(&new_lease.ref_id).unwrap().0;
                let hits_without_dual = hits_from_new_lease;

                hits_from_new_lease = hits_without_dual
                    - (hits_without_dual as f64 * (1.0 - long_lease_percentage)) as u64
                    + ((1.0 - long_lease_percentage) * hits_from_old_lease as f64) as u64;
            }
            println!(
                "Additional hits from assigned lease:{}",
                (hits_from_new_lease - hits_from_old_lease) * sample_rate
            );
        }
    }
}
