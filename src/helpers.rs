pub fn float_min(a: f64, b: f64) -> f64 {
    if a.lt(&b) {
        return a;
    }
    b
}

pub fn binary_search(vector: &Vec<(u64, u64)>, value: u64) -> Option<(u64, u64)> {
    let mut min = 0;
    let mut max = vector.len() - 1;

    while max >= min {
        let guess = (max + min) / 2;
        if vector[guess].0 == value {
            //on transition sample
            if guess < vector.len() - 1 {
                return Some(vector[guess + 1]);
            }
            //transition of last phase
            return None;
        }

        if vector[guess].0 < value {
            min = guess + 1;
        }
        if vector[guess].0 > value {
            if guess > 0 {
                max = guess - 1;
            } else {
                return Some(vector[guess]);
            }
        }
    }
    assert_eq!(max, min - 1);
    if min == vector.len() {
        return None;
    }

    Some(vector[min])
}
