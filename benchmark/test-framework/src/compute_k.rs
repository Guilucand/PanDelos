use std::collections::HashMap;
use std::path::Path;
use std::fs::File;
use std::io;
use std::io::BufRead;

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
    where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

pub fn compute_k(ifile: impl AsRef<Path>) -> u32 {

    let mut total_length = 0;
    let mut alphabet = HashMap::<char, usize>::new();

    let mut i = 0;

    for line in read_lines(ifile).unwrap() {
        let line = line.unwrap();
        if i % 2 != 0 {
            let trimmed = line.trim();
            total_length += trimmed.len();
            for s in trimmed.chars() {
                *alphabet.entry(s).or_insert(0) += 1;
            }
        }
        i += 1;
    }

    println!("total length {}", total_length);
    println!("alphabet {:?}", alphabet);
    println!("LG {}", (total_length as f64).log(alphabet.len() as f64));

    let mut k = 0.0;
    let size = alphabet.values().sum::<usize>() as f64;
    for count in alphabet.values() {
        let count = *count as f64;
        k += -(count / size).log(alphabet.len() as f64) * (count / size)
    }

    let uk = (total_length as f64).log(alphabet.len() as f64);
    let fk = uk / k;
    let k = fk.floor() as u32;

    println!("uk = {}", uk);
    println!("fk = {}", fk);
    println!("k = {}", k);
    k
}