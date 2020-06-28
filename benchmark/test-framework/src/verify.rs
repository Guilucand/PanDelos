use std::path::Path;
use std::fs::File;
use std::io::{BufReader, BufRead};
use std::collections::{BTreeSet, HashSet};
use std::cmp::Ordering;
use std::mem::swap;

pub struct Point {
    pub first: u64,
    pub second: u64,
    pub value: f32
}

impl PartialEq for Point {
    fn eq(&self, other: &Self) -> bool {
        self.first == other.first && self.second == other.second
    }
}

impl PartialOrd for Point {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.first == other.first  {
            self.second.partial_cmp(&other.second)
        }
        else {
            self.first.partial_cmp(&other.first)
        }
    }
}

impl Ord for Point {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

impl Eq for Point {
}


pub struct Network {
    pub points: BTreeSet<Point>
}

impl Network {
}

pub fn load_net(file: impl AsRef<Path>) -> Network {
    let reader = BufReader::new(File::open(file.as_ref()).unwrap());

    let mut result = Vec::new();

    for line in reader.lines() {

        let line = line.unwrap();

        let line: Vec<_> = line.split_whitespace().collect();

        let mut first = line[0].parse::<u64>().unwrap();
        let mut second = line[1].parse::<u64>().unwrap();
        let value = line[2].parse::<f32>().unwrap();

        if first > second {
            swap(&mut first, &mut second);
        }

        result.push(Point {
            first,
            second,
            value
        })
    }

    result.sort();
    result.dedup();
    Network {
        points: result.into_iter().collect()
    }
}