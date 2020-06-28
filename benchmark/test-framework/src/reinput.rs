use std::path::{Path, PathBuf};
use std::fs::File;
use std::io::{BufReader, BufRead, BufWriter};
use std::collections::{HashSet, HashMap};
use std::io::Write;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher, BuildHasherDefault};
use std::collections;

#[derive(Clone)]
pub struct Fasta {
    header: String,
    sequence: String
}

type InputDetMap = HashMap<String, Vec<Fasta>, BuildHasherDefault<collections::hash_map::DefaultHasher>>;

pub struct ReInput {
    seqs: InputDetMap
}


pub fn reinput(
    input_file: impl AsRef<Path>,
    sequences_num: Option<usize>,
    mut genomes_num: Option<usize>) -> (impl AsRef<Path>, usize, usize) {

    let mut s = DefaultHasher::new();
    sequences_num.hash(&mut s);
    genomes_num.hash(&mut s);
    let hash = s.finish();


    let mut sequences_num = match sequences_num {
        None => usize::MAX,
        Some(x) => x,
    };

    let mut genomes_num = match genomes_num {
        None => usize::MAX,
        Some(x) => x,
    };

    let mut gencount = 0;
    let mut seqcount = 0;

    let input = ReInput::new(&input_file);
    let mut new_input = ReInput {
        seqs: HashMap::default()
    };

    for (gen, seq) in input.seqs {
        if genomes_num == 0 || sequences_num == 0 {
            break;
        }

        if seq.len() > sequences_num {
            new_input.seqs.insert(gen, seq[0..sequences_num].to_vec());
            seqcount += sequences_num;
            break;
        }
        else {
            sequences_num -= seq.len();
            seqcount += seq.len();
            new_input.seqs.insert(gen, seq);
        }

        genomes_num -= 1;
        gencount += 1;

    }

    sequences_num.hash(&mut s);
    genomes_num.hash(&mut s);

    let out_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("temp-inputs")
        .join(String::from(input_file.as_ref().file_name().unwrap().to_str().unwrap()) + &format!(".{:x}.faa", hash));

    new_input.write(&out_path);

    (out_path, seqcount, gencount)
}


impl ReInput {
    pub fn new(input: impl AsRef<Path>) -> ReInput {

        let reader = BufReader::new(File::open(input).unwrap());

        let mut seqs = InputDetMap::default();

        let mut lines = reader.lines();

        while let Some(line) = lines.next() {

            let header = line.unwrap();
            let sequence = lines.next().unwrap().unwrap();

            let genome = header.split_whitespace().next().unwrap();

            seqs.entry(genome.to_string()).or_insert(Vec::new()).push(Fasta {
                header,
                sequence
            })
        }

        ReInput {
            seqs
        }
    }

    pub fn write(&self, output: impl AsRef<Path>) {
        let mut output = BufWriter::new(File::create(output).unwrap());
        for (_genome, fastas) in &self.seqs {
            for fasta in fastas {
                writeln!(output, "{}", fasta.header);
                writeln!(output, "{}", fasta.sequence);
            }
        }
    }
}
