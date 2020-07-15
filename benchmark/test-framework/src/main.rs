#![feature(vec_remove_item)]

mod manager;
mod compute_k;
mod bench;
mod verify;
mod reinput;

#[macro_use]
extern crate serde;

use structopt::StructOpt;
use crate::verify::load_net;
use crate::manager::build_output;
use crate::bench::{BenchmarkResults, BenchmarkReport};
use itertools::izip;
use std::fs::{File, OpenOptions};
use serde::Serialize;
use std::cmp::max;

#[derive(StructOpt, Debug, Clone)]
enum Args {
    Run {
        path: String,
        output_dir: String,
        test_files: Vec<String>,
        #[structopt(short, long)]
        genomes_max: Option<usize>,
        #[structopt(short, long)]
        sequences_max: Option<usize>,
        #[structopt(short, long)]
        nc: bool,
        #[structopt(short, long = "cmd-only")]
        cmd_only: bool
    },
    Check {
        first: String,
        second: String,
        #[structopt(short, long)]
        verbose: bool
    },
    Auto {
        input: Vec<String>,
        #[structopt(short, long)]
        genomes_max: Option<usize>,
        #[structopt(short, long)]
        sequences_max: Option<usize>,
        #[structopt(short, long)]
        jump_sequences: Option<usize>,
        #[structopt(short, long)]
        factor_jump_sequences: Option<f32>,
        #[structopt(short, long)]
        recompile: bool,
        #[structopt(short, long = "cmd-only")]
        cmd_only: bool
    }
}

#[derive(Serialize, Deserialize, Debug)]
pub struct CheckResults {
    equal: usize,
    missing_new: usize,
    missing_vanilla: usize,
    wrong_weight: usize
}

enum ProcessResults {
    None,
    Run(Vec<BenchmarkResults>),
    Check(CheckResults)
}

fn main() {
    let args = Args::from_args();
    process_args(args);
}

fn process_args(args: Args) -> ProcessResults {

    match args.clone() {
        Args::Run { path, output_dir, test_files, genomes_max, sequences_max, nc, cmd_only } => {

            if !nc && !cmd_only {
                manager::compile_pandelos(&path);
            }

            let mut res = Vec::new();

            for file in test_files {

                let (file, _seqc, _genc) = reinput::reinput(file,
                                            sequences_max, genomes_max);

                let k = compute_k::compute_k(&file);
                let bench = manager::execute_pandelos(
                    &path,
                    &file,
                    &output_dir,
                    k,
                    path.to_lowercase().contains("vanilla"),
                    cmd_only
                );
                res.push(bench);
            }
            ProcessResults::Run(res)
        },
        Args::Check { first, second, verbose } => {
            let a = load_net(first);
            let b = load_net(second);
            let mut missinga = 0;
            let mut missingb = 0;
            let mut present = 0;
            let mut diff_weight = 0;

            for x in &a.points {
                if !b.points.contains(&x) {
                    missinga += 1;
                    println!("MissingA {} <-> {} weight: {}", x.first, x.second, x.value);
                }
                else {

                    let val = b.points.get(&x).unwrap();
                    if (x.value - val.value).abs() > 0.001 {
                        println!("{} <-> {} = {} ~ {}", x.first, x.second, x.value, val.value);
                        diff_weight += 1;
                    }

                    present += 1;
                }
            }

            for x in &b.points {
                if !a.points.contains(&x) {
                    missingb += 1;
                    println!("MissingB {} <-> {} weight: {}", x.first, x.second, x.value);
                }
            }

            println!("Values {}+{}+({}) / {}", missinga, missingb, diff_weight, present);

            ProcessResults::Check(CheckResults {
                equal: present,
                missing_new: missingb,
                missing_vanilla: missinga,
                wrong_weight: diff_weight
            })
        }
        Args::Auto { mut input, mut genomes_max, mut sequences_max, jump_sequences, factor_jump_sequences, recompile, cmd_only } => {

            let mut work_done = true;

            while work_done {
                work_done = false;

                let bench = process_args(Args::Run {
                    path: dirs::home_dir().unwrap().join("PanDelos").to_str().unwrap().to_string(),
                    output_dir: "results/".to_string(),
                    test_files: input.clone(),
                    genomes_max,
                    sequences_max,
                    nc: !recompile,
                    cmd_only
                });

                let benchv = process_args(Args::Run {
                    path: dirs::home_dir().unwrap().join("PanDelos-vanilla").to_str().unwrap().to_string(),
                    output_dir: "results/".to_string(),
                    test_files: input.clone(),
                    genomes_max,
                    sequences_max,
                    nc: !recompile,
                    cmd_only
                });

                let bench = if let ProcessResults::Run(b) = bench { b } else { panic!("Error!") };
                let benchv = if let ProcessResults::Run(b) = benchv { b } else { panic!("Error!") };

                let logfile = OpenOptions::new()
                    .append(true).create(true).open("log.json").unwrap();


                let input_cpy = input.clone();

                for (origfile, b, bv) in izip!(input_cpy.into_iter(), bench.into_iter(), benchv.into_iter()) {

                    let (file, seqcount, gencount) = reinput::reinput(&origfile,
                                                                      sequences_max, genomes_max);

                    let dir = build_output("results/", &file, false);
                    let dir_vanilla = build_output("results/", &file, true);

                    if seqcount < sequences_max.unwrap_or(0) {
                        input.remove_item(&origfile).unwrap();
                    }
                    else {
                        work_done = true;
                    }


                    println!("AAA {} {}", dir.as_ref().display(), dir_vanilla.as_ref().display());

                    let checkres = process_args(Args::Check {
                        first: String::from(dir.as_ref().to_str().unwrap()),
                        second: String::from(dir_vanilla.as_ref().to_str().unwrap()),
                        verbose: false
                    });

                    let checkres = if let ProcessResults::Check(c) = checkres { c } else { panic!("Error!") };

                    let report = BenchmarkReport {
                        inpath: file.as_ref().to_str().unwrap().to_string(),
                        seqcount,
                        gencount,
                        new_bench: b,
                        old_bench: bv,
                        check: checkres
                    };

                    serde_json::to_writer(&logfile, &report);
                }
                match jump_sequences {
                    None => break,
                    Some(val) => {
                        match sequences_max {
                            None => break,
                            Some(seq) => {
                                if let Some(factor) = factor_jump_sequences {
                                    sequences_max = Some(max(seq + val, (seq as f32 * factor) as usize))
                                }
                                else {
                                    sequences_max = Some(seq + val)
                                }
                            },
                        }
                    },
                }
            }
            ProcessResults::None
        }
    }
}