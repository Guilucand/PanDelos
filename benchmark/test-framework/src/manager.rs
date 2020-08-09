use std::fs::File;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio, Child, exit};
use std::cell::UnsafeCell;
use std::time::Duration;
use std::io::{BufReader, BufRead, Read};
use std::os::raw::c_int;
use std::mem::MaybeUninit;
use crate::bench::BenchmarkResults;
use std::convert::TryInto;
use std::fmt::Arguments;

pub fn make_error(error: &str) -> ! {
    eprintln!("{}", error);
    std::process::exit(1);
}

pub trait UnwrapExitMessage {
    type Target;
    fn onerror_msg(self, message: &str) -> Self::Target;
    fn onerror_args(self, args: Arguments) -> Self::Target;
}

impl<T, U: ToString> UnwrapExitMessage for Result<T, U> {
    type Target = T;
    fn onerror_msg(self, message: &str) -> Self::Target {
        self.onerror_args(format_args!("{}", message))
    }

    fn onerror_args(self, args: Arguments) -> Self::Target {
        match self {
            Ok(x) => x,
            Err(err) => {
                eprintln!("Wargo command failed: {} [{}]!", args, err.to_string());
                exit(1);
            },
        }
    }
}



fn benchmark_command_get_output(name: &str,
                                args: &[&str],
                                dir: Option<&str>,
                                print: bool) -> (String, BenchmarkResults) {
    let mut command = Command::new(name);

    let start = std::time::Instant::now();

    println!("Calling pandelos: {}", String::from(name) + " " + &args.join(" "));



    let mut process = command.args(args)
        .current_dir(dir.map(|x| PathBuf::from(x)).unwrap_or(std::env::current_dir().unwrap()))
        .stderr(Stdio::inherit())
        .stdout(Stdio::piped())
        .spawn().unwrap();

    let pid = process.id();

    let mut stdout = Vec::new();
    match process.stdout.take() {
        None => {}
        (Some(mut out)) => {
            for line in BufReader::new(out).lines() {
                let line = line.unwrap();
                if print {
                    println!("{}", line);
                }
                stdout.extend_from_slice(line.as_str().as_bytes());
                stdout.push(b'\n');
            }
        }
    }

    let mut status = 0;
    let mut usage = unsafe { MaybeUninit::<libc::rusage>::zeroed().assume_init() };

    let elapsed;

    unsafe {
        let result = libc::wait4(pid as libc::pid_t, &mut status as *mut _, 0, &mut usage as *mut _);
        elapsed = start.elapsed();

        if result == -1 {
            make_error(&format!("Wait4 call failed process: {}!", name));
        }
    }

    let user = usage.ru_utime.tv_sec as f32 + usage.ru_utime.tv_usec as f32 / 1000000.0;
    let system = usage.ru_stime.tv_sec as f32 + usage.ru_stime.tv_usec as f32 / 1000000.0;
    let memory = usage.ru_maxrss as f32 / 1024.0;

    let results = BenchmarkResults::new(
        elapsed.as_secs_f32(), user, system, memory
    );

    println!("Results: {:#?}", results);

    // println!(r#"Rusage:
    // {{
    //     ru_real: {elapsed}s
    //     ru_utime: {ru_utime:?}s
    //     ru_stime: {ru_stime:?}s
    //     ru_maxrss: {ru_maxrss}MB
    //     ru_ixrss: {ru_ixrss}
    //     ru_idrss: {ru_idrss}
    //     ru_isrss: {ru_isrss}
    //     ru_minflt: {ru_minflt}
    //     ru_majflt: {ru_majflt}
    //     ru_nswap: {ru_nswap}
    //     ru_inblock: {ru_inblock}
    //     ru_oublock: {ru_oublock}
    //     ru_msgsnd: {ru_msgsnd}
    //     ru_msgrcv: {ru_msgrcv}
    //     ru_nsignals: {ru_nsignals}
    //     ru_nvcsw: {ru_nvcsw}
    //     ru_nivcsw: {ru_nivcsw}
    // }}
    // "#,
    // elapsed = elapsed.as_secs_f32(),
    // ru_utime = real,
    // ru_stime = user,
    // ru_maxrss = memory,
    // ru_ixrss = usage.ru_ixrss,
    // ru_idrss = usage.ru_idrss,
    // ru_isrss = usage.ru_isrss,
    // ru_minflt = usage.ru_minflt,
    // ru_majflt = usage.ru_majflt,
    // ru_nswap = usage.ru_nswap,
    // ru_inblock = usage.ru_inblock,
    // ru_oublock = usage.ru_oublock,
    // ru_msgsnd = usage.ru_msgsnd,
    // ru_msgrcv = usage.ru_msgrcv,
    // ru_nsignals = usage.ru_nsignals,
    // ru_nvcsw = usage.ru_nvcsw,
    // ru_nivcsw = usage.ru_nivcsw);

    // if !process.wait().unwrap().success() {
    //     make_error(&format!("'{}' command was not successful!", name));
    // }

    (String::from_utf8(stdout).unwrap(),
     results)
}

fn run_command(name: &Path, args: &[&str], dir: Option<PathBuf>) {
    let mut command = Command::new(name);
    let mut process = UnsafeCell::new(command.args(args)
        .current_dir(dir.unwrap_or(std::env::current_dir().unwrap()))
        .stderr(Stdio::inherit())
        .stdout(Stdio::inherit())
        .spawn()
        .unwrap());

    let status = unsafe { (*process.get()).wait().unwrap() };

    if !status.success() {
        make_error(&format!("'{:?}' command was not successful!", name));
    }
}

pub fn compile_pandelos<P: AsRef<Path>>(path: P) {



    let working_dir = if path.as_ref().join("pandelos-native").exists() {
        path.as_ref().join("pandelos-native")
    } else {
        path.as_ref().join("ig")
    };

    let exname = "bash";

    let wdir = working_dir.join("compile.sh");
    run_command(exname.as_ref(), &[wdir.to_str().unwrap()], Some(working_dir));
}

fn compute_hash(file: impl AsRef<Path>) -> u64 {
    let mut vec = Vec::new();
    File::open(&file).unwrap().read_to_end(&mut vec).unwrap();
    u64::from_be_bytes(md5::compute(vec.as_slice()).0[0..8].try_into().unwrap())
}

pub fn build_output(
        output_dir: impl AsRef<Path>,
        input: impl AsRef<Path>,
        vanilla: bool,
    ) -> impl AsRef<Path> {
    output_dir.as_ref().join(format!("{}{}-{:x}.out.net",
                                     input.as_ref().file_name().unwrap().to_str().unwrap(),
                                     if vanilla { "-vanilla" } else { "" }, compute_hash(&input)))
}

pub fn execute_pandelos(path: impl AsRef<Path>,
                        input: impl AsRef<Path>,
                        output_dir: impl AsRef<Path>,
                        k: u32, vanilla: bool,
                        print_only: bool) -> BenchmarkResults {

    let kvalue: String = k.to_string();
    let output_path = build_output(output_dir, &input, vanilla);

    // New version
    if path.as_ref().join("pandelos-native").exists() {
        let args = [
            "-i", input.as_ref().to_str().unwrap(),
            "-k", &kvalue,
            "-o", &output_path.as_ref().to_str().unwrap()
        ];

        let path = path.as_ref().join("pandelos-native/build/");

        return if !print_only {
            let (_output, bench) = benchmark_command_get_output(
                path.join("PanDelos").to_str().unwrap(), &args, None, true);
            bench
        } else {
            println!("{}/PanDelos {}", path.display(), args.join(" "));
            BenchmarkResults::new(
                0.0, 0.0, 0.0, 0.0
            )
        }
    }

    let args_new;
    let args_vanilla;
    let cp;
    let lpath;

    let args = if vanilla {
        cp = format!("{path}/ig/ext/commons-io-2.6.jar:{path}/ig/ig.jar",
                     path = path.as_ref().to_str().unwrap());

        args_vanilla = [
                "-Xmx2048g",
                "-cp",
                &cp,
                "infoasys.cli.pangenes.Pangenes",
                input.as_ref().to_str().unwrap(),
                &kvalue,
                &output_path.as_ref().to_str().unwrap()
        ];
        &args_vanilla[..]
    }  else {
        cp = format!("{path}/ig/ext/commons-io-2.6.jar:{path}/ig/ext/commons-cli-1.4.jar:{path}/ig/ig.jar",
                     path = path.as_ref().to_str().unwrap());

        lpath = format!("-Djava.library.path={path}/ig/native/build", path = path.as_ref().to_str().unwrap());

        args_new = [
            "-Xmx2048g",
            &lpath,
            "-cp",
            &cp,
            "infoasys.cli.pangenes.Pangenes",
            "-i", input.as_ref().to_str().unwrap(),
            "-k", &kvalue,
            "-o", &output_path.as_ref().to_str().unwrap()
        ];
        &args_new[..]
    };

    if !print_only {
        let (_output, bench) = benchmark_command_get_output("java".as_ref(), args, None, true);
        bench
    }
    else {
        println!("java {}", args.join(" "));
        BenchmarkResults::new(
            0.0, 0.0, 0.0, 0.0
        )
    }
}