use serde::{Serialize, Deserialize};
use chrono::{DateTime, Utc};
use chrono::serde::ts_milliseconds;
use crate::CheckResults;

#[derive(Serialize, Deserialize, Debug)]
pub struct BenchmarkResults {
    real: f32,
    user: f32,
    system: f32,
    memory_mb: f32,
    #[serde(with = "ts_milliseconds")]
    exdate: DateTime<Utc>
}

impl BenchmarkResults {
    pub fn new(
        real: f32,
        user: f32,
        system: f32,
        memory_mb: f32,
    ) -> BenchmarkResults {
        BenchmarkResults {
            real,
            user,
            system,
            memory_mb,
            exdate: Utc::now()
        }
    }
}

#[derive(Deserialize, Serialize, Debug)]
pub struct BenchmarkReport {
    pub inpath: String,
    pub seqcount: usize,
    pub gencount: usize,
    pub new_bench: BenchmarkResults,
    pub old_bench: BenchmarkResults,
    pub check: CheckResults
}