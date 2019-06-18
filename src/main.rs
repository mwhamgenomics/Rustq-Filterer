extern crate structopt;
extern crate flate2;

use std::fs::File;
use structopt::StructOpt;
use flate2::read::GzDecoder;
use std::io::{BufRead,BufReader,Write};
use std::collections::HashSet;


#[derive(StructOpt)]
struct Cli {
    #[structopt(long="i1")]
    i1: std::path::PathBuf,

    #[structopt(long="i2")]
    i2: std::path::PathBuf,

    #[structopt(long="f1")]
    f1: std::path::PathBuf,

    #[structopt(long="f2")]
    f2: std::path::PathBuf,

    #[structopt(long="o1")]
    o1: std::path::PathBuf,

    #[structopt(long="o2")]
    o2: std::path::PathBuf,

    #[structopt(long="threshold", default_value="36")]
    len_threshold: usize,

    #[structopt(long="stats_file", parse(from_os_str))]
    stats_file: Option<std::path::PathBuf>,

    #[structopt(long="remove_tiles")]
    remove_tiles: Vec<String>,

    #[structopt(long="remove_reads", parse(from_os_str))]
    remove_reads: Option<std::path::PathBuf>,

    #[structopt(long="trim_r1")]
    trim_r1: Option<i32>,

    #[structopt(long="trim_r2")]
    trim_r2: Option<i32>
}


struct FastqEntry {
    read_id: String,
    seq: String,
    strand: String,
    qual: String,
    tile_id: String
}


impl FastqEntry {
    fn new() -> FastqEntry {
        FastqEntry {
            read_id: String::new(),
            seq: String::new(),
            strand: String::new(),
            qual: String::new(),
            tile_id: String::new()
        }
    }

    fn to_string(&self) {
        println!("read: {}seq: {}strand: {}qual: {}tile: {}", self.read_id, self.seq, self.strand, self.qual, self.tile_id);
    }

    fn output(&mut self, file: &mut File) {
        file.write(self.read_id.as_bytes());
        file.write(self.seq.as_bytes());
        file.write(self.strand.as_bytes());
        file.write(self.qual.as_bytes());
    }
}


struct FastqChecker {
    reader: BufReader<GzDecoder<File>>,
}


impl FastqChecker {
    fn new(fp: &std::path::PathBuf) -> FastqChecker {
        FastqChecker {
            reader: BufReader::new(GzDecoder::new(File::open(fp).unwrap())),
        }
    }

    fn read_entry(&mut self) -> Option<FastqEntry> {
        let mut entry = FastqEntry::new();
        self.reader.read_line(&mut entry.read_id);
        self.reader.read_line(&mut entry.seq);
        self.reader.read_line(&mut entry.strand);
        self.reader.read_line(&mut entry.qual);

        if entry.read_id.is_empty() ||
           entry.seq.is_empty() ||
           entry.strand.is_empty() ||
           entry.qual.is_empty() {
            None
        } else {
            let parts = &mut entry.read_id.split(":");
            entry.tile_id = parts.nth(4).unwrap().to_string();

            Some(entry)
        }
    }
}


struct RuntimeInfo<'a> {
    args: &'a Cli,
    rm_tiles: HashSet<String>,
    rm_reads: HashSet<String>,
    criteria: Vec<&'a Fn(&FastqEntry, &FastqEntry, &RuntimeInfo) -> bool>,
    read_pairs_checked: i64,
    read_pairs_removed: i64,
    read_pairs_remaining: i64,
}

impl<'a> RuntimeInfo <'a>{
    fn new(args: &Cli) -> RuntimeInfo {
        let mut rm_tiles = HashSet::new();
        let mut rm_reads = HashSet::new();

        let mut criteria: Vec<&Fn(&FastqEntry, &FastqEntry, &RuntimeInfo) -> bool> = vec![&check_read];

        if !args.remove_tiles.is_empty() {
            RuntimeInfo::build_rm_tiles(&args.remove_tiles, &mut rm_tiles);
            criteria.push(&tile_check_read);
        }

        match &args.remove_reads {
            Some(file_path) => {
                RuntimeInfo::build_rm_reads(file_path.to_path_buf(), &mut rm_reads);
                criteria.push(&id_check_read);
            },
            None => {}
        }

        RuntimeInfo {
            args: args,
            rm_tiles: rm_tiles,
            rm_reads: rm_reads,
            criteria: criteria,
            read_pairs_checked: 0,
            read_pairs_removed: 0,
            read_pairs_remaining: 0,
        }
    }

    fn build_rm_tiles(input_tiles: &Vec<String>, output_tiles: &mut HashSet<String>) {
        for t in input_tiles {
            output_tiles.insert(t.to_string());
        }
    }

    fn build_rm_reads(input_reads: std::path::PathBuf, output_reads: &mut HashSet<String>) -> std::io::Result<()> {
        let f = File::open(input_reads)?;
        let f = BufReader::new(f);
        for line in f.lines() {
            output_reads.insert(line.unwrap());
        }
        Ok(())
    }
}


fn check_read(entry_1: &FastqEntry, entry_2: &FastqEntry, info: &RuntimeInfo) -> bool {
    entry_1.seq.chars().count() > info.args.len_threshold && entry_2.seq.chars().count() > info.args.len_threshold
}


fn tile_check_read(entry_1: &FastqEntry, _entry_2: &FastqEntry, info: &RuntimeInfo) -> bool {
    let tiles = &info.rm_tiles;
    if tiles.contains(&entry_1.tile_id) {
        false
    } else {
        true
    }
}


fn id_check_read(entry_1: &FastqEntry, _entry_2: &FastqEntry, info: &RuntimeInfo) -> bool {
    let reads = &info.rm_reads;
    if reads.contains(&entry_1.read_id) {
        false
    } else {
        true
    }
}
   

fn check_reads(entry_1: &FastqEntry, entry_2: &FastqEntry, info: &RuntimeInfo) -> bool {
    let mut result: bool = true;
    for check_func in &info.criteria {
        if !check_func(entry_1, entry_2, &info) {
            result = false
        }
    }
    result
}


fn write_stats_file(fp: std::path::PathBuf, info: RuntimeInfo) -> std::io::Result<()> {
    let mut f = File::create(&fp)?;
    println!("Stats: {:?}", &fp);
    let report = format!(
        "r1i {:?}\nr1o {:?}\nr2i {:?}\nr2o {:?}\nr1f {:?}\nr2f {:?}\nread_pairs_checked {}\nread_pairs_removed {}\nread_pairs_remaining {}\n",
        info.args.i1, info.args.o1, info.args.i2, info.args.o2, info.args.f1, info.args.f2,
        info.read_pairs_checked, info.read_pairs_removed, info.read_pairs_remaining
    );
    f.write(report.as_bytes());
    Ok(())
}


fn main() -> std::io::Result<()> {
    let args = Cli::from_args();

    let mut info = RuntimeInfo::new(&args);

    let mut checker_1 = FastqChecker::new(&args.i1);
    let mut checker_2 = FastqChecker::new(&args.i2);

    let mut o1 = File::create(&args.o1)?;
    let mut o2 = File::create(&args.o2)?;

    let mut f1 = File::create(&args.f1)?;
    let mut f2 = File::create(&args.f2)?;

    loop {
        let entry_1 = checker_1.read_entry();
        let entry_2 = checker_2.read_entry();

        match (entry_1, entry_2) {
            (Some(mut e1), Some(mut e2)) => {
                info.read_pairs_checked += 1;
                if !check_reads(&e1, &e2, &info) {
                    info.read_pairs_removed += 1;
                    e1.output(&mut f1);
                    e2.output(&mut f2);
                } else {
                    info.read_pairs_remaining += 1;
                    e1.output(&mut o1);
                    e2.output(&mut o2);
                }
            }
            _ => {
                println!("Finished");
                break
            }
        }
    }

    match &info.args.stats_file {
        Some(file_path) => {
            write_stats_file(file_path.to_path_buf(), info);
        },
        None => {}
    }

    Ok(()) 
}

// cargo run -- --f1 this --f2 that --i1 other --i2 another --o1 more --o2 things

