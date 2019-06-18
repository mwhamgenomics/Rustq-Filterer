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


struct RuntimeParams<'a> {
    len_threshold: usize,
    rm_tiles: HashSet<String>,
    rm_reads: HashSet<String>,
    trim_r1: Option<i32>,
    trim_r2: Option<i32>,
    criteria: Vec<&'a Fn(&FastqEntry, &FastqEntry, &RuntimeParams) -> bool>
}

impl<'a> RuntimeParams <'a>{
    fn new(args: &Cli) -> RuntimeParams {
        let mut rm_tiles = HashSet::new();
        let mut rm_reads = HashSet::new();

        let mut criteria: Vec<&Fn(&FastqEntry, &FastqEntry, &RuntimeParams) -> bool> = vec![&check_read];

        if !args.remove_tiles.is_empty() {
            RuntimeParams::build_rm_tiles(&args.remove_tiles, &mut rm_tiles);
            criteria.push(&tile_check_read);
        }

        match &args.remove_reads {
            Some(file_path) => {
                RuntimeParams::build_rm_reads(file_path.to_path_buf(), &mut rm_reads);
                criteria.push(&id_check_read);
            },
            None => {}
        }

        RuntimeParams {
            len_threshold: args.len_threshold,
            rm_tiles: rm_tiles,
            rm_reads: rm_reads,
            trim_r1: args.trim_r1,
            trim_r2: args.trim_r2,
            criteria: criteria
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


fn check_read(entry_1: &FastqEntry, entry_2: &FastqEntry, params: &RuntimeParams) -> bool {
    entry_1.seq.chars().count() > params.len_threshold && entry_2.seq.chars().count() > params.len_threshold
}


fn tile_check_read(entry_1: &FastqEntry, _entry_2: &FastqEntry, params: &RuntimeParams) -> bool {
    let tiles = &params.rm_tiles;
    if tiles.contains(&entry_1.tile_id) {
        false
    } else {
        true
    }
}


fn id_check_read(entry_1: &FastqEntry, _entry_2: &FastqEntry, params: &RuntimeParams) -> bool {
    let reads = &params.rm_reads;
    if reads.contains(&entry_1.read_id) {
        false
    } else {
        true
    }
}
   

fn check_reads(entry_1: &FastqEntry, entry_2: &FastqEntry, params: &RuntimeParams) -> bool {
    let mut result: bool = true;
    for check_func in &params.criteria {
        if !check_func(entry_1, entry_2, &params) {
            result = false
        }
    }
    result
}



fn main() -> std::io::Result<()> {
    let args = Cli::from_args();

    let params = RuntimeParams::new(&args);

    let mut i1_reader = FastqChecker::new(&args.i1);
    let mut i2_reader = FastqChecker::new(&args.i2);

    let mut o1 = File::create(&args.o1)?;
    let mut o2 = File::create(&args.o2)?;

    let mut f1 = File::create(&args.f1)?;
    let mut f2 = File::create(&args.f2)?;

    let mut check_func: fn(entry_1: &FastqEntry, entry_2: &FastqEntry, params: &RuntimeParams) -> bool = check_read;
    if !params.rm_tiles.is_empty() {
        check_func = tile_check_read;
    }

    loop {
        let entry_1 = i1_reader.read_entry();
        let entry_2 = i2_reader.read_entry();

        match (entry_1, entry_2) {
            (Some(mut e1), Some(mut e2)) => {
                if !check_reads(&e1, &e2, &params) {
                    e1.output(&mut f1);
                    e2.output(&mut f2);
                } else {
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

    Ok(()) 
}

// cargo run -- --f1 this --f2 that --i1 other --i2 another --o1 more --o2 things

