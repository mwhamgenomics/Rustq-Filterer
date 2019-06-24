extern crate env_logger;
extern crate flate2;
extern crate log;
extern crate structopt;

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead,BufReader,Write,BufWriter};
use std::path::PathBuf;
use flate2::read::GzDecoder;
use log::{info,debug};
use structopt::StructOpt;


#[derive(StructOpt)]
struct Cli {
    #[structopt(long="i1")]
    i1: PathBuf,

    #[structopt(long="i2")]
    i2: PathBuf,

    #[structopt(long="f1")]
    f1: Option<PathBuf>,

    #[structopt(long="f2")]
    f2: Option<PathBuf>,

    #[structopt(long="o1")]
    o1: Option<PathBuf>,

    #[structopt(long="o2")]
    o2: Option<PathBuf>,

    #[structopt(long="threshold", default_value="36")]
    len_threshold: usize,

    #[structopt(long="stats_file", parse(from_os_str))]
    stats_file: Option<PathBuf>,

    #[structopt(long="remove_tiles")]
    remove_tiles: Vec<String>,

    #[structopt(long="remove_reads", parse(from_os_str))]
    remove_reads: Option<PathBuf>,

    #[structopt(long="trim_r1")]
    trim_r1: Option<i32>,

    #[structopt(long="trim_r2")]
    trim_r2: Option<i32>
}


struct FastqEntry {
    id: String,
    seq: String,
    strand: String,
    qual: String,
    tile_id: String,
    read_id: String
}


impl FastqEntry {
    fn new() -> FastqEntry {
        FastqEntry {
            id: String::new(),
            seq: String::new(),
            strand: String::new(),
            qual: String::new(),
            tile_id: String::new(),
            read_id: String::new()
        }
    }

    fn to_string(&self) -> String {
        format!(
            "read: {}seq: {}strand: {}qual: {}tile: {}",
            self.id, self.seq, self.strand, self.qual, self.tile_id
        )
    }

    fn output(&mut self, file: &mut BufWriter<File>) {
        file.write(self.id.as_bytes());
        file.write(self.seq.as_bytes());
        file.write(self.strand.as_bytes());
        file.write(self.qual.as_bytes());
    }
}


struct FastqChecker {
    reader: BufReader<GzDecoder<File>>,
}


impl FastqChecker {
    fn new(fp: &PathBuf) -> FastqChecker {
        FastqChecker {
            reader: BufReader::new(GzDecoder::new(File::open(fp).unwrap())),
        }
    }

    fn read_entry(&mut self) -> Option<FastqEntry> {
        let mut entry = FastqEntry::new();
        self.reader.read_line(&mut entry.id);
        self.reader.read_line(&mut entry.seq);
        self.reader.read_line(&mut entry.strand);
        self.reader.read_line(&mut entry.qual);

        if entry.id.is_empty() ||
           entry.seq.is_empty() ||
           entry.strand.is_empty() ||
           entry.qual.is_empty() {
            None
        } else {
            let space = &entry.id.find(" ").unwrap();
            let read_id = &entry.id[0..*space];
            let parts = &mut read_id.split(":");
            let tile_id = parts.nth(4).unwrap().to_string();

            entry.tile_id = tile_id;
            entry.read_id = read_id.to_string();

            Some(entry)
        }
    }
}


struct RuntimeInfo<'a> {
    args: &'a Cli,
    rm_tiles: HashSet<String>,
    rm_reads: HashSet<String>,
    criteria: Vec<&'a Fn(&FastqEntry, &FastqEntry, &RuntimeInfo) -> bool>,
    f1: PathBuf,
    f2: PathBuf,
    o1: PathBuf,
    o2: PathBuf,
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
                RuntimeInfo::build_rm_reads(file_path.to_path_buf(), &mut rm_reads).expect("Could not build rm_reads from file");
                criteria.push(&id_check_read);
            },
            None => {}
        }

        RuntimeInfo {
            args,
            rm_tiles,
            rm_reads,
            criteria,
            f1: RuntimeInfo::infer_output_path(&args.f1, &args.i1, "_filtered_reads.fastq"),
            f2: RuntimeInfo::infer_output_path(&args.f2, &args.i2, "_filtered_reads.fastq"),
            o1: RuntimeInfo::infer_output_path(&args.o1, &args.i1, "_filtered.fastq"),
            o2: RuntimeInfo::infer_output_path(&args.o2, &args.i2, "_filtered.fastq"),
            read_pairs_checked: 0,
            read_pairs_removed: 0,
            read_pairs_remaining: 0,
        }
    }

    fn infer_output_path(fp: &Option<PathBuf>, input_file: &PathBuf, default_file_ext: &str) -> PathBuf {
        match fp {
            Some(file_path) => file_path.to_path_buf(),
            None => {
                let input_file_slice = input_file.to_str().unwrap();
                let base;
                if input_file_slice.ends_with(".fastq.gz") {
                    base = &input_file_slice[0..input_file_slice.len()-9];
                } else {
                    base = &input_file_slice[0..input_file_slice.len()-6];
                }

                let mut output_file = base.to_string();
                output_file.push_str(default_file_ext);
                let output_file = PathBuf::from(output_file);
                output_file
            }
        }
    }

    fn build_rm_tiles(input_tiles: &Vec<String>, output_tiles: &mut HashSet<String>) {
        debug!("Removing tiles: {:?}", input_tiles);
        for t in input_tiles {
            output_tiles.insert(t.to_string());
        }
    }

    fn build_rm_reads(input_reads: PathBuf, output_reads: &mut HashSet<String>) -> std::io::Result<()> {
        debug!("Removing reads in {:?}", input_reads);
        let f = File::open(input_reads)?;
        let f = BufReader::new(f);
        for line in f.lines() {
            let read_id = format!("@{}", line.unwrap().split(" ").nth(0).unwrap());
            output_reads.insert(read_id);
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


fn write_stats_file(fp: PathBuf, info: RuntimeInfo) -> std::io::Result<()> {
    let mut report = format!(
    "r1i {:?}\nr1o {:?}\nr1f {:?}\nr2i {:?}\nr2o {:?}\nr2f {:?}\n\
        read_pairs_checked {}\nread_pairs_removed {}\nread_pairs_remaining {}\nfilter_threshold {}\n",
    info.args.i1, info.o1, info.f1, info.args.i2, info.o2, info.f2,
    info.read_pairs_checked, info.read_pairs_removed, info.read_pairs_remaining, info.args.len_threshold
    );

    if !info.rm_tiles.is_empty() {
        let mut rm_tiles = Vec::new();
        for t in &info.rm_tiles {
            rm_tiles.push(t);
        }
        rm_tiles.sort();
        report = format!("{}remove_tiles {:?}\n", report, rm_tiles);
    }

    match &info.args.remove_reads {
        Some(file_path) => {
            report = format!("{}remove_reads {:?}\n", report, file_path.to_str());
        },
        None => {}
    }

    let mut f = File::create(&fp)?;
    f.write(report.as_bytes()).expect("Could not write stats file");

    Ok(())
}


fn main() -> std::io::Result<()> {
    env_logger::init();
    info!("Starting");
    let args = Cli::from_args();

    let mut info = RuntimeInfo::new(&args);

    let mut checker_1 = FastqChecker::new(&args.i1);
    let mut checker_2 = FastqChecker::new(&args.i2);

    let mut o1 = BufWriter::new(File::create(&info.o1)?);
    let mut o2 = BufWriter::new(File::create(&info.o2)?);

    let mut f1 = BufWriter::new(File::create(&info.f1)?);
    let mut f2 = BufWriter::new(File::create(&info.f2)?);

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
                info!("Finished");
                break
            }
        }
    }

    match &info.args.stats_file {
        Some(file_path) => {
            write_stats_file(file_path.to_path_buf(), info).expect("Could not write stats file");
        },
        None => {}
    }

    Ok(()) 
}
