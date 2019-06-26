extern crate env_logger;
extern crate flate2;
extern crate log;
extern crate structopt;

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead,BufReader,Write,BufWriter,Result};
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

    fn clear(&mut self) {
        self.id.clear();
        self.seq.clear();
        self.strand.clear();
        self.qual.clear();
        self.tile_id.clear();
        self.read_id.clear();
    }

    fn to_string(&self) -> String {
        format!(
            "read: {}seq: {}strand: {}qual: {}tile: {}",
            self.id, self.seq, self.strand, self.qual, self.tile_id
        )
    }
}


struct FastqHandler {
    reader: BufReader<GzDecoder<File>>,
    mask: FastqEntry,
    output_file: BufWriter<File>,
    filtered_file: BufWriter<File>
}


impl FastqHandler {
    fn new(input_file: &PathBuf, output_file: &Option<PathBuf>, filtered_file: &Option<PathBuf>) -> FastqHandler {
        let output_file = FastqHandler::infer_output_path(output_file, input_file, "_filtered.fastq");
        let filtered_file = FastqHandler::infer_output_path(filtered_file, input_file, "_filtered_reads.fastq");

        FastqHandler {
            reader: BufReader::new(GzDecoder::new(File::open(input_file).unwrap())),
            mask: FastqEntry::new(),
            output_file: BufWriter::new(File::create(&output_file).expect("Could not open output file")),
            filtered_file: BufWriter::new(File::create(&filtered_file).expect("Could not open filtered file"))

        }
    }

    fn is_empty(&self) -> bool {
        self.mask.id.is_empty()
    }

    fn read_entry(&mut self) -> bool {
        self.mask.clear();
        self.reader.read_line(&mut self.mask.id).expect("Could not read from fastq");
        self.reader.read_line(&mut self.mask.seq).expect("Could not read from fastq");
        self.reader.read_line(&mut self.mask.strand).expect("Could not read from fastq");
        self.reader.read_line(&mut self.mask.qual).expect("Could not read from fastq");

        if !self.mask.id.is_empty() {
            let space = &self.mask.id.find(" ").unwrap();
            let read_id = &self.mask.id[0..*space];
            let parts = &mut read_id.split(":");
            let tile_id = parts.nth(4).unwrap().to_string();

            self.mask.tile_id = tile_id;
            self.mask.read_id = read_id.to_string();

            true
        } else {
            false
        }
    }

    fn output_entry(&mut self) {
        self.output_file.write(self.mask.id.as_bytes());
        self.output_file.write(self.mask.seq.as_bytes());
        self.output_file.write(self.mask.strand.as_bytes());
        self.output_file.write(self.mask.qual.as_bytes());
    }

    fn filter_entry(&mut self) {
        self.filtered_file.write(self.mask.id.as_bytes());
        self.filtered_file.write(self.mask.seq.as_bytes());
        self.filtered_file.write(self.mask.strand.as_bytes());
        self.filtered_file.write(self.mask.qual.as_bytes());
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
}


struct FastqPairChecker<'a> {
    args: &'a Cli,
    r1: FastqHandler,
    r2: FastqHandler,
    rm_tiles: HashSet<String>,
    rm_reads: HashSet<String>,
    criteria: Vec<&'a Fn(&Self) -> bool>,
    read_pairs_checked: i64,
    read_pairs_removed: i64,
    read_pairs_remaining: i64,
}


impl<'a> FastqPairChecker <'a>{
    fn new(args: &'a Cli) -> FastqPairChecker<'a> {
        let mut rm_tiles = HashSet::new();
        let mut rm_reads = HashSet::new();

        let mut criteria: Vec<& Fn(&Self) -> bool> = vec![&FastqPairChecker::check_read];

        if !args.remove_tiles.is_empty() {
            FastqPairChecker::build_rm_tiles(&args.remove_tiles, &mut rm_tiles);
            criteria.push(&FastqPairChecker::tile_check_read);
        }

        match &args.remove_reads {
            Some(file_path) => {
                FastqPairChecker::build_rm_reads(file_path.to_path_buf(), &mut rm_reads).expect("Could not build rm_reads from file");
                criteria.push(&FastqPairChecker::id_check_read);
            },
            None => {}
        }

        FastqPairChecker {
            args,
            r1: FastqHandler::new(&args.i1, &args.o1, &args.f1),
            r2: FastqHandler::new(&args.i2, &args.o2, &args.f2),
            rm_tiles,
            rm_reads,
            criteria,
            read_pairs_checked: 0,
            read_pairs_removed: 0,
            read_pairs_remaining: 0,
        }
    }

    fn build_rm_tiles(input_tiles: &Vec<String>, output_tiles: &mut HashSet<String>) {
        debug!("Removing tiles: {:?}", input_tiles);
        for t in input_tiles {
            output_tiles.insert(t.to_string());
        }
    }

    fn build_rm_reads(input_reads: PathBuf, output_reads: &mut HashSet<String>) -> Result<()> {
        debug!("Removing reads in {:?}", input_reads);
        let f = File::open(input_reads)?;
        let f = BufReader::new(f);
        for line in f.lines() {
            let read_id = format!("@{}", line.unwrap().split(" ").nth(0).unwrap());
            output_reads.insert(read_id);
        }
        Ok(())
    }

    fn check_read(&self) -> bool {
        self.r1.mask.seq.chars().count() > self.args.len_threshold && self.r2.mask.seq.chars().count() > self.args.len_threshold
    }

    fn tile_check_read(&self) -> bool {
        let tiles = &self.rm_tiles;
        if tiles.contains(&self.r1.mask.tile_id) {
            false
        } else {
            true
        }
    }

    fn id_check_read(&self) -> bool {
        let reads = &self.rm_reads;
        if reads.contains(&self.r1.mask.read_id) {
            false
        } else {
            true
        }
    }

    fn check_reads(&self) -> bool {
        let mut result: bool = true;
        for check_func in &self.criteria {
            if !check_func(self) {
                result = false
            }
        }
        result
    }

    fn write_stats_file(&self) -> Result<()> {
        match &self.args.stats_file {
            Some(file_path) => {
                let mut report = format!(
                    "r1i {:?}\nr1o {:?}\nr1f {:?}\nr2i {:?}\nr2o {:?}\nr2f {:?}\n\
                    read_pairs_checked {}\nread_pairs_removed {}\nread_pairs_remaining {}\nfilter_threshold {}\n",
                    self.args.i1, self.args.o1, self.args.f1, self.args.i2, self.args.o2, self.args.f2,
                    self.read_pairs_checked, self.read_pairs_removed, self.read_pairs_remaining, self.args.len_threshold
                );

                if !self.rm_tiles.is_empty() {
                    let mut rm_tiles = Vec::new();
                    for t in &self.rm_tiles {
                        rm_tiles.push(t);
                    }
                    rm_tiles.sort();
                    report = format!("{}remove_tiles {:?}\n", report, rm_tiles);
                }

                match &self.args.remove_reads {
                    Some(file_path) => {
                        report = format!("{}remove_reads {:?}\n", report, file_path.to_str());
                    },
                    None => {}
                }

                let mut f = File::create(&file_path)?;
                f.write(report.as_bytes()).expect("Could not write stats file");
            },
            None => {}
        }
        Ok(())
    }

    fn run(&mut self) -> Result<()> {
        info!("Starting");
        loop {
            let read_1 = self.r1.read_entry();
            let read_2 = self.r2.read_entry();

            if read_1 && read_2 {
                self.read_pairs_checked += 1;
                if !self.check_reads() {
                    self.read_pairs_removed += 1;
                    self.r1.filter_entry();
                    self.r2.filter_entry();
                } else {
                    self.read_pairs_remaining += 1;
                    self.r1.output_entry();
                    self.r2.output_entry();
                }
            } else {
                info!("Finished");
                break
            }
        }
        self.write_stats_file();
        Ok(())
    }
}


fn main() -> Result<()> {
    env_logger::init();
    let args = Cli::from_args();
    let mut info = FastqPairChecker::new(&args);
    info.run()
}
