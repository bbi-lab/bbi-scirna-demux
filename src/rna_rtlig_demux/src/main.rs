#![allow(unused_parens)]
#![allow(unused_variables)]
//#![allow(warnings)] 

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

/*
** Maximum number of microtiter plates of any type, e.g. RT, Ligation, PCR.
*/
const MAX_NUM_PLATES: usize = 16;


/*
** Input reads from BAM
**   https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md#input-reads-from-bam-files
**   https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md#input-reads-from-bam-files
**   --soloInputSAMattrBarcodeSeq CR UR    --soloInputSAMattrBarcodeQual CY UY
**   --readFilesIn input.bam --readFilesType SAM SE
**   --readFilesCommand samtools view -F 0x100
**   To output BAM tags into SAM/BAM file, add them to the list of standard tags in
**   --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM
**   Any combinations of tags can be used.
**   CR/UR: raw (uncorrected) CellBarcode/UMI
**   CY/UY: quality score for CellBarcode/UMI
**   GX/GN: for gene ID/names
**   sS/sQ: for sequence/quality combined CellBarcode and UMI; sM for barcode match status.
**   CB/UB: corrected CellBarcode/UMI. Note, that these tags require sorted BAM output, i.e. we need to add the following option:
**   --outSAMtype BAM SortedByCoordinate
**
**   https://gatk.broadinstitute.org/hc/en-us/articles/360036351132-FastqToSam-Picard
**
**   https://docs.rs/rust-htslib/latest/rust_htslib/bam/struct.Writer.html
*/

/*
** Notes:
**   o well P01-A01 has index 1
**   o index 1 encoded as base4 has sequence AAAAAAC
**   o index 0 represents no barcode; e.g., no p5 sequence.
*/


use std::io::Write;
use std::error::Error;
use std::collections::HashMap;
use std::collections::HashSet;
use base_custom::BaseCustom;
use regex::Regex;
extern crate clap;
use clap::{Arg, Command};
use bio::io::fastq::{FastqRead, Record, Writer};
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use std::io::BufReader;
use std::io::BufWriter;
use serde::Deserialize;
use serde::Serialize;
use rna_rtlig_demux::barcode_utils;
use chrono;


/*
** Define command line arguments.
*/
fn set_cl_options() -> Result<clap::Command, Box<dyn Error>> {
  let cl_options = Command::new("process_hashes")
        .version(env!("CARGO_PKG_VERSION"))
        .about("Demultiplex fastq files on RT and ligation barcodes.")
        .arg(Arg::new("fastq_read_1")   //required=true, no default
                  .required(true)
                  .short('1')
                  .long("read_1")
                  .help("Input read 1 fastq file path."))
        .arg(Arg::new("fastq_read_2")   // required=true, no default
                  .required(true)
                  .short('2')
                  .long("read_2")
                  .help("Input read 2 fastq file path."))
        .arg(Arg::new("samplesheet_file")  // required=true, no default
                  .required(true)
                  .short('s')
                  .long("samplesheet")
                  .help("Input samplesheet JSON file path."))
        .arg(Arg::new("rt_barcode_file")   //required=true, no default
                  .required(true)
                  .short('r')
                  .long("rt_file")
                  .help("Input default RT barcode sequence file path."))
        .arg(Arg::new("ligation_barcode_file")   //required=true, no default
                  .required(true)
                  .short('l')
                  .long("ligation_file")
                  .help("Input default ligation barcode sequence file path."))
        .arg(Arg::new("uncompressed_fastqs")  // required=false, default false
                  .short('u')
                  .long("uncompressed_fastqs")
                  .number_of_values(0)
                  .help("Do not compress output fastq files."))
        .arg(Arg::new("number_bam_threads")  // required=false, default 1
                  .long("ncpu")
                  .number_of_values(1)
                  .default_value("1")
                  .value_parser(clap::value_parser!(usize))
                  .help("Number of threads for BAM compression."))
        .arg(Arg::new("output_file_format")     // required=false, default bam
                  .required(false)
                  .short('f')
                  .long("output_format")
                  .value_parser(["bam", "fastq"])
                  .default_value("bam")
                  .help("Output file format: fastq or bam."));
  Ok(cl_options)
}


/*
** Get the lane, P7, and P5 indices from the input
** fastq filename.
*/
fn get_file_indices(filename: &str) -> Result<HashMap<String, usize>, Box<dyn Error>> {
  let path = std::path::Path::new(filename);
  let filename = path.file_name().unwrap();
  let binding = filename.to_string_lossy();
  let parts: Vec<_> = binding.split("_").collect();

  let lane_index: usize = parts[0].parse().unwrap();
  let p7_index: usize = parts[1].parse().unwrap();
  let p5_index: usize = parts[2].parse().unwrap();

  let mut file_indices: HashMap<String, usize> = HashMap::new();
  file_indices.insert("lane".to_string(), lane_index);
  file_indices.insert("p7".to_string(), p7_index);
  file_indices.insert("p5".to_string(), p5_index);

  Ok(file_indices)
}


/*
** Make a HashMap that maps utiter plate wells to
** well indices. Well A01 has index 1.
*/
fn make_well_index_map(max_index: usize, across_row_first: bool, with_plate: bool) -> Result<HashMap<String, usize>, Box<dyn Error>> {
  let mut well_index_map: HashMap<String, usize> = HashMap::with_capacity(max_index);
  let mut ipl: usize;
  let mut well: String;
  for i in 1..(max_index+1) {
    (ipl, well) = barcode_utils::index_to_well(i, across_row_first).unwrap();
    if(with_plate) {
      let plate_well: String = format!("P{:02}-{}", ipl, well);
      let _ = well_index_map.insert(plate_well, i);
    }
    else {
      let _ = well_index_map.insert(well, i);
    }
  }
  Ok(well_index_map)
}


/*
** Make a HashMap that maps utiter plate well indices
** to well names. Well A01 has index 1.
*/
/*
fn make_index_well_map(max_index: usize, across_row_first: bool, with_plate: bool) -> Result<HashMap<usize, String>, Box<dyn Error>> {
  let mut index_well_map: HashMap<usize, String> = HashMap::with_capacity(max_index);
  let mut ipl: usize;
  let mut well: String;
  for i in 1..(max_index+1) {
    (ipl, well) = barcode_utils::index_to_well(i, across_row_first).unwrap();
    if(with_plate) {
      let plate_well: String = format!("P{:02}-{}", ipl, well);
      let _ = index_well_map.insert(i, plate_well);
    }
    else {
      let _ = index_well_map.insert(i, well);
    }
  }
  Ok(index_well_map)
}
*/

/// Convert base10 barcode index to a base4 index encoded
/// as a String with 'A'=0, 'C'=1, 'G'=2, and 'T'=3.
///
/// Arguments:
///
/// max_index: maximum index value in returned vector of strings.
///
/// Return:
///
/// Vector of Strings that encode base4 indices. Vector index 0
/// is 'A'.
///
/// Note:
///
/// Convert base4 encoded index to base10 int using
///  let es = String; // base4 encoded index
///  let mut digmap: HashMap<char, u32> = HashMap::new();
///  let _ = digmap.insert('A', 0);
///  let _ = digmap.insert('C', 1);
///  let _ = digmap.insert('G', 2);
///  let _ = digmap.insert('T', 3);
///  y = es.char().rev().enumerate().map(|(ii,k)| 4u32.pow((ii) as u32) as u32 * digmap[&k]).sum();
///
fn make_index_encoder(max_index: u64) -> Result<Vec<String>, Box<dyn Error>> {
  let base4 = BaseCustom::<char>::new("ACGT".chars().collect());

  let mut convert_vec = Vec::<String>::new();
  for iv in (0..max_index+1) {
    convert_vec.push(base4.gen(iv));
  }
  Ok(convert_vec)
}


/*
** Read the samplesheet json file.
*/
fn get_samplesheet_json(filename: String) -> Result<serde_json::Value, Box<dyn Error>>  {
  let reader = std::fs::File::open(filename).unwrap();
  let v: serde_json::Value = serde_json::from_reader(reader).unwrap();
  Ok(v)
}


/*
** Make a sorted vector of distinct indices as usize values from
** a String slice of index 'ranges' from the samplesheet.json file. 
** For example, &str = "1-4,5,7".
*/
fn make_index_vec(index_str: &str) -> Result<Vec<usize>, Box<dyn Error>> {
  let re = Regex::new(r"([0-9]+)([-]([0-9]+))?").unwrap();
  let mut index_vec: Vec<usize> = Vec::new();
  for index_range in index_str.split(",") {
    let re_mats = re.captures(index_range).unwrap();
    let index1: usize = re_mats.get(1).unwrap().as_str().parse().unwrap();
    let mut index2: usize = index1;
    if(re_mats.get(3) != None) {
      index2 = re_mats.get(3).unwrap().as_str().parse().unwrap();
    }
    for idx in index1..(index2+1) {
      index_vec.push(idx);
    }
  }

  /*
  ** Find distinct indices.
  */
  let mut index_set: HashSet<usize> = HashSet::new();
  for index in index_vec.iter() {
    index_set.insert(*index);
  }

  let mut index_vec_distinct: Vec<usize> = Vec::new();
  for index in index_set {
    index_vec_distinct.push(index);
  }

  /*
  ** Sort indices.
  */
  index_vec_distinct.sort();

  Ok(index_vec_distinct)
}


/*
** Make a vector of 'sample range' structs for the samples
** in the specified lane. Separate the barcode ranges and
** expand them into sorted lists of individual indices, not
** ranges.
**
** sample ranges hashmap
**   {
**       "sample_id": "sample.2",
**       "ranges": "32-35,100-104:1-12:17-24",
**       "lanes": "5-8",
**       "tissue": "lung",
**       "genome": "mouse",
**       "hash_file": "",
**       "sample_flags": "S",
**       "external_sample_name": "extsample.2",
**       "wrap_group": "Smyth",
**       "rt_file": "",
**       "ligation_file": "",
**       "p7_file": "/net/bbi/barcodes2/p7_file.txt",
**       "p5_file": "/net/bbi/barcodes2/p5_file.txt",
**       "library": "RNA3-088"
**      }
*/

#[allow(dead_code)]
#[derive(Deserialize, Debug, Clone)]
struct SampleMap {
  sample_id: String,
  ranges: String,
  lanes: String,
  tissue: String,
  genome: String,
  hash_file: String,
  sample_flags: String,
  external_sample_name: String,
  wrap_group: String,
  rt_file: String,
  ligation_file: String,
  p7_file: String,
  p5_file: String,
  library: String,
  process_group: String
}


/*
** Deserialize a SampleMap.
**
** See
**   https://docs.rs/serde_json/latest/serde_json/fn.from_value.html 
** and
**   https://docs.rs/serde_json/latest/serde_json/value/enum.Value.html
** 
*/
fn deserialize_sample_map(serialized_sample_map: serde_json::Value) -> Result<SampleMap, Box<dyn Error>> {

  let sample_map: SampleMap = serde_json::from_value(serialized_sample_map).expect("Error: deserialize_sample_map: unable to deserialize the sample map.\n  Perhaps the samplesheet sample_index_list map changed.\n  If so, update the SampleMapStrings struct in this program.");

  Ok(sample_map)
}


/*
** Deserialize an array of SampleMaps.
*/
fn deserialize_sample_map_vector(samplesheet_json: serde_json::Value) -> Result<Vec<SampleMap>, Box<dyn Error>> {
  let mut sample_map_vec: Vec<SampleMap> = Vec::new();

  /*
  ** Loop over each sample dictionary in samplesheet_json["sample_index_list"].
  */
  for sample_index_map in samplesheet_json["sample_index_list"].as_array().unwrap() {
    sample_map_vec.push(deserialize_sample_map(sample_index_map.clone()).unwrap());
  }

  Ok(sample_map_vec)
}


/*
** Make a tuple of lane, RT, P7, and P5 index vectors given
** a sample map that contains a ranges string.
*/
fn get_sample_index_vecs(sample_map: &SampleMap) -> Result<(Vec<usize>, Vec<usize>, Vec<usize>, Vec<usize>), Box<dyn Error>> {
  let lane_index_vec: Vec<usize> = Vec::new();
  let rt_index_vec: Vec<usize> = Vec::new();
  let p7_index_vec: Vec<usize> = Vec::new();
  let p5_index_vec: Vec<usize> = Vec::new();

  let lane_index_vec = make_index_vec(&sample_map.lanes).unwrap();

  let ranges_parts: Vec<&str> = sample_map.ranges.split(":").collect();
  let rt_index_vec = make_index_vec(ranges_parts[0]).unwrap();
  let p7_index_vec = make_index_vec(ranges_parts[1]).unwrap();
  let p5_index_vec = make_index_vec(ranges_parts[2]).unwrap();

  Ok((lane_index_vec, rt_index_vec, p7_index_vec, p5_index_vec))
}


/*
** Make a vector of SampleMaps specific to the fastq file, i.e., the
** same lane and p7 and p5 barcodes.
*/
fn get_lane_samples(sample_map_vec_all: &Vec<SampleMap>, lane_index: usize, p7_index: usize, p5_index: usize) -> Result<Vec<SampleMap>, Box<dyn Error>> {

  let mut sample_map_vec: Vec<SampleMap> = Vec::new();

  /*
  ** Loop over each sample dictionary in samplesheet_json["sample_index_list"].
  */
  for sample_map in sample_map_vec_all.iter() {
    let (lane_index_vec, rt_index_vec, p7_index_vec, p5_index_vec) = get_sample_index_vecs(sample_map).unwrap();

    /*
    ** Skip if this sample entry does not include 'lane_index'.
    */
    if(!lane_index_vec.iter().any(|i| *i == lane_index)) {
      println!("skip sample/not in lanes:\n{:#?}", sample_map);
      continue;
    }

    /*    
    ** Skip if this sample entry does not include 'p7_index'.
    */    
    if(!p7_index_vec.clone().iter().any(|i| *i == p7_index)) {
      println!("skip sample/not in p7:\n{:#?}", sample_map);
      continue;
    }

    /*    
    ** Skip if this sample entry does not include 'p5_index'.
    */    
    if(!p5_index_vec.clone().iter().any(|i| *i == p5_index)) {
      println!("skip sample/not in p5:\n{:#?}", sample_map);
      continue;
    }

    sample_map_vec.push((*sample_map).clone());
  }

  Ok(sample_map_vec)
}


/*
** Barcode correction and identifier structure.
*/
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
struct BarcodeIdentifier {
  sample_index: usize,
  well_index: usize,
  well_name: String,
  corrected_flag: bool,
  num_mismatch: u32,
  read_count: u64,
}


/*
** I get a fastq file with reads from one lane. I need to find
** the fastq file lane and the samplesheet data for that lane.
** The required samplesheet values for the lane are
**   o  rt indices
**   o  rt barcode file
**   o  ligation barcode file
**   o  sample_id
**   o  get rt barcode indices for reads in lane
**
** Notes:
**   o  find rt and ligation barcode file paths
**   o  read rt and ligation barcode files as HashMap with sequence as key
**   o  make rt and ligation corrected barcode hashmaps
**   o  make well name to index hashmap with plate identifiers
**   o  make the barcode_id_map
*/
fn make_barcode_id_map(sample_map_vec: &Vec<SampleMap>, barcode_type: &str, default_file: &str, recipe: &str) -> Result<HashMap<Vec<u8>, BarcodeIdentifier>, Box<dyn Error>> {
  /*
  ** Find barcode file path.
  ** Check the samplesheet json file for a file path. If zero length, use default path.
  */
  let mut file_name = "".to_string();
  let test_name = String::new();

  /*
  **  Set barcode file path using either the value in
  **  the samplesheet or the default.
  */
  if(barcode_type == "rt_file") {
    for sample_map in sample_map_vec.iter() {
      if(sample_map.rt_file.len() > 0) {
        file_name = sample_map.rt_file.clone();
        break
      }
    }
  }
  else if(barcode_type == "ligation_file") {
    for sample_map in sample_map_vec.iter() {
      if(sample_map.rt_file.len() > 0) {
        file_name = sample_map.ligation_file.clone();
        break
      }
    }
  }
  else {
    Err("Error: make_barcode_id_map: unrecognized barcode_type value.")?;
  }

  println!("file len: {}", file_name.len());
  if(file_name.len() == 0) {
    file_name = default_file.to_string();
  }
  println!("file_name: {}", file_name);


  /*
  ** Read the barcodes into a HashMap keyed by barcode sequence.
  */
  let barcode_whitelist = rna_rtlig_demux::barcode_utils::read_barcode_file(&file_name).unwrap();
  let num_barcode: usize = barcode_whitelist.keys().len();
  if(num_barcode > MAX_NUM_PLATES * 96) {
    eprintln!("The number of {} barcodes, {}, exceeds the current available storage", barcode_type, num_barcode);
    eprintln!("space. Increase the storage space by increasing the value of MAX_NUM_PLATES");
    eprintln!("in the rna_rtlig_demux source file main.rs.");
    std::process::exit(-1);
  }
  let plate_well_index_map_by_row = make_well_index_map(MAX_NUM_PLATES*96, true, true).unwrap();

  let max_num_mismatch: usize = 1_usize;
  let barcode_correct = barcode_utils::construct_mismatch_to_whitelist_map(barcode_whitelist.keys().map(|s| s.to_string()).collect(), max_num_mismatch, true).unwrap();

  let mut barcode_id_map: HashMap<Vec<u8>, BarcodeIdentifier> = HashMap::new();

  for num_mismatch in 0..(max_num_mismatch+1) {
    println!("num mismatch: {}", num_mismatch);
    let mut corrected_flag: bool = false;
    if(num_mismatch > 0) {
      corrected_flag = true;
    }

    if(barcode_type == "rt_file") {
      for barcode_seq in barcode_correct[num_mismatch].keys() {
        let barcode_corrected = barcode_correct[num_mismatch][barcode_seq].clone();
        let sample_index: usize = 0;
        let barcode_well = barcode_whitelist[&barcode_corrected].clone();
        let barcode_index = plate_well_index_map_by_row[&barcode_well];
        let read_count: u64 = 0;
        let barcode_identifier = BarcodeIdentifier { sample_index: sample_index,
                                                     well_index: barcode_index,
                                                     well_name: barcode_well,
                                                     corrected_flag: corrected_flag,
                                                     num_mismatch: num_mismatch as u32,
                                                     read_count: read_count};
        let _ = barcode_id_map.insert(barcode_seq.as_bytes().to_vec(), barcode_identifier);
      }
    }
    else if(barcode_type == "ligation_file") {
      for barcode_seq in barcode_correct[num_mismatch].keys() {
        let barcode_corrected = barcode_correct[num_mismatch][barcode_seq].clone();
        let sample_index: usize = 0;
        let barcode_well = barcode_whitelist[&barcode_corrected].clone();

        /*
        ** The sequences in the ligation barcode file have names
        ** in the form "LIG<index>" where index starts with 1. Trim
        ** off "LIG" in the name and convert to a usize integer.
        */
        let barcode_index = barcode_well.clone()[3..].parse::<usize>().unwrap();
        let read_count: u64 = 0;
        let barcode_identifier = BarcodeIdentifier { sample_index: sample_index,
                                                     well_index: barcode_index,
                                                     well_name: barcode_well,
                                                     corrected_flag: corrected_flag,
                                                     num_mismatch: num_mismatch as u32,
                                                     read_count: read_count};
        let _ = barcode_id_map.insert(barcode_seq.as_bytes().to_vec(), barcode_identifier);
      }
    }
    else {
      Err("Error: make_barcode_id_map: unrecognized barcode_type value.")?;
    }
  }

  /*
  for key in barcode_id_map.keys() {
    println!("barcode id map key: {:#?}  val: {:#?}", key, barcode_id_map[key].well_name);
  }
  */

  Ok(barcode_id_map)
}


#[derive(Debug, Clone)]
struct SampleIndices {
  sample_name: String, // key to HashMap
  sample_index: usize,
  rt_indices: Vec<usize>
}

/*
** Notes:
**   o  the first RT well index is 1 in well_index_to_sample_index_map.
**   o  the first sample index is 1 in sample_index_to_name_map. The
**      0 index is reserved for 'undetermined'.
*/
fn make_index_maps(sample_map_vec: &Vec<SampleMap>,
                   barcode_id_map: &mut HashMap<Vec<u8>, BarcodeIdentifier>) -> Result<(Vec<usize>, Vec<String>), Box<dyn Error>> {
  /*
  ** Make a HashMap of SampleIndices structs.
  */
  let mut sample_indices_map: HashMap<String, SampleIndices> = HashMap::new();
  let mut num_sample_indices = 1;
  for sample_map in sample_map_vec.iter() {
    let (lane_index_vec, rt_index_vec, p7_index_vec, p5_index_vec) = get_sample_index_vecs(sample_map).unwrap();
    let sample_id = sample_map.sample_id.clone();
    if(!sample_indices_map.contains_key(&sample_id)) {
      let rt_index_vec: Vec<usize> = Vec::new();
      let sample_indices = SampleIndices { sample_name: sample_id.clone(), sample_index: num_sample_indices, rt_indices: rt_index_vec};
      num_sample_indices += 1;
      let _ = sample_indices_map.insert(sample_id.clone(), sample_indices);
    }

    /*
    ** Add rt indices.
    */
    let sample_indices = sample_indices_map.get_mut(&sample_id).unwrap();
    for rt_index in rt_index_vec.into_iter() {
      sample_indices.rt_indices.push(rt_index);
    }
  }

  /*
  ** Find maximum RT well index value.
  */
  let mut max_well_index: usize = 0;
  for key in barcode_id_map.keys() {
    if(barcode_id_map[key].well_index > max_well_index) {
      max_well_index = barcode_id_map[key].well_index;
    }
  }

  /*
  ** Find number of distinct sample names.
  */
  let num_samples: usize = sample_indices_map.keys().len();
  println!("num_samples: {}", num_samples);

  /*
  ** Declare and initialize well index to sample index map.
  */ 
  let mut well_index_to_sample_index_map: Vec<usize> = Vec::with_capacity(max_well_index + 1);
  for i in 0..(max_well_index + 2) {
    well_index_to_sample_index_map.push(0_usize);
  }

  /*
  ** Declare and initialize sample_index_to_name_map.
  */
  let mut sample_index_to_name_map: Vec<String> = Vec::with_capacity(num_samples);
  sample_index_to_name_map.push("Undetermined".to_string()); 
  for i in 1..(num_samples + 1) {
    sample_index_to_name_map.push(String::from(""));
  }
  
  /*
  ** Fill in  well_index_to_sample_index_map and
  ** sample_index_to_name_map values from
  ** sample_index_map.
  */
  for key in sample_indices_map.keys() {
    let sample_index = sample_indices_map[key].sample_index;
    sample_index_to_name_map[sample_index] = sample_indices_map[key].sample_name.clone();
    for index in &sample_indices_map[key].rt_indices {
      well_index_to_sample_index_map[*index] = sample_index;
    }
  }

  Ok((well_index_to_sample_index_map, sample_index_to_name_map))
}


/*
** Make a file reader for either an uncompressed or
** compressed file. Use file extension to identify
** gzipped files.
*/
fn open_reader_file(filename: &str) -> Result<Box<dyn std::io::Read>, Box<dyn Error>> {
  let file_extension = std::path::Path::new(filename)
                                .extension()
                                .and_then(std::ffi::OsStr::to_str).unwrap();
  let reader = std::fs::File::open(filename).expect(&format!("Error: unable to open file {}", &filename));
  if(file_extension == "gz") {
     let reader_gzdecoder = MultiGzDecoder::new(reader);
     Ok(Box::new(reader_gzdecoder))
  }
  else {
    Ok(Box::new(reader))
  }
}


/*
** Make a file writer for either an uncompressed or
** compressed file. Use file extension to identify
** gzipped files.
*/
fn open_writer_file(filename: &str) -> Result<Box<dyn std::io::Write>, Box<dyn Error>> {
  let file_extension = std::path::Path::new(filename)
                                .extension()
                                .and_then(std::ffi::OsStr::to_str).unwrap();
  let writer = std::fs::File::create(filename).unwrap();
  if(file_extension == "gz") {
     let writer_gzencoder = GzEncoder::new(writer, flate2::Compression::new(6));
     Ok(Box::new(writer_gzencoder))
  }
  else {
    Ok(Box::new(writer))
  }
}


/*
** Open writers to the output fastq files. Store the writers in
** a vector. The writers may compress the fastq files.
**
** Notes:
**   o  sample_index_to_name_map:
**        o  index 0 is undetermined reads
**        o  index 1 is the first sample name
*/
fn open_fastq_writers(sample_index_to_name_map: &Vec<String>,
                      well_index_to_sample_index_map: &Vec<usize>,
                      barcode_id_map: &mut HashMap<Vec<u8>, BarcodeIdentifier>,
                      file_indices: &HashMap<String, usize>,
                      uncompress_flag: bool) -> Result<Vec<Vec<Writer<Box<dyn Write>>>> , Box<dyn Error>> {
  let mut fastq_out_vec: Vec<Vec<Writer<Box<dyn std::io::Write>>>> = Vec::new();
  fastq_out_vec.push(Vec::new());  // Read 1
  fastq_out_vec.push(Vec::new());  // Read 2

  if(uncompress_flag == true) {
    for i in 0..(sample_index_to_name_map.len()) {
      let filename: String = format!("{}-{:03}_{:03}_{:03}-L{:03}-R1.fastq", sample_index_to_name_map[i], &file_indices["lane"], &file_indices["p7"], &file_indices["p5"], &file_indices["lane"]);
      println!("index r1: {}  filename: {}", i, filename);
      let writer = open_writer_file(&filename).unwrap();
      let bufwriter = BufWriter::new(writer);
      fastq_out_vec[0].push(bio::io::fastq::Writer::from_bufwriter(bufwriter));

      let filename: String = format!("{}-{:03}_{:03}_{:03}-L{:03}-R2.fastq", sample_index_to_name_map[i], &file_indices["lane"], &file_indices["p7"], &file_indices["p5"], &file_indices["lane"]);
      println!("index r2: {}  filename: {}", i, filename);
      let writer = open_writer_file(&filename).unwrap();
      let bufwriter = BufWriter::new(writer);
      fastq_out_vec[1].push(bio::io::fastq::Writer::from_bufwriter(bufwriter));
    }
  }
  else {
    for i in 0..(sample_index_to_name_map.len()) {
      let filename: String = format!("{}-{:03}_{:03}_{:03}-L{:03}-R1.fastq.gz", sample_index_to_name_map[i], &file_indices["lane"], &file_indices["p7"], &file_indices["p5"], &file_indices["lane"]);
      println!("index r1: {}  filename: {}", i, filename);
      let writer = open_writer_file(&filename).unwrap();
      let bufwriter = BufWriter::new(writer);
      fastq_out_vec[0].push(bio::io::fastq::Writer::from_bufwriter(bufwriter));

      let filename: String = format!("{}-{:03}_{:03}_{:03}-L{:03}-R2.fastq.gz", sample_index_to_name_map[i], &file_indices["lane"], &file_indices["p7"], &file_indices["p5"], &file_indices["lane"]);
      let writer = open_writer_file(&filename).unwrap();
      let bufwriter = BufWriter::new(writer);
      fastq_out_vec[1].push(bio::io::fastq::Writer::from_bufwriter(bufwriter));
    }
  }

  /*
  ** Set sample_index in barcode_id_map entries.
  */
  let mut key_vec: Vec<Vec<u8>> = Vec::new();
  for key in barcode_id_map.keys() {
    key_vec.push(key.clone());
  }
  for key in key_vec.into_iter() {
    let well_index = barcode_id_map.entry(key.clone()).or_default().well_index;
    barcode_id_map.entry(key.clone()).or_default().sample_index = well_index_to_sample_index_map[well_index];
  }

  Ok(fastq_out_vec)
}


/*
** Open writers to the output bam files. Store the writers in
** a vector.
**
** Notes:
**   o  sample_index_to_name_map:
**        o  index 0 is undetermined reads
**        o  index 1 is the first sample name
*/
fn open_bam_writers(sample_index_to_name_map: &Vec<String>,
                    well_index_to_sample_index_map: &Vec<usize>,
                    barcode_id_map: &mut HashMap<Vec<u8>, BarcodeIdentifier>,
                    file_indices: &HashMap<String, usize>,
                    num_threads: usize) -> Result<Vec<Box<rust_htslib::bam::Writer>> , Box<dyn Error>> {
  let mut bam_out_vec: Vec<Box<rust_htslib::bam::Writer>> = Vec::new();
  let date = chrono::offset::Local::now();

  for i in 0..(sample_index_to_name_map.len()) {
      let filename: String = format!("{}-{:03}_{:03}_{:03}-L{:03}.bam", sample_index_to_name_map[i], &file_indices["lane"], &file_indices["p7"], &file_indices["p5"], &file_indices["lane"]);
      let mut header = rust_htslib::bam::Header::new();
      header.push_comment(b"Made by rt_lig_demux");
      header.push_comment(format!("rt_lig_demux run date: {}", date.to_string()).as_bytes());
      bam_out_vec.push(Box::new(rust_htslib::bam::Writer::from_path(filename, &header, rust_htslib::bam::Format::Bam).unwrap()));
  }

  if(num_threads > 1) {
    for i in 0..(sample_index_to_name_map.len()) {
      let _ = bam_out_vec[i].set_threads(4_usize);
    }
  }

  /* 
  ** Set sample_index in barcode_id_map entries.
  */ 
  let mut key_vec: Vec<Vec<u8>> = Vec::new();
  for key in barcode_id_map.keys() {
    key_vec.push(key.clone());
  }  
  for key in key_vec.into_iter() {
    let well_index = barcode_id_map.entry(key.clone()).or_default().well_index;
    barcode_id_map.entry(key.clone()).or_default().sample_index = well_index_to_sample_index_map[well_index];
  } 

  Ok(bam_out_vec)
}


/*
** process_reads needs:
**   o  two readers (done)
**   o  many writers
**   o  barcode decoder HashMaps: rt and ligation
**   o  PCR p7 and p5 well names
**
** fastq read header looks like
**   @SeahubZ01.fq.part1-P5F02-P7E03_1|SeahubZ01.fq.part1|F02|E03|P05-G10_LIG76|AAAATCAA
*/
fn process_reads(fastq1_file: &str,
                 fastq2_file: &str,
                 rt_barcode_id_map: &mut HashMap<Vec<u8>, BarcodeIdentifier>,
                 ligation_barcode_id_map: &mut HashMap<Vec<u8>, BarcodeIdentifier>,
                 file_indices: &HashMap<String, usize>,
                 recipe: &str,
                 index_encoder: &Vec<String>,
                 output_file_format: &str,
                 uncompress_flag: bool,
                 num_threads: usize,
                 well_index_to_sample_index_map: &Vec<usize>,
                 sample_index_to_name_map: &Vec<String>,
                 log_read_counts: &mut HashMap<String, u64>) -> Result<String , Box<dyn Error>> {

  let log = String::new();

  let (ipl, p7_well_name) = barcode_utils::index_to_well(file_indices["p7"], true).unwrap();
  let (ipl, p5_well_name) = barcode_utils::index_to_well(file_indices["p5"], false).unwrap();

  let p7_index_encoded = index_encoder[file_indices["p7"]].clone();
  let p5_index_encoded = index_encoder[file_indices["p5"]].clone();

/*
  println!("p7 well: {}", p7_well_name);
  println!("p5 well: {}", p5_well_name);

  println!("p7 index ({}) encoded: {:A>7}", file_indices["p7"], p7_index_encoded);
  println!("p5 index ({}) encoded: {:A>7}", file_indices["p5"], p5_index_encoded);
*/

  let mut fastq_record1 = Record::new();
  let mut fastq_record2 = Record::new();

  let mut reader1 = bio::io::fastq::Reader::new(BufReader::new(open_reader_file(fastq1_file).unwrap()));
  let mut reader2 = bio::io::fastq::Reader::new(BufReader::new(open_reader_file(fastq2_file).unwrap()));

  /*
  ** Get output file writers.
  */
  let mut fastq_writer_vec: Vec<Vec<Writer<Box<dyn std::io::Write>>>> = Vec::new();
  let mut bam_writer_vec: Vec<Box<rust_htslib::bam::Writer>> = Vec::new();
  if(output_file_format == "fastq") {
    fastq_writer_vec = open_fastq_writers(&sample_index_to_name_map, &well_index_to_sample_index_map, rt_barcode_id_map, file_indices, uncompress_flag).unwrap();
  }
  else
  if(output_file_format == "bam") {
    bam_writer_vec = open_bam_writers(&sample_index_to_name_map,
                                      &well_index_to_sample_index_map,
                                      rt_barcode_id_map,
                                      &file_indices,
                                      num_threads).unwrap();
  }

  let mut nrecord = 0_u64;
  let mut count_rt_9 = 0_u64;
  let mut count_lg_9 = 0_u64;
  let mut count_rt_10 = 0_u64;
  let mut count_lg_10 = 0_u64;
  let mut count_rtlg_9 = 0_u64;
  let mut count_rtlg_10 = 0_u64;
  let mut count_rtlg_9_and_rtlg_10 = 0_u64;
  let mut count_matches_failed = 0_u64;

  let mut lig_9_slice: std::ops::Range<usize> = 0..0;
  let mut umi_9_slice: std::ops::Range<usize> = 0..0;
  let mut rt_9_slice: std::ops::Range<usize> = 0..0;
  let mut lig_10_slice: std::ops::Range<usize> = 0..0;
  let mut umi_10_slice: std::ops::Range<usize> = 0..0;
  let mut rt_10_slice: std::ops::Range<usize> = 0..0;

  let mut qual_bam: Vec<u8> = Vec::with_capacity(1024);

  if(recipe == "std_1") {
    lig_9_slice = 0..9;
    umi_9_slice = 15..23;
    rt_9_slice = 23..33;
    lig_10_slice = 0..10;
    umi_10_slice = 16..24;
    rt_10_slice = 24..34;
  }
  else {
    Err("Error: unrecognized recipe value.")?;
  }

  let mut records_written_index = 0;

  /*
  ** Process the reads from the input fastq file pair.
  */
  loop {
    nrecord += 1;

/*
    if(nrecord >= 10) {
      break;
    }
*/

    reader1.read(&mut fastq_record1).expect(&format!("oh-no 1: nrecord: {}", nrecord));
    reader2.read(&mut fastq_record2).expect(&format!("oh-no 2: nrecord: {}", nrecord));
    if(fastq_record1.is_empty() || fastq_record2.is_empty()) {
      break;
    }

    /*
    ** Check that names of reads 1 and 2 match.
    */
    if(fastq_record1.id() != fastq_record2.id()) {
      eprintln!("Error: read 1 and read 2 ids differ: id 1: {} id 2: {}",
                fastq_record1.id(),
                fastq_record2.id());
      std::process::exit(-1);
    }

/*
    let seq1 = fastq_record1.seq().clone();
    let seq2 = fastq_record2.seq().clone();
*/
    let seq1 = fastq_record1.seq();
    let seq2 = fastq_record2.seq();

    /*
    ** 9 base ligation barcode candidate.
    */
    let lig_read_9  = &seq1[lig_9_slice.clone()];
    let umi_read_9  = &seq1[umi_9_slice.clone()];
    let rt_read_9   = &seq1[rt_9_slice.clone()];

    /*
    ** 10 base ligation barcode candidate.
    */
    let lig_read_10 = &seq1[lig_10_slice.clone()];
    let umi_read_10 = &seq1[umi_10_slice.clone()];
    let rt_read_10  = &seq1[rt_10_slice.clone()];

    /*
    ** Check for matches to known barcodes with and
    ** without mismatches.
    */
    #[warn(unused_assignments)]
    let mut sample_index: usize       = 0_usize;
    #[warn(unused_assignments)]
    let mut rt_well_index: usize      = 0_usize;
    #[warn(unused_assignments)]
    let mut rt_well_name: String      = String::new();
    #[warn(unused_assignments)]
    let mut lig_well_name: String     = String::new();
    #[warn(unused_assignments)]
    let mut rt_index_encoded: String  = String::new();
    #[warn(unused_assignments)]
    let mut lig_index_encoded: String = String::new();
    #[warn(unused_assignments)]
    let mut umi_seq: &[u8];

    /*
    ** Interior blocks limit scope of mutable borrows.
    */
    #[warn(unused_assignments)]
    let mut rtlg_9_and_rtlg_10_flag = false;
    #[warn(unused_assignments)]
    let mut rt_match_9_flag: bool   = false;
    #[warn(unused_assignments)]
    let mut lig_match_9_flag: bool  = false;
    #[warn(unused_assignments)]
    let mut rtlg_9_flag: bool       = false;
    {
      let mut rt_match_9_bi: &mut BarcodeIdentifier = &mut Default::default();
      let rt_match_9_option  = rt_barcode_id_map.get_mut(rt_read_9);
      match rt_match_9_option {
        Some(rt_match_bi) => { rt_match_9_bi    = rt_match_bi;
                               rt_match_9_flag  = true;
                               count_rt_9      += 1;
                               rt_match_9_bi.read_count += 1 },
        None => rt_match_9_flag = false,
      }

      let mut lig_match_9_bi: &mut BarcodeIdentifier = &mut Default::default();
      let lig_match_9_option  = ligation_barcode_id_map.get_mut(lig_read_9);
      match lig_match_9_option {
        Some(lig_match_bi) => { lig_match_9_bi = lig_match_bi;
                                lig_match_9_flag = true;
                                count_lg_9 += 1;
                                lig_match_9_bi.read_count += 1 },
        None => lig_match_9_flag = false,
      }

      if(rt_match_9_flag && lig_match_9_flag) {
        umi_seq           = umi_read_9;
        sample_index      = rt_match_9_bi.sample_index;
        rt_well_index     = rt_match_9_bi.well_index;
        rt_well_name      = rt_match_9_bi.well_name.clone();
        lig_well_name     = lig_match_9_bi.well_name.clone();
        rt_index_encoded  = index_encoder[rt_match_9_bi.well_index].clone();
        lig_index_encoded = index_encoder[lig_match_9_bi.well_index].clone();
        rtlg_9_flag       = true;
        count_rtlg_9     += 1;
      }
      else {
        umi_seq = b"";
      }
    }
    

    let mut rt_match_10_flag: bool  = false;
    let mut lig_match_10_flag: bool = false;
    let mut rtlg_10_flag: bool      = false;
    {
      let mut rt_match_10_bi: &mut BarcodeIdentifier = &mut Default::default();
      let rt_match_10_option = rt_barcode_id_map.get_mut(rt_read_10);
      match rt_match_10_option {
        Some(rt_match_bi) => { rt_match_10_bi = rt_match_bi;
                               rt_match_10_flag = true;
                               count_rt_10 += 1;
                               rt_match_10_bi.read_count += 1 },
        None => rt_match_10_flag = false,
      }

      let mut lig_match_10_bi: &mut BarcodeIdentifier = &mut Default::default();
      let lig_match_10_option = ligation_barcode_id_map.get_mut(lig_read_10);
      match lig_match_10_option {
        Some(lig_match_bi) => { lig_match_10_bi = lig_match_bi;
                                lig_match_10_flag = true;
                                count_lg_10 += 1;
                                lig_match_10_bi.read_count += 1 },
        None => lig_match_10_flag = false,
      }

      if(rt_match_10_flag && lig_match_10_flag) {
        if(rtlg_9_flag == false) {
          umi_seq           = umi_read_10;
          sample_index      = rt_match_10_bi.sample_index;
          rt_well_index     = rt_match_10_bi.well_index;
          rt_well_name      = rt_match_10_bi.well_name.clone();
          lig_well_name     = lig_match_10_bi.well_name.clone();
          rt_index_encoded  = index_encoder[rt_match_10_bi.well_index].clone();
          lig_index_encoded = index_encoder[lig_match_10_bi.well_index].clone();
          rtlg_10_flag      = true;
          count_rtlg_10    += 1;
        }
        else {
          rtlg_9_and_rtlg_10_flag = true; 
          umi_seq = b"";
        }
      }
    }

    /*
    ** Reject if both rtlg9_flag and rtlg10_flag are true.
    */
    if(rtlg_9_and_rtlg_10_flag) {
      count_rtlg_9_and_rtlg_10 += 1;
      continue;
    }

    if(rtlg_9_flag == false && rtlg_10_flag == false)
    {
      count_matches_failed += 1;
      continue;
    }

    /*
    ** Write the reads.
    **   o  header w/p5 @24.0133-P5E02-P7B04_3|24.0133|E02|B04|P06-B10_LIG374|ATTGCTAA
    **   o  header wo/p5 @24.0053-P5none-P7G08_2|24.0053|none|G08|P07-A03_LIG43|TCTGTTGT
    **
    ** What do I need make before writing reads?
    **   o  read1 sequence
    **        o  rt/lig/p7/p5 index encodings
    **   o  read header string
    **        o  sample name
    **        o  p7 and p5 well names
    **        o  rt well name
    **        o  ligation well name
    */

    /*
    ** Limit scope of mutable borrows.
    */
    let sample_name: String     = String::from(sample_index_to_name_map[sample_index].clone());
    let umi_seq_string: String  = std::str::from_utf8(umi_seq).unwrap().to_string();
    records_written_index += 1;
    let read_out_header: String = format!("{}-P5{}-P7{}_{}|{}|{}|{}|{}_{}|{}", sample_name,
                                                                               p5_well_name,
                                                                               p7_well_name,
                                                                               records_written_index,
                                                                               sample_name,
                                                                               p5_well_name,
                                                                               p7_well_name,
                                                                               rt_well_name,
                                                                               lig_well_name,
                                                                               umi_seq_string);

    let barcode_umi_seq: String = format!("{:A>7}{:A>7}{:A>7}{:A>7}{}",
                                          rt_index_encoded,
                                          lig_index_encoded,
                                          p7_index_encoded,
                                          p5_index_encoded,
                                          umi_seq_string);

    if(output_file_format == "fastq") {
      {
        let writer1 = &mut fastq_writer_vec[0][sample_index];
        writer1.write(&read_out_header,
                      None,
                      barcode_umi_seq.as_bytes(),
                      b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC").expect("Error: unable to write sequence to fastq file.");
      }

      {
        let writer2 = &mut fastq_writer_vec[1][sample_index];
        writer2.write(&read_out_header, 
                      None,
                      fastq_record2.seq(),
                      fastq_record2.qual()).expect("Error: unable to write quality values to fastq file.");
      }
    }
    else
    if(output_file_format == "bam") {
    let cell_barcode: String = format!("{:A>7}{:A>7}{:A>7}{:A>7}",
                                       rt_index_encoded,
                                       lig_index_encoded,
                                       p7_index_encoded,
                                       p5_index_encoded);
      let umi_seq: String = format!("{}", std::str::from_utf8(umi_seq).unwrap().to_string());

      let mut record: rust_htslib::bam::Record = rust_htslib::bam::Record::new();

      /*
      ** The BAM record quality values have no offset so subtract
      ** the offset from the input fastq quality values.
      */
      qual_bam.clear();
      for i in (0..fastq_record2.qual().len()) {
        qual_bam.push(fastq_record2.qual()[i] - 33);
      }

      record.set(read_out_header.as_bytes(),
                 None,
                 fastq_record2.seq(),
                 &qual_bam);
      record.push_aux("sS".as_bytes(), rust_htslib::bam::record::Aux::String(&barcode_umi_seq)).expect("Error: unable to add barcode sequence to BAM record tags.");
      record.push_aux("sQ".as_bytes(), rust_htslib::bam::record::Aux::String("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC")).expect("Error: unable to add barcode quality values to BAM record tags.");
      bam_writer_vec[sample_index].write(&record).expect("Error: unable to write record to BAM file.");
    }

  } /* loop end */

  /*
  ** Flush buffers to files.
  */
  if(output_file_format == "fastq") {
    for iread in 0..2 {
      for isample in (0..fastq_writer_vec[0].len()) {
        fastq_writer_vec[iread][isample].flush().expect("Error: unable to flush BAM file to disk.");
      }
    }
  }

  log_read_counts.insert("input_read_pairs".to_string(), nrecord);
  log_read_counts.insert("rt_9_matched".to_string(), count_rt_9);
  log_read_counts.insert("rt_10_matched".to_string(), count_rt_10);
  log_read_counts.insert("lg_9_matched".to_string(), count_lg_9);
  log_read_counts.insert("lg_10_matched".to_string(), count_lg_10);
  log_read_counts.insert("rtlg_9_matched".to_string(), count_rtlg_9);
  log_read_counts.insert("rtlg_10_matched".to_string(), count_rtlg_10);
  log_read_counts.insert("rtlg_9_and_rtlg_10_matched".to_string(), count_rtlg_9_and_rtlg_10);
  log_read_counts.insert("matches_failed".to_string(), count_matches_failed);

  Ok(log)
}


#[derive(Debug, Clone, Serialize)]
struct BarcodeCounter {
  well_index: usize,
  well_name: String,
  whitelist_sequence: String,
  mismatch_sequences: Vec<String>,
  whitelist_read_counts: u64,
  mismatch_read_counts: u64,
  sample_index: usize,
  sample_name: String,
  lanes: Vec<usize>,
  p7_indices: Vec<usize>,
  p5_indices: Vec<usize>
}


#[derive(Debug, Clone, Serialize)]
struct JsonLog {
  input_read_pairs: u64,
  rt_9_matched: u64,
  lg_9_matched: u64,
  rtlg_9_matched: u64,
  rt_10_matched: u64,
  lg_10_matched: u64,
  rtlg_10_matched: u64,
  rtlg_9_and_rtlg_10_matched: u64,
  matches_failed: u64,
  rt_counter_vec: Vec<BarcodeCounter>,
}


/*
** Write log file.
*/
fn write_log_file(rt_barcode_id_map: &HashMap<Vec<u8>, BarcodeIdentifier>,
                  log_read_counts: &HashMap<String, u64>,
                  well_index_to_sample_index_map: &Vec<usize>,
                  sample_index_to_name_map: &Vec<String>,
                  file_indices: &HashMap<String, usize>,
                  log_file_name: &str) -> Result<(), Box<dyn Error>> {

  /*
  ** Gather whitelist sequences.
  */
  let mut rt_whitelist_map: HashMap<String, String> = HashMap::new();

  for key in rt_barcode_id_map.keys() {
    let bi = rt_barcode_id_map.get(key).unwrap();
    if(bi.corrected_flag == false) {
      let well_name = bi.well_name.clone();
      let seq = std::str::from_utf8(key).unwrap().to_string();
      rt_whitelist_map.insert(well_name, seq);
    }
  }

  /*
  ** Gather mismatch sequences and read counts.
  */
  let mut rt_mismatch_map: HashMap<String, Vec<String>> = HashMap::new();

  for key in rt_barcode_id_map.keys() {
    let bi = rt_barcode_id_map.get(key).unwrap();
    if(bi.corrected_flag == true) {
      let well_name = bi.well_name.clone();
      let read_count = bi.read_count;
      let entry = rt_mismatch_map.entry(well_name).or_insert_with(Vec::new);
      entry.push(format!("{}:{}", std::str::from_utf8(key).unwrap().to_string(), read_count));
    }
  }

  let mut barcode_counter_map: HashMap<String, BarcodeCounter> = HashMap::new();
  for key in rt_barcode_id_map.keys() {
    let bi = rt_barcode_id_map.get(key).unwrap();
    let well_index = bi.well_index;
    let well_name = bi.well_name.clone();
    if(!barcode_counter_map.contains_key(&well_name)) {
      let sample_index = well_index_to_sample_index_map[bi.well_index];
      let sample_name = sample_index_to_name_map[sample_index].clone();
      let _ = barcode_counter_map.insert(well_name.clone(), BarcodeCounter {well_index: bi.well_index,
                                                                            well_name: well_name.clone(),
                                                                            whitelist_sequence: rt_whitelist_map[&well_name].clone(),
                                                                            mismatch_sequences: rt_mismatch_map[&well_name].clone(),
                                                                            whitelist_read_counts: 0_u64,
                                                                            mismatch_read_counts: 0_u64,
                                                                            sample_index: sample_index,
                                                                            sample_name: sample_name,
                                                                            lanes: vec![file_indices["lane"]],
                                                                            p7_indices: vec![file_indices["p7"]],
                                                                            p5_indices: vec![file_indices["p5"]]});
    }
    let bc = barcode_counter_map.get_mut(&well_name).unwrap();
    if(!bi.corrected_flag) {
      bc.whitelist_read_counts += bi.read_count;
    }
    else {
      bc.mismatch_read_counts += bi.read_count;
    }
  }

  /*
  ** Serialize data.
  */
  let mut rt_counter_vec: Vec<BarcodeCounter> = Vec::new();
  let mut keys: Vec<String> = barcode_counter_map.keys().map(|i| i.clone()).collect();
  keys.sort();
  for key in keys.iter() {
    rt_counter_vec.push(barcode_counter_map[key].clone());
  }

  let json_log = JsonLog { input_read_pairs:           log_read_counts["input_read_pairs"],
                           rt_9_matched:               log_read_counts["rt_9_matched"],
                           lg_9_matched:               log_read_counts["lg_9_matched"],
                           rtlg_9_matched:             log_read_counts["rtlg_9_matched"],
                           rt_10_matched:              log_read_counts["rt_10_matched"],
                           lg_10_matched:              log_read_counts["lg_10_matched"],
                           rtlg_10_matched:            log_read_counts["rtlg_10_matched"],
                           rtlg_9_and_rtlg_10_matched: log_read_counts["rtlg_9_and_rtlg_10_matched"],
                           matches_failed:             log_read_counts["matches_failed"],
                           rt_counter_vec:             rt_counter_vec};

  let writer_json = std::fs::File::create(log_file_name)?;
  serde_json::to_writer_pretty(writer_json, &json_log)?;

  Ok(())
}


fn main() {
  /*
  ** Process command line options.
  */
  let cl_options = set_cl_options().unwrap();
  let cl_arg = cl_options.get_matches();

  /*
  ** Get command line argument values.
  */
  let fastq1_file: String           = cl_arg.get_one::<String>("fastq_read_1").unwrap().to_string();
  let fastq2_file: String           = cl_arg.get_one::<String>("fastq_read_2").unwrap().to_string();
  let samplesheet_file: String      = cl_arg.get_one::<String>("samplesheet_file").unwrap().to_string();
  let default_rt_file: String       = cl_arg.get_one::<String>("rt_barcode_file").unwrap().to_string();
  let default_ligation_file: String = cl_arg.get_one::<String>("ligation_barcode_file").unwrap().to_string();
  let output_file_format: String    = cl_arg.get_one::<String>("output_file_format").unwrap().to_string();
  let num_threads: usize            = cl_arg.get_one::<usize>("number_bam_threads").unwrap().clone();
  let uncompress_flag: bool         = cl_arg.get_flag("uncompressed_fastqs");

  println!("");
  println!("fastq1:           {}", fastq1_file);
  println!("fastq1:           {}", fastq2_file);
  println!("samplesheet:      {}", samplesheet_file);
  println!("rt default:       {}", default_rt_file);
  println!("ligation default: {}", default_ligation_file);
  println!("out format:       {}", output_file_format);
  println!("num threads:      {}", num_threads);
  println!("uncompressed:     {}", uncompress_flag);
  println!("");

  
  /*
  ** Get lane index from fastq filename.
  */
  let file_indices: HashMap<String, usize> = get_file_indices(&fastq1_file).unwrap();
  println!("file indexes: {:#?}", file_indices);

  /*
  ** For each barcode make a whitelist of barcode sequences, which includes
  ** a pointer, or equivalent, to the barcode index.
  **   o  determine required files: use either names in
  **      samplesheet file or the default names on the
  **      command line.
  **   o  load barcode data from file
  **   o  sort barcode data by indices (the indices are determined using
  **      the well to index map)
  **   o  expand the whitelist to include discrepancies
  **   o  notes:
  **        o  I want fast access to the sequences, especially the original
  **           sequences
  **        o  I want to find the barcode index and well name quickly:
  **        o  I need an rt barcode to sample_id map in order to find
  **           the fastq file to which the read is written.
  **             o  potentially store reads in fastq files by either sample_id
  **                or well index
  **        o  I need the rt and ligation well names for the read name
  **        o  store the barcode sequences in hashmaps in order to
  **           find barcode sequences from read subsequences
  **        o  sequence test order
  **             o  long ligation vs long ligation barcodes, if pass
  **                  o  test matched long ligation with rt barcode, if pass
  **                       o  done
  **                     if fail
  **                  o  test match long ligation with rt barcode w/discreps, if pass
  **                       o  done
  **             o  long ligation vs long ligaton barcodes w/discrep, if pass
  **                  o  test matched long ligation w/discrep with rt barcode, if pass
  **                       o  done
  **                  o  test matched long ligatin w/discrep with rt barcdoe w/discrep, if pass
  **                       o  done
  **             o if not pass, repeat above with short ligation barcodes
  **        o  Andrew/Hannah demux
  **             o  tests 9 base ligation and rt barcodes, and 10 base ligation
  **                and rt barcode
  **                  o  if both correct, count as ambiguous
  **                  o  else if one or the other correct, accept
  **                  o  else bad barcode, reject as undetermined
  **             o  python script selects read1 vs read2 by spec per read.
  **                This is necessary for the python versions because it
  **                demultiplexes by PCR indices as well but not for the
  **                new demultiplexer, which does not demultiplex PCR
  **                indices.
  **             o  the spec dictionary provides greater generality than is
  **                used by the demultiplexing script so I am inclined to
  **                reduce the generality for the benefit of simpler
  **                implementation and greated execution speed.
  **             o  how do I find a reasonably 'safe' balance between flexibility,
  **                simplicity, and execution speed?
  **             o  the spec provides specifying the following
  **                  o  barcode sequence name (e.g., rt_9, rt_10, ligation_9, ...)
  **                  o  read 1/2 where barcode is found
  **                  o  barcode sequence start and end in read sequence
  **                  o  whether the read sequence is checked by comparing
  **                     to a barcode sequence whitelist.
  **                       o  RT, ligation, P7, and P5 sequences are compared
  **                          with barcode whitelists
  **                       o  UMI sequence is not compared and corrected
  **             o  the fastq iterator returns the (corrected) barcode
  **                or None if uncorrectable) and UMI sequences
  **
  **  o  data structure for barcode 'lookup' using read substring as HashMap key
  **       o  notes:
  **            o  use bio::io::fastq for reading and writing fastq files (this
  **               works with multiline sequences and quality values).
  **            o  read seq and qual are &[u8] slices when using reader.read(&record)
  **            o  read substrings are &[u8] slices
  **
  **       
  */
  
  /*
  **  Make index encoder.
  **  Use the index encoder to construct the R1 read sequences from
  **  the rt, ligation, p7, p5 indices. This is necessary to get
  **  around STARsolo's 31 barcode sequence limit.
  */
  if(MAX_NUM_PLATES * 96 >= 4_usize.pow(7)) {
    eprintln!("Error: the maximum number of wells exceeds the largest well");
    eprintln!("       index encodable as a string of bases. You must");
    eprintln!("       decrease the value of MAX_NUM_PLATES in the");
    eprintln!("       rna_rtlig_demux program source file main.rs." );
  }
  let index_encoder = make_index_encoder((MAX_NUM_PLATES*96) as u64).unwrap();
/*
  for (ii, idx) in index_encoder.into_iter().enumerate() {
    println!("{}  {:A>7}", ii, idx);
  }
*/
  // Read samplesheet JSON file.
  let samplesheet_json: serde_json::Value = get_samplesheet_json(samplesheet_file).unwrap();

  // Deserialize JSON sample_map_list from serde_json Value.
  let sample_map_vec_all: Vec<SampleMap> = deserialize_sample_map_vector(samplesheet_json.clone()).unwrap();

  let sample_map_vec: Vec<SampleMap> = get_lane_samples(&sample_map_vec_all, file_indices["lane"], file_indices["p7"], file_indices["p5"]).unwrap();

  /*
  ** Make a lane to barcode file map with fastq output file indices.
  */
  let recipe = "std_1";
  let mut rt_barcode_id_map = make_barcode_id_map(&sample_map_vec, "rt_file", &default_rt_file, recipe).unwrap();
  let mut ligation_barcode_id_map = make_barcode_id_map(&sample_map_vec, "ligation_file", &default_ligation_file, recipe).unwrap();

  /*
  ** Make sample-related index maps.
  **   o  well_index_to_sample_index_map: Vec<usize>: maps RT well index to sample index.
  **   o  sample_index_to_name_map: Vec<String>: mapes sample index to sample name string.
  */
  let (well_index_to_sample_index_map, sample_index_to_name_map) = make_index_maps(&sample_map_vec, &mut rt_barcode_id_map).unwrap();

/*
  println!("well_index_to_sample_index_map:");
  for (i, index) in well_index_to_sample_index_map.iter().enumerate() {
    println!("{}:  {:#?}", i, well_index_to_sample_index_map[i]);
  }
  println!("");

  println!("sample_index_to_name_map:");
  println!("{:#?}", sample_index_to_name_map);
  println!("");
*/

  /*
  ** Process fastq reads.
  */
  let mut log_read_counts: HashMap<String, u64> = HashMap::new();
  let _ = process_reads(&fastq1_file,
                        &fastq2_file,
                        &mut rt_barcode_id_map,
                        &mut ligation_barcode_id_map,
                        &file_indices,
                        &recipe,
                        &index_encoder,
                        &output_file_format,
                        uncompress_flag,
                        num_threads,
                        &well_index_to_sample_index_map,
                        &sample_index_to_name_map,
                        &mut log_read_counts);

  /*
  ** Write log file.
  */
  let log_file_name: String = String::from(format!("rna_rtlig_demux_log.{}_{}_{}.json", file_indices["lane"], file_indices["p7"], file_indices["p5"]));
  write_log_file(&rt_barcode_id_map,
                 &log_read_counts,
                 &well_index_to_sample_index_map,
                 &sample_index_to_name_map,
                 &file_indices,
                 &log_file_name).unwrap();
}


