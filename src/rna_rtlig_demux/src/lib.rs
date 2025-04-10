/*
** Turn off warnings about unused parentheses.
*/
#![allow(unused_parens)]


//! Reproduce Andrew Hill's barcode_utils.py in Rust, a bit at a time.
pub mod barcode_utils {

  use std::io;
  use std::error::Error;
  use std::collections::HashMap;
  use itertools::Itertools;
  use csv::{ReaderBuilder, Trim};
  use serde::Deserialize;


  /*
  ** Convert between UTF-8 and String, etc.
  ** https://stackoverflow.com/questions/41034635/how-do-i-convert-between-string-str-vecu8-and-u8
  */

  /// Convert an &str to &[u8].
  ///
  /// Arguments:
  ///- in_rstr: a string slice to convert.
  ///
  /// Return:
  ///
  /// A vector of ASCII character codes.
  ///
  #[inline]
  pub fn str_to_u8(in_rstr: &str) -> &[u8] {
    return(in_rstr.as_bytes()) 
  }


  /// Convert &[u8] to &str.
  ///
  /// Arguments:
  ///- in_ru8: a u8 vector slice.
  ///
  /// Return:
  ///
  /// An &str.
  ///
  #[inline]
  pub fn u8_to_str(in_ru8: &[u8]) -> Result<&str, std::str::Utf8Error> {
    return(std::str::from_utf8(in_ru8));
  }


  /// Read a barcode TSV file that has the format
  ///
  /// barcode name\tbarcode sequence
  ///
  /// Arguments:
  ///- file_path: a reference to the input file path
  ///
  /// Return:
  ///
  /// A hash map: the key is the barcode sequence and the value is
  /// the barcode name. The barcode sequence may be String or Vec<u8>,
  /// choose below.
  ///
  /// Based on URL: https://stackoverflow.com/questions/78639668/fast-reading-from-a-tsv-file-in-rust
  #[derive(Deserialize)]
  struct Record {
    hash_name: String,
    hash_barcode: String,
  }

  pub fn read_barcode_file(file_path: &str) -> Result<HashMap<String, String>, Box<dyn Error>> {
    println!("read barcode file: {}", file_path);
    let fp = std::fs::File::open(file_path)?;
    let buf_reader = io::BufReader::new(fp);
    let mut tsv_reader = ReaderBuilder::new()
                           .has_headers(false)
                           .trim(Trim::Fields)
                           .delimiter(b'\t')
                           .comment(Some(b'#'))
                           .from_reader(buf_reader);

    let mut hash_map: HashMap<String, String> = HashMap::new();

    /*
    ** Read the file and store the names and sequences.
    */
    for result in tsv_reader.deserialize() {
       let record: Record = result?;
       if(hash_map.contains_key(&record.hash_barcode)) {
         panic!("Error: read_barcode_file: duplicate barcode string in file {:#?}", file_path);
       }
       hash_map.entry(record.hash_barcode.to_owned()).or_insert(record.hash_name.to_string());
    }

    Ok(hash_map)
  }


  /// Generate a vector of mismatched sequences to a given sequence. Must contain only ACGT.
  /// This is based heavily on a Biostars answer.
  /// Curiously, it returns the input sequence when num_mismatches is zero; otherwise, it does not.
  ///
  /// Arguments:
  ///- sequence: a reference to the input string.
  ///- num_mismatches: the number of mismatches in each output sequence.
  ///- allow_n: include Ns in the output sequences.
  ///
  /// Return:
  ///
  /// A Result enum containing a vector of mismatched sequence Strings.
  ///
  pub fn generate_mismatches(sequence: &str, num_mismatches: usize, allow_n: bool) -> Result<Vec<String>, Box<dyn Error>> {

    let sequence_upper: String = sequence.to_uppercase();

    if(num_mismatches == 0) {
      return(Ok(vec!(sequence_upper)));
    }

    let mut letters = String::from("ACGT");
    if(allow_n == true) {
      letters.push('N');
    }

    /*
    ** Return a Result of vector of strings.
    */
    let mut mismatches: Vec<String> = Vec::new();

    let seq_len = sequence_upper.len();

    /*
    ** Make iterators of vectors of indices where mismatches
    ** may show up in the input sequence, and iterator through
    ** these vectors.
    */
    let combination_iter = (0..seq_len).combinations(num_mismatches);
    for locs in combination_iter {
      /*
      ** sequence_vector begins with a vector of vectors where each base
      ** in the sequence is an inner vector. That is, for sequence ATGCTA,
      ** each loop iteration gives the save vector of vectors, at this point.
      **   sequence_vector: [['A'], ['T'], ['G'], ['C'], ['T'], ['A']]
      **   ...
      */
      let mut sequence_vector: Vec<Vec<char>> = sequence_upper.chars().map(|c| vec!(c)).collect();

      /*
      ** Here inner vectors are replaced with the three bases that are not
      ** in the original sequence. For example,
      **   [['A'], ['T'], ['G'], ['C'], ['T'], ['A']]
      ** becomes
      **   [['A'], ['A','C','G'], ['G'], ['C'], ['T'], ['A']]
      ** at mismatch index 1, which is the second inner
      ** vector.
      */
      for loc in locs {
        let orig_char = sequence_upper.chars().nth(loc as usize).unwrap();
        sequence_vector[loc] = letters.chars().filter(|c| *c != orig_char).collect();
      }

      /*
      ** And the cartesian product of the set of inner vectors expands
      ** to the desired mismatch sequences when the product vectors are
      ** converted to strings.
      */
      for vector_set in sequence_vector.iter().multi_cartesian_product() {
        let primer_string: String = vector_set.into_iter().collect();
        mismatches.push(primer_string);
      }
    }

    Ok(mismatches)
  }


  /// Construct a precomputed set of all mismatches within a specified
  /// edit distance and the barcode whitelist.
  ///
  /// Arguments:
  ///- whitelist: sequences to expand as a vector of strings.
  ///- edit_distance: maximum number of substitutions per sequence.
  ///- allow_n: include Ns as substitutions.
  ///
  /// Return:
  ///
  /// A result enum containing a vector of maps of mismatched sequences to
  /// their whitelist sequences.
  ///
  pub fn construct_mismatch_to_whitelist_map(whitelist: Vec<String>, edit_distance: usize, allow_n: bool) -> Result<Vec<HashMap<String, String>>, Box<dyn Error>> {

    /*
    ** Set whitelist sequences to upper-case.
    */
    let whitelist_upper = whitelist.iter().map(|s| s.to_uppercase()).collect::<Vec<String>>();

    /*
    ** mismatch_to_whitelist_map is a vector of hash maps where the vector
    ** index is the number of substitutions in the mismatch sequence, and
    ** the hash maps are keyed by the sequences with substitutions and the
    ** values are the original whitelist sequences (no mismatches).
    */
    let mut mismatch_to_whitelist_map: Vec<HashMap<String, String>> = Vec::with_capacity(edit_distance+1);
    for _i in (0..edit_distance+1) {
      mismatch_to_whitelist_map.push(HashMap::new());
    }

    /*
    ** Set the zero mismatch sequence maps where the key and
    ** value are the same.
    */
    mismatch_to_whitelist_map[0] = whitelist_upper.iter().map(|s| (s.to_owned(), s.to_owned())).collect::<HashMap<String, String>>();


    /*
    ** Track  conflicts where mismatches map to different sequences.
    */
    let mut conflicting_mismatches: Vec<String> = Vec::new();

    /*
    ** Doesn't really matter as correction function will never see it,
    ** but exclude any perfect matches to actual seqs by mismatches.
    */
    conflicting_mismatches.extend(whitelist_upper.clone());


    for mismatch_count in (1..edit_distance+1) {

      for sequence in &whitelist_upper {
        /*
        ** Generate all possible mismatches in range.
        */
        let mismatches = generate_mismatches(&sequence, mismatch_count, allow_n)?;

        for mismatch in mismatches.iter() {
          if(mismatch_to_whitelist_map[mismatch_count].contains_key::<str>(&mismatch) == true) {
            conflicting_mismatches.push(mismatch.to_string());
          }
          mismatch_to_whitelist_map[mismatch_count].insert(mismatch.to_string(), sequence.clone());
        }
      }

      /*
      ** Go back and remove any conflicting mismatches.
      */
      for mismatch in conflicting_mismatches.clone().into_iter().unique() {
        if(mismatch_to_whitelist_map[mismatch_count].contains_key::<str>(&mismatch) == true) {
          mismatch_to_whitelist_map[mismatch_count].remove(&mismatch);
        }
      }
    }

    Ok(mismatch_to_whitelist_map)
  }


  /// Convert the output sequences of construct_mismatch_to_whitelist_map()
  /// from String to Vec<u8> in order to avoid UTF-8 checking. My tests
  /// suggest that using Vec<u8> rather than String may reduce the run
  /// time by about 8% to 9% where the run time includes reading the fastq
  /// input time and building the hashdict. My impression is that there is
  /// no compelling reason to use references because Vecs and Strings are
  /// stored as fat pointers with the strings/sequences being stored on the
  /// heap.
  ///
  /// Argument:
  ///- whitelist_map_string: the map list using Strings for sequences. This is returned by construct_mismatch_to_whitelist_map(), above.
  ///
  /// Return:
  ///
  /// A map list using Vec<u8> for sequences.
  ///
  pub fn mismatch_to_whitelist_map_as_u8(whitelist_map_string: Vec<HashMap<String, String>>) -> Result<Vec<HashMap<Vec<u8>, Vec<u8>>>, Box<dyn Error>> {
    let vec_len: usize = whitelist_map_string.len();
    let mut whitelist_map_u8: Vec<HashMap<Vec<u8>, Vec<u8>>> = Vec::with_capacity(vec_len);
    for i in 0..vec_len {
      whitelist_map_u8.push(HashMap::new());
      for key in whitelist_map_string[i].keys() {
        whitelist_map_u8[i].insert(key.as_bytes().to_vec(), whitelist_map_string[i][key].as_bytes().to_vec());
      }
    }

    Ok(whitelist_map_u8)
  }


  /// Convert a well index to a well name string. The string does not
  /// include the plate substring; however, an integer identifying the
  /// plate is returned.
  ///
  /// Argument:
  ///- well_index: (usize) the index to convert to a well name starts
  ///              with 1 (A01) and ends with 96 (H12) for plate 1.
  ///- across_row_first: (bool) do the indices increase by one
  ///  across the row (true) or down the column (false).
  ///
  /// Return:
  ///
  ///   A Result containing a tuple with the values (<iplate>, <well name>)
  ///   where iplate is one for the first plate.
  ///
  pub fn index_to_well(mut well_index: usize, across_row_first: bool) -> Result<(usize, String), Box<dyn Error>> {

    /*
    ** well_index == 0 is none.
    */
    if(well_index < 1) {
      return( Ok((0_usize, String::from("none"))) );
    }

    well_index = well_index - 1;
    let nrow: usize = 8;
    let ncol: usize = 12;
    let ipl: usize = well_index / 96;
    let i96: usize = well_index - ipl * 96;
    let well_row: char;
    let well_col: usize;
    if(across_row_first) {
      well_row = char::from_u32((65 + i96 / ncol).try_into().unwrap()).expect("unable to convert usize to u8");
      well_col = (i96 % ncol) + 1;
    }
    else {
      well_row = char::from_u32((65 + (i96 % nrow)).try_into().unwrap()).expect("unable to convert usize to u8");
      well_col = i96 / nrow + 1;
    }

    let well_id: String = format!("{}{:02}", well_row, well_col);

    /*
    ** Return tuple of ipl and well_id
    */
    Ok((ipl + 1, well_id))
  }


  /// Convert a well name string to a well index.
  ///
  /// Argument:
  ///- plate: (usize) the plate number beginning with 1.
  ///- row: (char) the row identifier a-h in lower- or
  ///       upper-case.
  ///- column: (usize) the column identifer in the range 1 to 12.
  ///- across_row_first: (bool) do the indices increase by one
  ///  across the row (true) or down the column (false).
  ///
  /// Return:
  ///   index where well A01 is index 1.
  ///
  pub fn well_to_index(iplate: usize, row: char, column: usize, across_row_first: bool) -> Result<usize, String> {
    let row_values: Vec<char> = vec!['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'];
    let crow: char = row.to_lowercase().next().expect("failed conversion to lowercase");
    if(!row_values.iter().any(|&c| c == crow)) {
      //return(Err(String::from("not a valid row")));
      return(Err(format!("{} is not a valid row", row)));
    }
    if(column < 1 || column > 12) {
      return(Err(format!("{} is not a valid column", column)));
    }
    let irow = (crow as usize) - ('a' as usize);
    let icol: usize = column - 1;
    let well_index: usize;
    if(across_row_first) {
      well_index = irow * 12 + icol + 1 + (iplate - 1)  * 96;
    }
    else {
      well_index = icol * 8 + irow + 1 + (iplate - 1) * 96;
    }

    Ok(well_index)
  }

  /// Reverse complement a DNA sequence with bit-twiddling
  ///
  /// # Example
  ///
  /// ```rust
  /// let seq = b"ATGCTTCCAGAA";
  ///
  /// let actual = revcomp(seq);
  /// let expected = b"TTCTGGAAGCAT";///
  /// assert_eq!(actual, expected)
  /// ```
  ///
  /// # Notes
  ///   Implementation is taken from https://doi.org/10.1101/082214
  ///   The sequence may be in upper or lower case and consists of
  ///   only the four possible bases.
  ///
  pub fn revcomp_2(seq: &[u8]) -> Vec<u8> {
      seq.iter()
         .rev()
         .map(|c| if c & 2 != 0 { c ^ 4 } else { c ^ 21 })
         .collect()
  }


  // Set DNA complement hashmap once.
  static GCMP: std::sync::LazyLock<HashMap<u8, u8>> = std::sync::LazyLock::new(|| HashMap::from([
        (97,  116),
        (99,  103),
        (103, 99),
        (116, 97),
        (110, 110),
        (65,  84),
        (67,  71),
        (71,  67),
        (84,  65),
        (78,  78)]));

  /// Reverse complement a DNA sequence given as a slice
  /// of u8 ASCII values.
  ///
  ///
  /// Argument:
  ///- seq: a u8 vector slice.
  ///
  /// Return:
  ///   A u8 vector with the reverse complement
  ///   of the input.
  ///
  pub fn revcomp(seq: &[u8]) -> Vec<u8> { 
      seq.iter()
         .rev()
         .map(|c| {GCMP[c]} )
       .collect()
  }


/*
  /// Validate a barcode spec.

  /// Load a barcode spec.

  /// Validate a whitelist.

  /// class WriteBuffer

  /// class BarcodeCorrector
  ///   method init
  ///   method correct
  ///   method get_min_hamming: returns a list of the minimum N hamming distances observed
  ///   get barcode length

  /// get index coordinates

  /// validate barcode read pair

  /// generate well ids

  /// parse fastq barcodes
*/


}

