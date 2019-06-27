#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <map>
#include <getopt.h>
#include <math.h>
#include <omp.h>
#include "kmer_library.h"
#include "cpas_debug.h"

using namespace std;

void error(const char* msg) {
  fprintf(stderr, "ERROR: %s\n", msg);
  exit(2);
}

void printUsageAndExit() {
  fprintf(stderr, "Usage: count_kmer [options] <FASTA file name> <SAM file name>\n\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t--kmer,-k <SIZE>\tSpecifies the k-mer size\n");
  fprintf(stderr, "\t--csv           \tOutput csv to STDOUT\n");
  fprintf(stderr, "\n");
  exit(2);
}

template<uint KMERSIZE>
struct FrequencyTable {
  static const size_t contextsize = ipow(5, KMERSIZE * 2 + 1);
  typedef int Frequency;
  typedef int Score;
  vector<Frequency> err_table;
  vector<Frequency> err_context_table;
  vector<double>    err_context_normalized_table;
  //vector<Score>     score_table;
  inline Frequency& ect(size_t r, size_t q) {
    MYASSERT_WMD("Out of range (r)", r < 5, DUMP(r));
    MYASSERT_WMD("Out of range (q)", q < contextsize, DUMP(q));
    return err_context_table[r * contextsize + q];
  }
  inline double&   ecnt(size_t r, size_t q) {
    MYASSERT_WMD("Out of range (r)", r < 5, DUMP(r));
    MYASSERT_WMD("Out of range (q)", q < contextsize, DUMP(q));
    return err_context_normalized_table[r * contextsize + q];
  }
/*
  void outputAsBinaryTable(const string& outputFileName)
  {
    const char *fn = outputFileName.c_str();
    FILE *fle = fopen(fn, "wb");
    if(fle == NULL) {
      fprintf(stderr, "ERROR: Cannot open binary table file '%s'\n", fn);
      exit(2);
    }
    //// FILE FORMAT:
    ////  offset  0: 'BINFREQT' (8 bytes)
    ////  offset  8: k-mer size (8 bytes)
    ////  offset 16: sizeof(Frequency) (8 bytes)
    ////  offset 24: table_size (8 bytes)
    ////  offset 32: table (8 x table_size bytes)
    static const char MAGIC_STRING[] = "BINFREQT";
    MYASSERT_WMD("MAGIC_STRING must be 8 bytes", sizeof(MAGIC_STRING) - 1 == 8, DUMP(MAGIC_STRING));
    fwrite(MAGIC_STRING, sizeof(MAGIC_STRING) - 1, 1, fle);
    size_t kmer_size = KMERSIZE;
    MYASSERT_WMD("sizeof(size_t) must be 8", sizeof(size_t) == 8, DUMP(sizeof(size_t)));
    fwrite(&kmer_size, sizeof(kmer_size), 1, fle);
    size_t size_of_frequency = sizeof(Frequency);
    fwrite(&size_of_frequency, sizeof(size_of_frequency), 1, fle);
    size_t var_table_size = contextsize * contextsize;
    fwrite(&var_table_size, sizeof(var_table_size), 1, fle);
    //MYASSERT_WMD("var_table_size == kmer_kmer_table.size()", var_table_size == kmer_kmer_table.size(), DUMP(var_table_size, contextsize, kmer_kmer_table.size()));
    fwrite(&*score_table.begin(), sizeof(Frequency) * score_table.size(), 1, fle);
    fclose(fle);
  }
*/
  void countAlignment(
    const BString& ras,
    const BString& qas
  )
  {
    MYASSERT_WMD("ras.size() must be qas.size()", ras.size() == qas.size(), DUMP(ras.size(), qas.size()));
    const size_t aligned_len = ras.size();
    if(aligned_len < KMERSIZE) return;
    KInt<2 * KMERSIZE + 1> con_idx;
    KInt<1>                err_idx;
    KInt<KMERSIZE>         pre;
    KInt<KMERSIZE>         nxt;
    for(size_t i = 0; i < KMERSIZE - 1; ++i) {
      con_idx.ShiftIn(ras[i]);
      err_idx.ShiftIn(qas[i]);
    }
    for(size_t i = KMERSIZE - 1; i <= aligned_len - KMERSIZE; i++){
      if(ras[i] != qas[i]){
        err_idx.ShiftIn(qas[i]);
        //pre
        size_t diff = 0;
        size_t shiftcount = 1;
        while(shiftcount <= KMERSIZE){
          if(ras[i - shiftcount - diff] != char2Base('-')){
            pre.unshift(ras[i - shiftcount - diff]);
            shiftcount++;
          }else{
            diff++;
          }
        }
        diff = 0;
        shiftcount = 1;
        //nxt
        while(shiftcount <= KMERSIZE){
          if(ras[i + shiftcount + diff] != char2Base('-')){
            nxt.ShiftIn(ras[i + shiftcount + diff]);
            shiftcount++;
          }else{
            diff++;
          }
        }
        con_idx = pre * ipow(5, KMERSIZE + 1) + ras[i] * ipow(5,KMERSIZE) + nxt;
        err_table[con_idx]++;
        ect(err_idx, con_idx)++;
      }
    }
  }

  size_t score_count = 1;

  void scorerize(const int diff){
    #pragma omp parallel num_threads(4)
    {
      #pragma omp for
      for(int i = 0; i < 5; i++){
        for(int j = 0; j < contextsize; j++){
          if(score_count % 10000 == 0) fprintf(stderr, "first loop : %'d / %'d\r", score_count, contextsize * 24);
          #pragma omp atomic
          score_count++;
          if(err_table[j] != 0){
            ecnt(i, j) = static_cast<double>( ect(i, j) ) / err_table[j];
          }else{
            ecnt(i, j) = 0.0;
          }
          MYASSERT_WMD("prob must be in [0, 1]", ecnt(i, j) <= 1.0 && ecnt(i, j) >= 1.0, DUMP(ecnt(i, j)));
          }
        }
      }
    fprintf(stderr, "first loop done                                                \n");
    }

  void printet(){
    KInt<KMERSIZE * 2 + 1> tmp;
    for(int j = 0; j < contextsize; j++){
      tmp = j;
      fprintf(stdout, "%11s", tmp.str().c_str());
      if(j != contextsize - 1){
        fprintf(stdout, ", ");
      }else{
        fprintf(stdout, "\n");
      }
    }
    for(int i = 0; i < contextsize; i++){
      fprintf(stdout, "%11d", err_table[i]);
      if(i != contextsize - 1){
        fprintf(stdout, ", ");
      }else{
        fprintf(stdout, "\n");
        }
    }
  }
  void printect(){
    KInt<KMERSIZE * 2 + 1> tmp;
    for(int j = 0; j < contextsize; j++){
      tmp = j;
      fprintf(stdout, "%11s", tmp.str().c_str());
      if(j != contextsize - 1){
        fprintf(stdout, ", ");
      }else{
        fprintf(stdout, "\n");
      }
    }
    for(int i = 0; i < 5; i++){
      for(int j = 0; j < contextsize; j++){
        fprintf(stdout, "%11d", ect(i, j));
        if(j != contextsize - 1){
          fprintf(stdout, ", ");
        }else{
          fprintf(stdout, "\n");
        }
      }
    }
  }

  void printecnt(){
    KInt<KMERSIZE * 2 + 1> tmp;
    for(int j = 0; j < contextsize; j++){
      tmp = j;
      fprintf(stdout, "%11s", tmp.str().c_str());
      if(j != contextsize - 1){
        fprintf(stdout, ", ");
      }else{
        fprintf(stdout, "\n");
      }
    }
    for(int i = 0; i < 5; i++){
      for(int j = 0; j < contextsize; j++){
        fprintf(stdout, "%11lf", ecnt(i, j));
        if(j != contextsize - 1){
          fprintf(stdout, ", ");
        }else{
          fprintf(stdout, "\n");
        }
      }
    }
  }


  FrequencyTable() : err_table(contextsize, 0), err_context_table(5 * contextsize, 0), err_context_normalized_table(5 * contextsize, 0) {}
  void countKmerFrequencies (
    const char* FASTAFileName,
    const char* SAMFileName,
    uint KmerSize,
    const bool outputInCSV,
    const string binaryOutputFileName ///< empty() if --binary is not given
  )
  {
    fprintf(stderr, "\n===Parameters===\n");
    fprintf(stderr, "Reference FASTA: %s\n", FASTAFileName);
    fprintf(stderr, "Input SAM file : %s\n", SAMFileName);
    fprintf(stderr, "K-mer size     : %d\n", KmerSize);
    fprintf(stderr, "================\n");
    fprintf(stderr, "\n");

    fprintf(stderr, "Loading from the FASTA file ...\n");
    const MultiFASTA multiFASTA = loadFromFASTA(FASTAFileName);
    fprintf(stderr, "Done                           \n");
    cerr << "reference record are as below\n---------" << endl;
    for(auto itr = multiFASTA.begin(); itr != multiFASTA.end(); itr++) {
      fprintf(stderr, "%s\n", itr->first.c_str());
    }
    cerr << "---------" << endl;
    FastTSVParse ftp(SAMFileName);
    if(!ftp) {
      fprintf(stderr, "ERROR: Cannot open SAM file '%s'\n", SAMFileName);
      exit(2);
    }

    SAMRecord record;
    size_t recordCount = 0;
    while(ftp.readNextLine()){
      const char* line = ftp.c_str();
      const bool isEmptyLine = line[0] == '\0';
      if(isEmptyLine) continue;
      const bool isCommentLine = line[0] == '@';
      if(isCommentLine) continue;
      if(!record.fill(ftp)) {
        cerr << "SAM Parsing ERROR at line " << ftp.getLineNumber();
        exit(2);
      }
      if(record.rname == "*") continue;
      if(record.seq == "*") continue;
      // use only primary alignment and supplementary alignment and reverse compliment of them.
      // sample has many supplimentary alignment.
      if((record.flag & 2064) != record.flag) continue;


      if(!multiFASTA.count(record.rname.c_str())) {
        cerr << "SAM record says RNAME = '" << record.rname << "', but the reference genome does not have '" << record.rname << "'" << endl;
        continue;
        //exit(2);
      }


      // get aligned sequence at here
      const CIGAROPS cops     = parseCIGARString(record.cigar);

      BString refBS           = multiFASTA.at(record.rname);
      BString queryBS         = String2BString(record.seq);
      const int refStartPos   = record.pos;
      const int queryStartPos = 0; // generateAlignmentSequencesFromCIGARAndSeqs() will manege first Softclip / Hardclip
      BString ras, qas;
      generateAlignmentSequencesFromCIGARAndSeqs(refBS, queryBS, cops, refStartPos, queryStartPos, ras, qas);

      // if record flag has revcomp flag, modify aligned reference sequence and read sequence.
      if((record.flag & 16) != 16){
        revCompBString(ras);
        revCompBString(qas);
      }
      countAlignment(ras, qas);
      ++recordCount;
      if(recordCount % 100 == 0) {
        cerr << recordCount << " processed\r" << flush;
      }
    }
    cerr << recordCount << " processed\n";
    scorerize(100);
    if(!binaryOutputFileName.empty()) {
      //outputAsBinaryTable(binaryOutputFileName);
    }
    if(outputInCSV) {
      //printet();
      //printect();
      printecnt();
    }
  }
};


int main(int argc, char *argv[]){
  GDB_On_SEGV g(argv[0]);

  struct option longopts[] = {
    { "csv"       , no_argument       , NULL , 'c' } ,
    // { "delete" , optional_argument , NULL , 'd' } ,
    { "kmer"      , required_argument , NULL , 'k' } ,
    { "binary"    , required_argument , NULL , 'b' } ,
    { 0           , 0                 , 0    , 0  }  ,
  };

  /// PARAMETERS ///
  int kmer_size = 1;
  bool output_in_csv = false;
  string binary_output_file_name;
  //////////////////
  int opt;
  int longindex;
  while ((opt = getopt_long(argc, argv, "k:", longopts, &longindex)) != -1) {
    switch (opt) {
    case 'b':
      binary_output_file_name = optarg;
      break;
    case 'c':
      output_in_csv = true;
      break;
    case 'k':
      kmer_size = atoi(optarg);
      if(kmer_size < 1 || 10 < kmer_size) {
        error("kmer size must be 1-10");
      }
      break;
    default:
      MYASSERT_NEVERREACH();
    }
  }
  const int NUM_REQUIRED_ARGUMENTS = 2;
  if(optind + NUM_REQUIRED_ARGUMENTS != argc) {
    printUsageAndExit();
  }
  const char* fasta_file_name = argv[optind + 0];
  const char* sam_file_name   = argv[optind + 1];

  #define FT() ft.countKmerFrequencies(fasta_file_name, sam_file_name, kmer_size, output_in_csv, binary_output_file_name)
  switch(kmer_size) {
    case  1: { FrequencyTable< 1> ft; FT(); } break;
    case  2: { FrequencyTable< 2> ft; FT(); } break;
    case  3: { FrequencyTable< 3> ft; FT(); } break;
    case  4: { FrequencyTable< 4> ft; FT(); } break;
    case  5: { FrequencyTable< 5> ft; FT(); } break;
    case  6: { FrequencyTable< 6> ft; FT(); } break;
    case  7: { FrequencyTable< 7> ft; FT(); } break;
    case  8: { FrequencyTable< 8> ft; FT(); } break;
    case  9: { FrequencyTable< 9> ft; FT(); } break;
    case 10: { FrequencyTable<10> ft; FT(); } break;
    default: MYASSERT_NEVERREACH();
  }
  #undef FT
  return 0;
}

