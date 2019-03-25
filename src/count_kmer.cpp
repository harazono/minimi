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
  fprintf(stderr, "\n");
  exit(2);
}

template<uint KMERSIZE>
struct FrequencyTable {
  static const size_t context_tablesize = ipow(5, 2 * (KMERSIZE - 1));
  static const size_t next_tablesize    = ipow(5, 2);
  typedef int Frequency;
  typedef int Score;
  vector<Frequency> context_table; ///< context table.
  vector<Frequency> kmer_kmer_table; ///< The first dimension is reference, the second is query.
  vector<double>    kmer_kmer_prob_table;
  vector<Score>     score_table;
  inline Frequency& c_table(size_t c) {
    MYASSERT_WMD("Out of range (c)", c < context_tablesize, DUMP(c));
    return context_table[c];
  }
  inline Frequency& kk(size_t c, size_t n) {
    MYASSERT_WMD("Out of range (c)", c < context_tablesize, DUMP(c));
    MYASSERT_WMD("Out of range (n)", n < next_tablesize, DUMP(n));
    return kmer_kmer_table[c * next_tablesize + n];
  }
  inline double& kkp(size_t c, size_t n) {
    MYASSERT_WMD("Out of range (c)", c < context_tablesize, DUMP(c));
    MYASSERT_WMD("Out of range (n)", n < next_tablesize, DUMP(n));
    return kmer_kmer_prob_table[c * next_tablesize + n];
  }
  inline Score& scr(size_t c, size_t n) {
    MYASSERT_WMD("Out of range (c)", c < context_tablesize, DUMP(c));
    MYASSERT_WMD("Out of range (n)", n < next_tablesize, DUMP(n));
    return score_table[c * next_tablesize + n];
  }

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
    size_t var_table_size = context_tablesize * next_tablesize;
    fwrite(&var_table_size, sizeof(var_table_size), 1, fle);
    //MYASSERT_WMD("var_table_size == kmer_kmer_table.size()", var_table_size == kmer_kmer_table.size(), DUMP(var_table_size, tablesize, kmer_kmer_table.size()));
    fwrite(&*score_table.begin(), sizeof(Frequency) * score_table.size(), 1, fle);
    fclose(fle);
  }
  void countAlignment(
    const BString& ras,
    const BString& qas
  ){
    MYASSERT_WMD("ras.size() must be qas.size()", ras.size() == qas.size(), DUMP(ras.size(), qas.size()));
    const size_t aligned_len = ras.size();
    if(aligned_len < KMERSIZE) return;
    KInt<2 * (KMERSIZE - 1)> cki;// context kmer index
    KInt<2>                  nki;// next    kmer index

    if(KMERSIZE == 1){// there is no context.
      //cki.ShiftIn(char2Base('A'));
      for(size_t i = KMERSIZE - 1; i <= aligned_len - KMERSIZE; i++){
        nki.ShiftIn(ras[i]);
        nki.ShiftIn(qas[i]);
        c_table(cki)++;
        kk(cki, nki)++;
      }

    }else{

      if(KMERSIZE > 2){
        for(size_t i = 0; i < KMERSIZE - 2; i++){
          cki.ShiftIn(ras[i]);
          cki.ShiftIn(qas[i]);
        }
      }
      for(size_t i = KMERSIZE - 1; i <= aligned_len - KMERSIZE; i++){
        cki.ShiftIn(ras[i - 1]);
        cki.ShiftIn(qas[i - 1]);
        nki.ShiftIn(ras[i]);
        nki.ShiftIn(qas[i]);
        c_table(cki)++;
        kk(cki, nki)++;
      }
    }
  }

  void printCTable(){
    for(int i = 0; i < context_tablesize; i++){
      fprintf(stdout, "%10d, ", c_table(i));
    }
    fprintf(stdout, "\n");
  }

  void printKKTable(){
    for(int j = 0; j < next_tablesize; j++){
      for(int i = 0; i < context_tablesize; i++){
        fprintf(stdout, "%10d", kk(i, j));
        if(i != context_tablesize) fprintf(stdout, ", ");
      }
      fprintf(stdout, "\n");
    }
  }

  void printProbTable(){
    for(int j = 0; j < next_tablesize; j++){
      for(int i = 0; i < context_tablesize; i++){
        fprintf(stdout, "%10f", kkp(i, j));
        if(i != context_tablesize) fprintf(stdout, ", ");
      }
      fprintf(stdout, "\n");
    }
  }
  void printScoreTable(){
    for(int j = 0; j < next_tablesize; j++){
      for(int i = 0; i < context_tablesize; i++){
        fprintf(stdout, "%10d", scr(i, j));
        if(i != context_tablesize) fprintf(stdout, ", ");
      }
      fprintf(stdout, "\n");
    }
  }

  size_t score_count = 1;// for printing progress report
  void scorerize(const int diff){
    /*vector<double> kmer_kmer_prob_table(context_tablesize * next_tablesize, 0);
    auto kkp = [&kmer_kmer_prob_table](size_t r, size_t q) -> double& {
      MYASSERT_WMD("Out of range (r)", r < context_tablesize, DUMP(r));
      MYASSERT_WMD("Out of range (q)", q < next_tablesize, DUMP(q));
      return kmer_kmer_prob_table[r * context_tablesize + q];
    };*/

    #pragma omp parallel num_threads(16)
    {
      #pragma omp for
      for(int c = 0; c < context_tablesize; c++){
        for(int n = 0; n < next_tablesize; n++){
          if(score_count % 10000 == 0) fprintf(stderr, "first roop : %'d / %'d\r", score_count, context_tablesize * next_tablesize);
          #pragma omp atomic
          score_count++;
          if(c_table(c) != 0){
            kkp(c, n) = static_cast<double>(kk(c, n)) / c_table(c);//divide by context k-mer frequency.
          }else{
            kkp(c, n) = 0.0;
          }
          //fprintf(stderr, "kk(c, n) = %d, c_table(c) = %d, kkp(c, n) = %f\n", kk(c, n), c_table(c), kkp(c, n));
          MYASSERT_WMD("prob must be in [0, 1]", kkp(c, n) <= 1.0 && kkp(c, n) >= 0.0, DUMP(kkp(c, n)));
        }
      }
      fprintf(stderr, "first roop : %'d / %'d                         \r", score_count, context_tablesize * next_tablesize);
      score_count = 1;
    }
    score_count = 0;
    for(int c = 0; c < context_tablesize; c++){
      for(int n = 0; n < next_tablesize; n++){
        if(score_count % 10000 == 0) fprintf(stderr, "second  roop : %'d / %'d                                               \r", score_count, context_tablesize * next_tablesize) ;
        #pragma omp atomic
        score_count++;
        if(kkp(c, n) > 0.0){
          scr(c, n) = diff + static_cast<int>(round(100 * log10(kkp(c, n))));
        }else{
          scr(c, n) = diff + -1 * ipow(2, 10);
        }
      }
    }
    fprintf(stderr, "Done                                                                             \r");
  }


public:
  FrequencyTable() : context_table(context_tablesize, 0), kmer_kmer_table(context_tablesize * next_tablesize, 0), kmer_kmer_prob_table(context_tablesize * next_tablesize, 0), score_table(context_tablesize * next_tablesize, 0) {}
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

    fprintf(stderr, "Loading from the FASTA file ...\r");
    const MultiFASTA multiFASTA = loadFromFASTA(FASTAFileName);
    fprintf(stderr, "Done                           \n");
    FastTSVParse ftp(SAMFileName);
    if(!ftp) {
      fprintf(stderr, "ERROR: Cannot open SAM file '%s'\n", SAMFileName);
      exit(2);
    }
/*
    fprintf(stderr, "context_tablesize = %d, next_tablesize = %d\n", context_tablesize, next_tablesize);
    fprintf(stderr, "kmer_table.size()           = %d\n", kmer_table.size());
    fprintf(stderr, "kmer_kmer_table.size()      = %d\n", kmer_kmer_table.size());
    fprintf(stderr, "kmer_kmer_prob_table.size() = %d\n", kmer_kmer_prob_table.size());
    */

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
        exit(2);
      }

      // get aligned sequence at here
      //cerr << record.cigar.c_str() << endl;
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
    cerr << "Done." << endl;

    scorerize(100);
    if(!binaryOutputFileName.empty()) {
      outputAsBinaryTable(binaryOutputFileName);
    }
    if(outputInCSV) {
      printCTable();
      fprintf(stdout, "\n");
      printKKTable();
      fprintf(stdout, "\n");
      printProbTable();
      fprintf(stdout, "\n");
      printScoreTable();
      fprintf(stdout, "\n");
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
      if(kmer_size < 1 || 6 < kmer_size) {
        error("kmer size must be 1-6");
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
    case 1: { FrequencyTable<1> ft; FT(); } break;
    case 2: { FrequencyTable<2> ft; FT(); } break;
    case 3: { FrequencyTable<3> ft; FT(); } break;
    case 4: { FrequencyTable<4> ft; FT(); } break;
    case 5: { FrequencyTable<5> ft; FT(); } break;
    case 6: { FrequencyTable<6> ft; FT(); } break;
    default: MYASSERT_NEVERREACH();
  }
  #undef FT
  return 0;
}

