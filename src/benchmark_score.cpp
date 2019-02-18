#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <getopt.h>
#include <math.h>
#include "../../ms_thesis_scripts/kmer_library.h"

#define TESTCOUNT 1000

//"score_benchmark.rb <test_data.fa> <int kmer_size> <int seq_length(<=30)> <binary_scoretable_path>\n"
using namespace std;

int main(int argc, char* argv[]){
  const string fastaFileName = argv[1];
  const int kmerSize         = atoi(argv[2]);
  const int seqLen           = atoi(argv[3]) - 1;
  string binaryPath = "";
  if(argc > 4){
    binaryPath = argv[4];
  }
  FILE* fastaFp;
/*  fastaFp = fopen(fastaFileName , "r" );
  if(fastaFp == NULL){
    fprintf(stderr, "could not open %s\n", fastaFileName);
    exit(2);
  }
  */
  char refseq[32];
  char qryseq[32];
  string tmp;
  int lineCount = 0;
  vector <string> reads;
  ifstream ifs(fastaFileName.c_str());
  while(getline(ifs, tmp)){
      if(tmp.empty()) continue;
      if(tmp[0] == '>') continue;
      lineCount++;
      strcpy(qryseq, tmp.c_str());
      reads.push_back(qryseq);
  }
  for(int i = 0; i < reads.size(); i+=2){
    const string head = "~/ms_thesis/minimi/src/NWalignment ";
    const string spacer = " ";
    const string cmd = head + reads[i].substr(0, seqLen).c_str() + spacer + reads[i + 1].substr(0, seqLen).c_str() + spacer + to_string(kmerSize) + spacer + binaryPath;
    printf("%s\n", cmd.c_str());
    //system(cmd);
  }


  //fclose(fastaFp);
  return 0;
}

