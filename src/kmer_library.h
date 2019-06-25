#ifndef _KMER_LIBRARY_HEADER
#define _KMER_LIBRARY_HEADER

#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cctype>
#include <map>
#include "stackdump.h"
#include "cpas_debug.h"
#include "cpas_tsv.h"


typedef unsigned long ulong;
typedef unsigned int uint;
typedef unsigned char Base;
typedef std::vector<Base> BString; ///< string of bases

const uint NUM_CHARS_FOR_BASES = 5;

constexpr ulong ipow(ulong i, uint N)
{
  return N == 0 ? 1 : N <= 1 ? i : i * ipow(i, N - 1);
}

/*
ulong ipow(ulong i, uint N){
  ulong retval = 1;
  if(N <= 0) return retval;
  if(N >  0){
    retval = ipow(i, N - 1);
  }
}
*/
const Base GAP_BASE = 4;
inline Base char2Base(char inch)
{
  switch(inch) {
    case 'A':
    case 'a':
      return 0;
    case 'C':
    case 'c':
      return 1;
    case 'G':
    case 'g':
      return 2;
    case 'T':
    case 't':
      return 3;
    case '-':
      return GAP_BASE;
    default:
      MYASSERT_NEVERREACH_WD(DUMP(int(inch)));
      exit(2);
  }
}

inline char base2Char(Base b)
{
  MYASSERT_WMD("Base must be 0 to 4", 0 <= b && b < 5, DUMP(b));
  return "ACGT-"[b];
}

inline std::string BString2String(const BString& b)
{
  std::string retval;
  retval.resize(b.size());
  for(uint i = 0; i < b.size(); i++) {
    retval[i] = base2Char(b[i]);
  }
  return retval;
}

inline BString String2BString(const std::string& s)
{
  BString retval;
  retval.resize(s.size());
  for(uint i = 0; i < s.size(); i++) {
    retval[i] = char2Base(s[i]);
  }
  return retval;
}



std::string splitBy1stSpace(const std::string s) {
  std::string retval;
  int len = 0;
  for (int i = 0; i < s.size(); i++) {
    if (s[i] != ' ') {
        len++;
      }else{
        break;
      }
  }
  retval.resize(len);
  for (int i = 0; i < len; i++){
    retval[i] = s[i];
  }
  retval[len] = '\0';
  return retval;
}




typedef std::string SequenceName;
typedef std::map<SequenceName, BString> MultiFASTA;

inline MultiFASTA loadFromFASTA(const std::string& inputFASTAFileName)
{
  MultiFASTA retval;
  std::ifstream ifs(inputFASTAFileName.c_str());
  if(ifs.fail()) {
    std::cerr << "ERROR: Cannot open file '" << inputFASTAFileName << "'\n";
    exit(2);
  }
  std::string tmp;
  std::string currentSequenceName;
  BString currentBString;
  while(std::getline(ifs, tmp)) {
    if(tmp.empty()) continue;
    if(tmp[0] == '>') {
      if(!currentSequenceName.empty()) {
        retval[currentSequenceName] = currentBString;
      }
      currentSequenceName = splitBy1stSpace(tmp.substr(1));
      currentBString.resize(0);
    } else {
      for(uint i = 0; i < tmp.size(); i++) {
        currentBString.push_back(char2Base(tmp[i]));
      }
    }
  }
  if(!currentSequenceName.empty()) {
    retval[currentSequenceName] = currentBString;
  }
  return retval;
}

/**
 * k-mer integer class
 */
template<uint K>
class KInt {
  uint64_t kint;
  static const size_t KMERINT_COUNT = ipow(NUM_CHARS_FOR_BASES, K);

public:
  /// Put a new base at LSB, and lift the other stuff toward MSB.
  /// The overflown digit is thrown away.
  inline void ShiftIn(Base b) {
    kint = (kint * NUM_CHARS_FOR_BASES + b) % KMERINT_COUNT;
  }
  inline void unshift(Base b) {
    kint = kint / NUM_CHARS_FOR_BASES + b * KMERINT_COUNT / NUM_CHARS_FOR_BASES;
  }
  inline bool hasgap() {
    bool flag = false;
    uint64_t tmp = kint;
    for(int a = 1; a <= KMERINT_COUNT; a++){
      if(tmp % 5 == 4){
        flag = true;
        return flag;
      }else{
        tmp /= 5;
      }
    }
    return flag;
  }
  inline KInt() : kint(0) {}
  inline KInt(uint64_t v) : kint(v) {}
  /// Construct by string
  /// e.g.) KInt k("ACG-");
  inline KInt(const char* s) : kint(0) {
    MYASSERT_WMD("KInt must be constructed with a string of size K", std::strlen(s) == K, DUMP(s));
    for(ulong i = 0; i < K; i++) {
      ShiftIn(char2Base(s[i]));
    }
  }

  /// Output string for debug
  inline std::string str() const {
    std::string retval;
    uint64_t v = kint;
    for(ulong i = 0; i < K; i++) {
      retval += base2Char(v % NUM_CHARS_FOR_BASES);
      v /= NUM_CHARS_FOR_BASES;
    }
    for(size_t i = 0; i < K / 2; i++) {
      std::swap(retval[i], retval[K - 1 - i]);
    }
    return retval;
  }
  inline operator uint64_t () const {
    return kint;
  }
};

/*
inline std::string KInt2String(const int input_kint, const int k)
{
  std::string retval;
  BString tmp;
  size_t cnt = 0;
  int i = input_kint;
  for(cnt = k - 1; cnt > 0; cnt--){
    i /= ipow(5, cnt);
    tmp.ShiftIn(i);
  }
  tmp.ShiftIn(i % 5);
  retval = BString2String(tmp);
  return retval;
}
*/

struct CIGAROp {
  char op;
  int len;
  inline CIGAROp() {}
  inline CIGAROp(char op, int len) : op(op), len(len) {}
};

typedef std::vector<CIGAROp> CIGAROPS;

inline std::string CIGAROps2String(const CIGAROPS& cops) {
  std::string retval;
  char buf[8];
  for(auto& x: cops) {
    std::sprintf(buf, "%d%c", x.len, x.op);
    retval += buf;
  }
  return retval;
}

inline CIGAROPS parseCIGARString(const std::string& cigarString) {
  CIGAROPS retval;
  const std::string& s = cigarString;
  for(size_t i = 0; i < s.size(); ++i) {
    // Parse digits
    int num = 0;
    while(std::isdigit(s[i])) {
      num = num * 10 + s[i] - '0';
      i++;
    }
    if(num <= 0) {
      std::cerr << "CIGAR op size must be positive (num = " << num << ", cigar = '" << s << "'" << std::endl;
      std::exit(2);
    }
    char op = s[i];
    //std::cerr << num << op << std::endl;
    switch(op) {
    case 'I':
    case 'D':
    case 'M':
    case 'X':
    case '=':
    case 'H':
    case 'S':
    case 'P':
    case 'N':
      break;
    default:
      std::cerr << "CIGAR type '" << op << "' is not valid." << std::endl;
      exit(2);
    }
    retval.push_back(CIGAROp(op, num));
  }
  //std::cerr << std::endl;
  return retval;
}

inline void generateAlignmentSequencesFromCIGARAndSeqs(
  const BString& refbs,
  const BString& querybs,
  const CIGAROPS& ops,
  int refStartPos,
  int queryStartPos,
  BString& ras,
  BString& qas
) {
  size_t alignmentSize = 0;
  for(auto& x: ops) {
    switch(x.op) {
    case 'M':
    case 'I':
    case 'D':
    case 'X':
    case '=':
    case 'P':
      alignmentSize += x.len;
    }
  }
  ras.resize(alignmentSize);
  qas.resize(alignmentSize);
  size_t ai = 0;
  size_t rp = refStartPos;
  size_t qp = queryStartPos;
  for(auto& x: ops) {
    switch(x.op) {
    case 'M':
    case '=':
    case 'X':
      for(int i = 0; i < x.len; ++i) {
        ras[ai] = refbs[rp++];
        qas[ai] = querybs[qp++];
        ai++;
      }
      break;
    case 'I':
      for(int i = 0; i < x.len; ++i) {
        ras[ai] = GAP_BASE;
        qas[ai] = querybs[qp++];
        ai++;
      }
      break;
    case 'D':
      for(int i = 0; i < x.len; ++i) {
        ras[ai] = refbs[rp++];
        qas[ai] = GAP_BASE;
        ai++;
      }
      break;
    case 'S':
      qp += x.len;
      break;
    case 'H':
      break;
    case 'P':
    case 'N':
      // not implemented
      std::cout << "ERROR: P/N in CIGAR are not supported." << std::endl;
      std::exit(2);
    default:
      std::cout << "ERROR: Logic error." << std::endl;
      std::exit(2);
    }
  }
}


struct SAMRecord {
  std::string qname;
  int flag;
  std::string rname;
  int pos;   ///< 0-origin position
  int mapQ;
  std::string cigar;
  std::string rnext;
  int pnext;
  int tlen;
  std::string seq;
  std::string qual;

  /// fill SAMRecord from the current line
  /// @return true when the parsing succeeds
  bool fill(FastTSVParse& ftp) {
    // See: https://samtools.github.io/hts-specs/SAMv1.pdf for the detailed specification for SAM.
    // Col Field Type
    // 1 QNAME String
    // 2 FLAG Int
    // 3 RNAME String
    // 4 POS Int
    // 5 MAPQ Int
    // 6 CIGAR String
    // 7 RNEXT String
    // 8 PNEXT Int
    // 9 TLEN Int
    // 10 SEQ String
    // 11 QUAL String
    if(ftp.size() < 11) return false;
    qname = ftp.getString(0);
    flag  = ftp.getInteger(1);
    rname = ftp.getString(2);
    pos   = ftp.getInteger(3) - 1;
    mapQ  = ftp.getInteger(4);
    cigar = ftp.getString(5);
    rnext = ftp.getString(6);
    pnext = ftp.getInteger(7);
    tlen  = ftp.getInteger(8);
    seq   = ftp.getString(9);
    qual  = ftp.getString(10);
    return true;
  }
};


inline Base Base2CompBase(Base inb)
{
  switch(inb) {
    case 0:     //when got 'A'
      return 3; //return   'T'

    case 1:     //when got 'C'
      return 2; //retern   'G'

    case 2:     //when got 'G'
      return 1; //return   'C'

    case 3:     //when got 'T'
      return 0; //return   'A'
    case 4:
      return GAP_BASE;
    default:
      MYASSERT_NEVERREACH_WD(DUMP(int(inb)));
  }
}


inline void revCompBString(BString &b)
{
  BString tmp;
  size_t len = b.size();
  tmp.resize(len);
  for(size_t idx = 0; idx < len; idx++){
    tmp[idx] = Base2CompBase(b[len - idx - 1]);
  }
  for(size_t idx = 0; idx < len; idx++){
    b[idx] = tmp[idx];
  }
}



inline double* lacalNormalization(const int* mtx)
{

#define ind(i, j) ( (j * 5 + i) )
  double *retval = (double*)malloc(sizeof(double) * 25);
  uint64_t sum_all    = 0;
  for(int i = 0; i < 25; i++){
    retval[i] = 0;
  }
  for(int i = 0; i < 25; i++){
    sum_all += mtx[i];
  }

  if(sum_all == 0){ // need to be well concidered.
    for(int i = 0; i < 24; i++){ retval[i] = 1.0 / 24;}
    return retval;
  }
  uint64_t sum_ins = 0;
  for(int j = 0; j < 4; j++){
    sum_ins += mtx[ind(4, j)];
  }
  for(int j = 0; j < 4; j++){
    retval[ind(4, j)] = (double)mtx[ind(4, j)] / sum_all;
  }


  if(sum_all == sum_ins){//in case that only insertion line had values.
    for(int j = 0; j < 4; j++){
      retval[ind(4, j)] = mtx[ind(4, j)] / sum_ins;
    }
    return retval;
  }

  double not_ins_ratio = (double)(sum_all - sum_ins) / (sum_all * 4);

  for(int i = 0; i < 4; i++){
    uint64_t sum_line = 0;
    for(int j = 0; j < 5; j++){
      sum_line += mtx[ind(i, j)];
    }
    for(int j = 0; j < 5; j++){
      if(sum_line != 0){
        retval[ind(i, j)] = (double)mtx[ind(i, j)] * not_ins_ratio / sum_line;
      }else{
        retval[ind(i, j)] = 0;
      }
    }
  }

#undef index_25_mtx
  retval[24] = 0;
  return retval;
}


#endif // #ifndef _KMER_LIBRARY_HEADER
