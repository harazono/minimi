#include "gtest/gtest.h"
#include <string>
#include "kmer_library.h"

using namespace std;

class MyLibraryTest : public ::testing::Test {
  protected:
        // データメンバーの初期化
  virtual void SetUp() {
  }
};


TEST_F(MyLibraryTest, KInt_Simple1_Test) {
  KInt<1> k1;
  EXPECT_EQ(k1, 0);
  EXPECT_STREQ(k1.str().c_str(), "A");
  KInt<1> k2("A");
  EXPECT_EQ(k2, 0);
  KInt<1> k3("C");
  EXPECT_EQ(k3, 1);
  KInt<1> k4("G");
  EXPECT_EQ(k4, 2);
  KInt<1> k5("T");
  EXPECT_EQ(k5, 3);
  KInt<1> k6("-");
  EXPECT_EQ(k6, 4);
}

TEST_F(MyLibraryTest, KInt_Simple2_Test) {
  KInt<2> k1;
  EXPECT_EQ(k1, 0);
  KInt<2> k2("AA");
  EXPECT_EQ(k2, 0);
  KInt<2> k3("AC");
  EXPECT_EQ(k3, 1);
  KInt<2> k4("CA");
  EXPECT_EQ(k4, 5);
  KInt<2> k5("CC");
  EXPECT_EQ(k5, 6);
  KInt<2> k6("--");
  EXPECT_EQ(k6, 24);
}

TEST_F(MyLibraryTest, KInt_Simple3_Test) {
  KInt<3> k1;
  EXPECT_EQ(k1, 0);
  KInt<3> k4("GGC");
  EXPECT_EQ(k4, 61);
  KInt<3> k5("G-G");
  EXPECT_EQ(k5, 72);
  KInt<3> k6("-AG");
  EXPECT_EQ(k6, 102);
}

TEST_F(MyLibraryTest, KInt_ShiftIn_Test) {
  int gap = 4;
  KInt<3> k4("GGC");
  KInt<3> k5("GC-");
  KInt<3> k6("C--");
  k4.ShiftIn(gap);
  EXPECT_EQ(k4, k5);
}

TEST_F(MyLibraryTest, KInt_unshift_Test) {
  int gap = 4;
  KInt<3> k4("GGC");
  KInt<3> k5("-GG");
  KInt<3> k6("--G");
  KInt<3> k7("---");
  k4.unshift(gap);
  EXPECT_EQ(k4, k5);
  k4.unshift(gap);
  EXPECT_EQ(k4, k6);
  k4.unshift(gap);
  EXPECT_EQ(k4, k7);
}

TEST_F(MyLibraryTest, KInt_hasgap_Test) {
  int gap = 4;
  KInt<3> k4("GGG");
  KInt<3> k5("-GG");
  KInt<3> k6("--G");
  KInt<3> k7("---");
  KInt<3> k8("G--");
  KInt<3> k9("GG-");
  EXPECT_EQ(k4.hasgap(), false);
  EXPECT_EQ(k5.hasgap(), true);
  EXPECT_EQ(k6.hasgap(), true);
  EXPECT_EQ(k7.hasgap(), true);
  EXPECT_EQ(k8.hasgap(), true);
  EXPECT_EQ(k9.hasgap(), true);


}

// 2019/06/10
// need to add a test case in which multi fasta file should be tested.
TEST_F(MyLibraryTest, LoadMultiFASTA_Test) {
  MultiFASTA mf = loadFromFASTA("10.ref");
  EXPECT_EQ(mf.size(), 4);
  EXPECT_STREQ(BString2String(mf.at("apple")).c_str(), "CGTCTATTCC");
  EXPECT_STREQ(BString2String(mf.at("banana")).c_str(), "CGTCTATTCCCGTCTATTCC");
  EXPECT_STREQ(BString2String(mf.at("citron")).c_str(), "ATTTAGAAGGAGGTAGATAA");
  EXPECT_STREQ(BString2String(mf.at("dragon")).c_str(), "AAAAAAAAAAAAAA");
  MultiFASTA mf2 = loadFromFASTA("/home/harazono/r10_chemistry/ZymoBIOMICS.STD.refseq.v2/test.fa");
  //EXPECT_EQ(mf2.size(), 4);
  /*
  EXPECT_STREQ(mf.begin()->first.c_str(), "apple");
  EXPECT_STREQ(BString2String(mf.begin()->second).c_str(), "CGTCTATTCC");
  EXPECT_STREQ(mf.begin()->first.c_str(), "banana");
  EXPECT_STREQ(BString2String(mf.begin()->second).c_str(), "CGTCTATTCC");
  */
  }




TEST_F(MyLibraryTest, SplitBy1stSpace_Test) {
  std::string head = splitBy1stSpace("AAA WWW");
  EXPECT_STREQ(head.c_str(), "AAA");
  EXPECT_EQ(head.size(), 3);
  head = splitBy1stSpace("i;oarehwv wir efw;ioh");
  EXPECT_STREQ(head.c_str(), "i;oarehwv");
  EXPECT_EQ(head.size(), 9);
}

TEST_F(MyLibraryTest, ParseSAM_Test) {
  FastTSVParse ftp("10.sam");
  ASSERT_FALSE(!ftp);
  while(true) {
    ASSERT_TRUE(ftp.readNextLine());
    if(ftp.c_str()[0] != '@') break;
  }
  SAMRecord r;
  //cout << ftp.c_str() << endl;
  //cout << ftp.size() << endl;
  ASSERT_TRUE(r.fill(ftp));
  EXPECT_STREQ(r.qname.c_str(), "testread");
  EXPECT_STREQ(r.rname.c_str(), "testref");
  EXPECT_EQ(r.flag, 0);
  EXPECT_STREQ(r.cigar.c_str(), "10M");
  EXPECT_STREQ(r.seq.c_str(), "CGTCTATTCC");
  EXPECT_STREQ(r.qual.c_str(), "*");
  EXPECT_STREQ(r.rnext.c_str(), "*");
  EXPECT_EQ(r.tlen, 0);
  EXPECT_EQ(r.pnext, 0);
  EXPECT_EQ(r.mapQ, 60);
  EXPECT_EQ(r.pos, 0);
}

TEST_F(MyLibraryTest, BString_Test) {
  const string s = "ACCGGGTGGTC--G-";
  const BString bs = String2BString(s);
  const string rev_s = BString2String(bs);
  EXPECT_STREQ(s.c_str(), rev_s.c_str());
  EXPECT_EQ(bs[0], 0);
  EXPECT_EQ(bs[1], 1);
  EXPECT_EQ(bs.back(), 4);
}

TEST_F(MyLibraryTest, ParseCIGAR_Simple1_Test) {
  CIGAROPS ops = parseCIGARString("1M");
  ASSERT_EQ(ops.size(), 1);
  EXPECT_EQ(ops[0].op, 'M');
  EXPECT_EQ(ops[0].len, 1);
}

TEST_F(MyLibraryTest, ParseCIGAR_Simple2_Test) {
  CIGAROPS ops = parseCIGARString("4S21M412D");
  ASSERT_EQ(ops.size(), 3);
  EXPECT_EQ(ops[0].op, 'S');
  EXPECT_EQ(ops[0].len, 4);
  EXPECT_EQ(ops[1].op, 'M');
  EXPECT_EQ(ops[1].len, 21);
  EXPECT_EQ(ops[2].op, 'D');
  EXPECT_EQ(ops[2].len, 412);
}

TEST_F(MyLibraryTest, GenerateAlignmentFromCIGAR_Test) {
  const string refSeq     = "ACGTGCGT";
  const string querySeq   = "CGGCAG";
  const BString refBSeq   = String2BString(refSeq);
  const BString queryBSeq = String2BString(querySeq);
  // ACGTGC-GT
  //  || || |
  //  CG-GCAG
  CIGAROPS ops = parseCIGARString("1H2M1D2M1I1M1H");
  BString ras, qas;
  generateAlignmentSequencesFromCIGARAndSeqs(refBSeq, queryBSeq, ops, 1, 0, ras, qas);
  const string rass = BString2String(ras);
  const string qass = BString2String(qas);
  EXPECT_STREQ(rass.c_str(), "CGTGC-G");
  EXPECT_STREQ(qass.c_str(), "CG-GCAG");
}

TEST_F(MyLibraryTest,revComp1_Test){
  string  preSeq  = "AAAC";
  BString preBSeq = String2BString(preSeq);
  revCompBString(preBSeq);
  string ansSeq   = BString2String(preBSeq);
  EXPECT_STREQ(ansSeq.c_str(), "GTTT");
}

TEST_F(MyLibraryTest,revComp2_Test){
  string  preSeq  = "A";
  BString preBSeq = String2BString(preSeq);
  revCompBString(preBSeq);
  string ansSeq   = BString2String(preBSeq);
  EXPECT_STREQ(ansSeq.c_str(), "T");
}

TEST_F(MyLibraryTest,revComp3_Test){
  string  preSeq  = "AAACG";
  BString preBSeq = String2BString(preSeq);
  revCompBString(preBSeq);
  string ansSeq   = BString2String(preBSeq);
  EXPECT_STREQ(ansSeq.c_str(), "CGTTT");
}

TEST_F(MyLibraryTest,revComp4_Test){
  string  preSeq  = "A-ACG";
  BString preBSeq = String2BString(preSeq);
  revCompBString(preBSeq);
  string ansSeq   = BString2String(preBSeq);
  EXPECT_STREQ(ansSeq.c_str(), "CGT-T");
}


TEST_F(MyLibraryTest,localNormalization_Test1){
  int *localTable = (int*)malloc(sizeof(int) * 25);
  for(int i = 0; i < 25; i++){
    localTable[i] = 10000;
  }
  localTable[24] = 0;
  double *localprob = (double*)malloc(sizeof(double) * 25);
  localprob = lacalNormalization(localTable);
  double sum = 0;
  for(int i = 0; i < 25; i++){
    sum += localprob[i];
  }
  ASSERT_TRUE(sum - 1.0 < 0.001 && 1.0 - sum < 0.001);
}


TEST_F(MyLibraryTest,localNormalization_Test2){
  //int *localTable = (int*)malloc(sizeof(int) * 25);
  int localTable[] = {7671605, 1878030, 409216, 1946895, 8931149, 7685734, 6046752, 4981902, 5911312, 6719582, 1143531, 6398762, 225281, 53472, 4588717, 2709910, 3803638, 6196264, 4410145, 6005817, 7383021, 7429599, 1907415, 4834045, 0};
  double *localprob = (double*)malloc(sizeof(double) * 25);
  localprob = lacalNormalization(localTable);
  double sum = 0;
  for(int i = 0; i < 25; i++){
    sum += localprob[i];
  }
  ASSERT_TRUE(sum - 1.0 < 0.001 && 1.0 - sum < 0.001);
}

TEST_F(MyLibraryTest,localNormalization_Test3){
  //int *localTable = (int*)malloc(sizeof(int) * 25);
  int localTable[] = {0, 11, 9, 2, 0, 0, 12, 31, 0, 0, 31, 45, 111, 0, 0, 65, 12, 32, 0, 0, 11, 11, 11, 11, 0};
  double *localprob = (double*)malloc(sizeof(double) * 25);
  localprob = lacalNormalization(localTable);
  double sum = 0;
  for(int i = 0; i < 25; i++){
    sum += localprob[i];
  }
  ASSERT_TRUE(sum - 1.0 < 0.001 && 1.0 - sum < 0.001);
}

TEST_F(MyLibraryTest,localNormalization_Test4){
  //int *localTable = (int*)malloc(sizeof(int) * 25);
  int localTable[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  double *localprob = (double*)malloc(sizeof(double) * 25);
  localprob = lacalNormalization(localTable);
  double sum = 0;
  for(int i = 0; i < 25; i++){
    sum += localprob[i];
  }
  ASSERT_TRUE(sum - 1.0 < 0.001 && 1.0 - sum < 0.001);
}

TEST_F(MyLibraryTest,localNormalization_Test5){
  //int *localTable = (int*)malloc(sizeof(int) * 25);
  int localTable[] = {1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0};
  double *localprob = (double*)malloc(sizeof(double) * 25);
  localprob = lacalNormalization(localTable);
  double sum = 0;
  for(int i = 0; i < 25; i++){
    sum += localprob[i];
  }
  ASSERT_TRUE(sum - 1.0 < 0.001 && 1.0 - sum < 0.001);
}
TEST_F(MyLibraryTest,localNormalization_Test6){
  //int *localTable = (int*)malloc(sizeof(int) * 25);
  int localTable[] = {125243265, 831388, 1594635, 990496, 2042198, 1311240, 127379226, 968296, 6075575, 2546654, 1830296, 601625, 128516423, 1074008, 3062296, 1009615, 5800435, 1248973, 117309637, 2979892, 8888301, 8244568, 10193576, 12716559, 0};
  double *localprob = (double*)malloc(sizeof(double) * 25);
  localprob = lacalNormalization(localTable);
  double sum = 0;
  for(int i = 0; i < 25; i++){
    sum += localprob[i];
  }
  ASSERT_TRUE(sum - 1.0 < 0.001 && 1.0 - sum < 0.001);
}
  int main(int argc, char **argv) {
      GDB_On_SEGV g(argv[0]);
        ::testing::InitGoogleTest(&argc, argv);
          return RUN_ALL_TESTS();
}

