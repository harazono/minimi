// vim:ff=unix ft=cpp ts=4 sw=4 sts=4 si et :
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <string>
#include <vector>
#include <limits.h>
#include "kmer_library.h"
#include "score_1mer.h"
#include "score_1mer_mod_3digit.h"
#include "score_2mer_mod_3digit.h"
#include "score_3mer_mod_3digit.h"
#include "score_4mer_mod_3digit.h"

#define _pos(i, j, ilen)	( (j) * ((ilen) + 1) + (i) )
#define MAX2(p, q)			( ((p) < (q)) ? (q) : (p) )
#define MAX3(p, q, r)		( MAX2(p, MAX2(q, r)) )
#define MAX4(p, q, r, s)	( MAX2(MAX2(p, q), MAX2(r, s)) )

#define ASSERT(cond, message) do { if(!(cond)){ fprintf(stderr, "ASSERTION ERROR: %s\nCONDITION: %s\n", message, #cond); exit(2); } } while(0)

using namespace std;


uint32_t shift_digit(uint32_t pre, uint32_t add, size_t kmer_size){
    //fprintf(stderr, "%d, %d, %d\n", pre, add, kmer_size);
    ASSERT(0 <= add && add < 5, "add must be 0-4");
    size_t num_kmers = 1;
    for(size_t i = 0; i < kmer_size; i++) {
        num_kmers *= 5;
    }
    ASSERT(0 <= pre && pre < num_kmers, "pre is out of the range");
    return pre * 5 % num_kmers + add;
}

uint32_t base2int(char base){
    //fprintf(stderr, "base2int got %d\n", base);
	if(base == 'A' || base == 'a') return 0;
	if(base == 'C' || base == 'c') return 1;
	if(base == 'G' || base == 'g') return 2;
	if(base == 'T' || base == 't') return 3;
	if(base == '-')                return 4;
    ASSERT(0, "Logic error in base2int");
}
// for debug. use someday
int int2base(uint32_t intBase) {
    ASSERT(0 <= intBase && intBase < 5, "out of range");
    return "ACGT-"[intBase];
}

struct Context {
    size_t con;
    Context(){
        con = 0;
    }
    void slide(char ref_char, char query_char, int kmer_size){
        con = shift_digit(con, base2int(ref_char),   2 * (kmer_size - 1));
        con = shift_digit(con, base2int(query_char), 2 * (kmer_size - 1));
    }
};

enum TraceBackDirection {
    TB_MATCH = 0, ///< Traceback to the upperleft direction (Match/Mismatch)
    TB_LEFT  = 1, ///< Trackback to the left  (insertion to the reference)
    TB_UP    = 2, ///< Traceback to the upper (deletion in the reference)
    TB_START = 3  ///< This indicates the beginning of the trace path
};

struct Matrix {
	int32_t score;
	TraceBackDirection tb_dir;   ///< Traceback direction.(32bit)
	int64_t ipos;
	int64_t jpos;
	Context context;
    Matrix() : context() {
        score = ipos = jpos = 0;
        tb_dir = TB_MATCH;
    }
    Matrix(int16_t score, int64_t ipos, int64_t jpos)
        : score(score), ipos(ipos), jpos(jpos), context() {}
};


/***
 *
 *                       reference
 *              0   1   2   ...   i   ...   a
 *
 *          0
 *
 *          1
 *
 *          2
 *          .
 *          .
 *          .
 *   query  j
 *          .
 *          .
 *          .
 *
 *
 *          b
 *
 * S(i-1,j-1)  |  S(i, j-1)
 *                         v INS
 * S(i, j-1)   |  S(i, j)
 *            >
 *           DEL
 *
 * S(i, j) = max(0,
 *               S(i - 1, j - 1) + score(context, r[i], q[j]),
 *               S(i - 1, j    ) + score(context, '-', q[j] ), DEL
 *               S(i    , j - 1) + score(context, r[i], '-' )  INS
 *               )
 * then, update context(i, j)
 *
 */

int str2idx(char *str){
	size_t i = 0;
	int sum  = 0;
	size_t len = strlen(str);
	while(str[i] != '\0'){
		sum += base2int(str[i]) * ipow(5, len - i - 1);
		i++;
	}
	return sum;
}
string idx2str(int idx, size_t kmer_size) {
    string tmp;
    for(size_t i = 0; i < kmer_size; i++) {
        tmp += int2base(idx % 5); idx /= 5;
    }
    string n;
    for(size_t i = 0; i < tmp.size(); i++) {
        n += tmp[tmp.length() - 1 - i];
    }
    return n;
}
/*
void showMatrix(const vector<vector<Matrix> >& mtx)
{
    fprintf(stderr, "======\n");
    fprintf(stderr, "score\n");
    for(size_t y = 0; y < mtx.size(); y++) {
        const vector<Matrix>& line = mtx[y];
        for(size_t x = 0; x < line.size(); x++) {
            fprintf(stderr, "%5d ", line[x].score);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "context.ref_seq\n");
    for(size_t y = 0; y < mtx.size(); y++) {
        const vector<Matrix>& line = mtx[y];
        for(size_t x = 0; x < line.size(); x++) {
            fprintf(stderr, "%5zu ", line[x].context.ref_seq);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "context.query_seq\n");
    for(size_t y = 0; y < mtx.size(); y++) {
        const vector<Matrix>& line = mtx[y];
        for(size_t x = 0; x < line.size(); x++) {
            fprintf(stderr, "%5zu ", line[x].context.query_seq);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "tb_dir\n");
    for(size_t y = 0; y < mtx.size(); y++) {
        const vector<Matrix>& line = mtx[y];
        for(size_t x = 0; x < line.size(); x++) {
            fprintf(stderr, "%5d ", line[x].tb_dir);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "======\n");
}
*/

Matrix needlman_wunsch(
        const char* const reference_sequence,
        const char* const query_sequence,
        const size_t kmer_size,
        const int32_t kmer_score_table[]
) {
    const char* const rs = reference_sequence;
    const char* const qs = query_sequence;
	const uint32_t reference_length = strlen(reference_sequence);
	const uint32_t query_length = strlen(query_sequence);
    const uint32_t rl  = reference_length;
    const uint32_t ql  = query_length;
    //fprintf(stderr, "finish allocation working memory\n");
    int all_gap_index = ipow(5, 2 * (kmer_size - 1)) - 1;
    if(kmer_size == 1){
        all_gap_index = base2int('-');
    }
    /*{
        char gaps[kmer_size + 1];
        for(size_t i = 0; i < 2 * (kmer_size - 1); i++){
            gaps[i] = '-';
        }
        gaps[kmer_size] = '\0';
        all_gap_index = str2idx(gaps);
    }*/

    const size_t context_idx_size = ipow(5, 2 * (kmer_size - 1));
    #define table_index(i, j) ( (i) * 25 + (j) )
    #define ti(i, j) table_index(i, j)
    const int32_t *st = kmer_score_table;

    vector<vector<Matrix> > mtx(rl + 1, vector<Matrix>(ql + 1));
    // initialize
    //fill context, score and tb_dir

    for(size_t rsi = 0; rsi <=rl; rsi++){
        for(size_t qsi = 0; qsi <= ql; qsi++){
            mtx[rsi][qsi].context.con = all_gap_index;
        }
    }

    Context home_context;
    Context hc   = home_context;
    hc.con       = all_gap_index;
    hc.slide(rs[0], qs[0], kmer_size);
    mtx[0][0].context = hc;
    mtx[0][0].score   = 0;
    mtx[0][0].tb_dir  = TB_START;

    for(size_t rsi = 1; rsi < rl; rsi++){
        Matrix pre      = mtx[rsi - 1][0];
        Matrix nxt      = mtx[rsi][0];
        nxt.context     = pre.context;
        nxt.context.slide(rs[rsi], '-', kmer_size);
        size_t left_con = base2int(rs[rsi]) * 5 + base2int('-');
        nxt.score       = pre.score + st[ti(pre.context.con, left_con)];
        nxt.tb_dir      = TB_UP;
    }

    for(size_t qsi = 1; qsi < ql; qsi++){
        Matrix pre      = mtx[0][qsi - 1];
        Matrix nxt      = mtx[0][qsi];
        nxt.context     = pre.context;
        nxt.context.slide('-', qs[qsi - 1], kmer_size);
        size_t upper_con= base2int('-') * 5 + base2int(qs[qsi]);
        nxt.score       = pre.score + st[ti(pre.context.con, upper_con)];
        nxt.tb_dir      = TB_LEFT;
    }

    for(size_t rmi = 1; rmi <= rl; rmi++){
		for(size_t qmi = 1; qmi <= ql; qmi++){
            //fprintf(stderr, "proccecing (rmi, qmi) = (%d, %d)\n", rmi, qmi);
            //showMatrix(mtx);
            const size_t rsi = rmi - 1; // reference index on reference_sequence, i.e.   rs[rp]
            const size_t qsi = qmi - 1; // query index on query_sequence,         i.e.   qs[qp]
            Matrix  upper_left_matrix = mtx[rmi - 1][qmi - 1];
            Matrix        left_matrix = mtx[rmi - 1][qmi    ];
            Matrix       upper_matrix = mtx[rmi    ][qmi - 1];
            size_t         match_case = base2int(rs[rsi]) * 5 + base2int(qs[qsi]);//
            size_t          left_case = base2int(rs[rsi]) * 5 + base2int(    '-');//next pair's index
            size_t         upper_case = base2int(    '-') * 5 + base2int(qs[qsi]);//
            const int32_t     m_score = upper_left_matrix.score + st[ti(upper_left_matrix.context.con, match_case)];
            const int32_t  left_score =       left_matrix.score + st[ti(      left_matrix.context.con,  left_case)];
            const int32_t upper_score =      upper_matrix.score + st[ti(     upper_matrix.context.con, upper_case)];

            TraceBackDirection current_best_traceback_direction = TB_START;
            TraceBackDirection& cbt = current_best_traceback_direction;
            //fprintf(stderr, "%d, %d, %d\n", m_score, left_score, upper_score);
			int32_t score_max = m_score;
            cbt = TB_MATCH;

            if(score_max < left_score){
				score_max = left_score;
                cbt = TB_LEFT;
            }
            if(score_max < upper_score){
				score_max = upper_score;
                cbt = TB_UP;
            }
            //fprintf(stderr, "%d\n", score_max);
            Matrix& current_cell = mtx[rmi][qmi];
            Matrix& cc = current_cell;
            cc.score = score_max;
            cc.tb_dir = cbt;
            if(cbt == TB_MATCH){
                upper_left_matrix.context.slide(rs[rsi], qs[qsi], kmer_size);
                cc.context = upper_left_matrix.context;
            } else if(cbt == TB_LEFT){
                left_matrix.context.slide(rs[rsi],     '-', kmer_size);
                cc.context = left_matrix.context;
            } else if(cbt == TB_UP){
                upper_matrix.context.slide(    '-', qs[qsi], kmer_size);
                cc.context = upper_matrix.context;
            }else if(cbt == TB_START){
                /* new_reference_subseq contains previous k-mer string of reference sequence.
                 * new_query_subseq     contains previous k-mer string of query     sequence.
                 *
                 *                      rmi
                 *     . . . A   C   T   G   A
                 *     .                 .
                 *     .
                 *     C                 .
                 *
                 *     C                 .
                 *
                 * qmi T  .  .  .  .  .  0 <----result of recursive equivalation.
                 *
                 *     A
                 *
                 *
                 *
                 * In this case, 3-mer of reference and query sequence is "CTG" and "CCT", respectively.
                 *
                 */

                Context dummy_context;//when smith-waterman algorithm end with low score, use this context. need to be filled with something
                dummy_context.con   = all_gap_index;
                for(size_t i = 0; i < kmer_size; i++){
                    int dif = kmer_size - i;
                    char ref_buf, qry_buf;
                    if(rmi - dif >= 0 && rmi - dif < rl){
                        ref_buf = rs[rmi - dif];
                    }else{
                        ref_buf = '-';
                    }
                    if(qmi - dif >= 0 && qmi - dif < ql){
                        qry_buf = qs[qmi - dif];
                    }else{
                        qry_buf = '-';
                    }
                    dummy_context.slide(ref_buf, qry_buf, kmer_size);
                }
                cc.context = dummy_context;
            } else {
                //fprintf(stderr, "current_best_traceback_direction = %d\n", cbt);
                ASSERT(0, "You should not come here!");
            }
            //fprintf(stderr, "(rmi, qmi) = (%d, %d), (cc.context.ref_seq, cc.context.query_seq) = (%s, %s)\n", rmi, qmi, cc.context.ref_seq, cc.context.query_seq);
		}
	}

    //showMatrix(mtx);
    return mtx[rl][ql];
}

void load_from_score_table(
        int32_t **loaded_score_table,
        const char *score_table_file_name,
        size_t *kk
) {
    FILE *fle = fopen(score_table_file_name, "rb");
    if(fle == nullptr) {
        fprintf(stderr, "ERROR: cannot open the score table file '%s'\n", score_table_file_name);
        exit(2);
    }
    //// FILE FORMAT:
    ////  offset  0: 'BINFREQT' (8 bytes)
    ////  offset  8: k-mer size (8 bytes)
    ////  offset 16: sizeof(Frequency) (8 bytes)
    ////  offset 24: table_size (8 bytes)
    ////  offset 32: table (8 x table_size bytes)
    static const char MAGIC_STRING[] = "BINFREQT";
    char buffer[256];
    fread(buffer, sizeof(MAGIC_STRING) - 1, 1, fle);
    if(strncmp(buffer, MAGIC_STRING, sizeof(MAGIC_STRING) - 1)) {
        fprintf(stderr, "ERROR: wrong file. could not find MAGIC_HEADER\n");
        exit(2);
    }
    if(sizeof(size_t) != 8) {
        fprintf(stderr, "ERROR: sizeof(size_t) must be 8\n");
        exit(2);
    }
    fread(kk, sizeof(*kk), 1, fle);
    if(*kk < 1 || 6 < *kk) {
        fprintf(stderr, "ERROR: kk is out of range (kk = %zu)\n", *kk);
        exit(2);
    }
    size_t size_of_frequency = 0;
    fread(&size_of_frequency, sizeof(size_of_frequency), 1, fle);
    if(size_of_frequency == 8) {
        fprintf(stderr, "Frequency must be 8-byte integer\n");
        exit(2);
    }
    size_t var_table_size;
    fread(&var_table_size, sizeof(var_table_size), 1, fle);
    if(var_table_size != ipow(5, 2 * (*kk))) {
        fprintf(stderr, "ERROR: table_size is not valid. (actual = %zu, expected = %d)\n", var_table_size, ipow(5, 2 * (*kk)));
        exit(2);
    }

    *loaded_score_table = (int32_t*)malloc(size_of_frequency * var_table_size);
    if(*loaded_score_table == nullptr) {
        fprintf(stderr, "ERROR: out of memory\n");
        exit(2);
    }
    fread(*loaded_score_table, size_of_frequency * var_table_size, 1, fle);
    fclose(fle);
}

#ifdef MAIN
int main(int argc, char *argv[]){
	if(argc < 4){
		printf("few args\n");
		return 1;
	}
    //fprintf(stderr, "begin to set values\n");
    const char* const reference_sequence = argv[1];
    const char* const query_sequence = argv[2];
    const int k = std::stoi(argv[3]);
    const char* const score_table_file_name = argv[4];
    int32_t *loaded_score_table = nullptr;

    if(score_table_file_name != nullptr) {
        size_t kk;
        load_from_score_table(&loaded_score_table, score_table_file_name, &kk);
        if(kk != k) {
            fprintf(stderr, "ERROR: the score table has K = %zu, but you specified K = %d\n", kk, k);
            exit(2);
        }
    }

    if(loaded_score_table != nullptr) {
        Matrix result = needlman_wunsch(reference_sequence, query_sequence, k, loaded_score_table);
        printf("%d\n", result.score);
    } else {
        switch(k) {
            case 0:
                {
                    Matrix result = needlman_wunsch(reference_sequence, query_sequence, 1, score_1mer);
                    printf("%d\n", result.score);
                }
                break;
            case 1:
                {
                    Matrix result = needlman_wunsch(reference_sequence, query_sequence, 1, score_1mer_mod_3digit);
                    printf("%d\n", result.score);
                }
                break;
            case 2:
                {
                    Matrix result = needlman_wunsch(reference_sequence, query_sequence, 2, score_2mer_mod_3digit);
                    printf("%d\n", result.score);
                }
                break;
            case 3:
                {
                    Matrix result = needlman_wunsch(reference_sequence, query_sequence, 3, score_3mer_mod_3digit);
                    printf("%d\n", result.score);
                }
                break;
            case 4:
                {
                    Matrix result = needlman_wunsch(reference_sequence, query_sequence, 4, score_4mer_mod_3digit);
                    printf("%d\n", result.score);
                }
                break;
            default:
                fprintf(stderr, "Should never reach here!\n");
                exit(2);
                break;
        }
    }
	return 0;
}
#endif


#ifdef TEST
void test(
    const char    *reference_sequence,
    const char    *query_sequence,
    const size_t  kmer_size,
    const int32_t *score_table,
    const int32_t best_score,
    const int64_t best_score_reference_sequence_index,//1 origin
    const int64_t best_score_query_sequence_index// 1 origin
    ){
    const char    *refseq = reference_sequence;
    const char    *qryseq = query_sequence;
    const int32_t *st     = score_table;
    const int32_t bs      = best_score;
    const int32_t bsr     = best_score_reference_sequence_index;
    const int32_t bsq     = best_score_query_sequence_index;
    fprintf(stderr, "ref   : %s\nquery : %s\n", refseq, qryseq);
    Matrix test_result = smith_waterman(refseq, qryseq, kmer_size, st);
    Matrix tr = test_result;
    ASSERT( bs  == tr.score, "score mismatch");
    ASSERT( bsr == tr.ipos , "reference position mitchmatch");
    ASSERT( bsq == tr.jpos,  "query position mismatch");
    fprintf(stderr, "finish checking kmer_size : %zu, \"%s\", \"%s\", %d, %d, %d \n\n", kmer_size, refseq, qryseq, bs, bsr, bsq);
}
int main(void){
    fprintf(stderr, "=== BEGIN TEST ===\n");
    //       ref, query, kmer_size, score_table, best score, ref_pos, query_pos

    test("A", "A",     1, score_1mer, 1, 1, 1);
    test("AA", "AA",   1, score_1mer, 2, 2, 2);
    test("AC", "AAC",  1, score_1mer, 2, 2, 3);
    test("TT", "GGGG", 1, score_1mer, 0, 0, 0);
    test("GGGAAAGGG", "TTTTTTAAATTTTTT", 1, score_1mer, 3, 6, 9);

    test("A", "A",     1, score_1mer_mod, 3, 1, 1);
    test("AA", "AA",   1, score_1mer_mod, 6, 2, 2);
    test("ACG", "ACGT",  1, score_1mer_mod, 9, 3, 3);
    test("AC", "AAC",  1, score_1mer_mod, 6, 2, 3);
    test("TT", "GGGG", 1, score_1mer_mod, 2, 2, 4);

    test("A", "A",     2, score_2mer, 1, 1, 1);
    test("AA", "AA",   2, score_2mer, 2, 2, 2);
    test("TC", "TTC",  2, score_2mer, 1, 2, 3);
    test("TT", "GGGG", 2, score_2mer, 0, 0, 0);
    test("GGGAAAGGG", "TTTTTTAAATTTTTT", 2, score_2mer, 2, 6, 9);
/*
    test("A", "A",     3, score_3mer, 1, 1, 1);
    test("AA", "AA",   3, score_3mer, 2, 2, 2);
    test("TC", "TTC",  3, score_3mer, 1, 2, 3);
    test("TT", "GGGG", 3, score_3mer, 0, 0, 0);
    test("GGGAAAGGG", "TTTTTTAAATTTTTT", 3, score_3mer, 2, 6, 9);
*/

    test("A", "A",     2, score_2mer_mod, 0, 0, 0);
    test("AA", "AA",   2, score_2mer_mod, 2, 2, 2);
    test("TC", "TTC",  2, score_2mer_mod, 2, 2, 3);
    test("TT", "GGGG", 2, score_2mer_mod, 0, 0, 0);
    test("ACG", "ACGT", 2, score_2mer_mod, 4, 3, 4);

    fprintf(stderr, "===  END TEST  ===\n");
    return 0;
}

#endif
