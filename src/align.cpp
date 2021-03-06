
// vim:ff=unix ft=cpp ts=4 sw=4 sts=4 si et :
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <string>
#include <vector>
#include "score_1mer.h"
#include "score_2mer.h"
#include "score_3mer.h"
#include "score_1mer_mod.h"
#include "score_2mer_mod.h"
#include "score_2mer_mod_3digit.h"

#define _pos(i, j, ilen)	( (j) * ((ilen) + 1) + (i) )
#define MAX2(p, q)			( ((p) < (q)) ? (q) : (p) )
#define MAX3(p, q, r)		( MAX2(p, MAX2(q, r)) )
#define MAX4(p, q, r, s)	( MAX2(MAX2(p, q), MAX2(r, s)) )

#define ASSERT(cond, message) do { if(!(cond)){ fprintf(stderr, "ASSERTION ERROR: %s\nCONDITION: %s\n", message, #cond); exit(2); } } while(0)


using namespace std;

uint32_t shift_digit(uint32_t pre, uint32_t add, size_t kmer_size){
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
	size_t ref_seq;  ///< Reference context sequence (integer of 5-ord numbers)
    size_t query_seq;///< Query context sequence  (integer of 5-ord numbers)
    Context() {
        ref_seq = query_seq = 0;
    }

    /**
     * Slide context
     *
     * for example:
     *   when context = ref: ACGT,
     *                  qry: CGGC,
     *   calling slide('C', 'A') will lead to the following result:
     *        context = ref: CGTC,
     *                  qry: GGCA
     */
    void slide(char ref_char, char query_char, int kmer_size){
        ref_seq   = shift_digit(ref_seq,   base2int(ref_char), kmer_size);
        query_seq = shift_digit(query_seq, base2int(query_char), kmer_size);
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
	TraceBackDirection   tb_dir;   ///< Traceback direction.(32bit)
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
	int sum = 0;
	size_t len = strlen(str);
	while(str[i] != '\0'){
		sum += base2int(str[i]) * (int)pow(5, len - i - 1);
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

void showMatrix(const vector<vector<Matrix> >& mtx)
{
    fprintf(stderr, "======\n");
    fprintf(stderr, "score\n");
    for(size_t y = 0; y < mtx.size(); y++) {
        const vector<Matrix>& line = mtx[y];
        for(size_t x = 0; x < line.size(); x++) {
            fprintf(stderr, "%3d ", line[x].score);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "context.ref_seq\n");
    for(size_t y = 0; y < mtx.size(); y++) {
        const vector<Matrix>& line = mtx[y];
        for(size_t x = 0; x < line.size(); x++) {
            fprintf(stderr, "%3zu ", line[x].context.ref_seq);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "context.query_seq\n");
    for(size_t y = 0; y < mtx.size(); y++) {
        const vector<Matrix>& line = mtx[y];
        for(size_t x = 0; x < line.size(); x++) {
            fprintf(stderr, "%3zu ", line[x].context.query_seq);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "tb_dir\n");
    for(size_t y = 0; y < mtx.size(); y++) {
        const vector<Matrix>& line = mtx[y];
        for(size_t x = 0; x < line.size(); x++) {
            fprintf(stderr, "%3d ", line[x].tb_dir);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "======\n");
}


Matrix smith_waterman(
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
    int all_gap_index;
    {
        char gaps[kmer_size + 1];
        for(size_t i = 0; i < kmer_size; i++){
            gaps[i] = '-';
        }
        gaps[kmer_size] = '\0';
        all_gap_index = str2idx(gaps);
    }

    vector<vector<Matrix> > mtx(rl + 1, vector<Matrix>(ql + 1));
    // initialize
    mtx[0][0].context.ref_seq   = all_gap_index;
    mtx[0][0].context.query_seq = all_gap_index;

    for(size_t rmi = 0; rmi <= rl; rmi++){
        mtx[rmi][0].tb_dir = TB_LEFT;
	}
    for(size_t rsi = 0; rsi < rl; rsi++){
        mtx[rsi + 1][0].context.slide(rs[rsi], '-', kmer_size);
    }

	for(size_t qmi = 0; qmi <= ql; qmi++){
        mtx[0][qmi].tb_dir = TB_UP;
	}
    for(size_t qsi = 0; qsi < ql; qsi++){
        mtx[0][qsi + 1].context.slide('-', qs[qsi], kmer_size);
    }

    mtx[0][0].tb_dir = TB_START;
    //fprintf(stderr, "finish to initialize\n");
    int32_t current_best_score = 0;
    int32_t current_best_rmi   = 0;
    int32_t current_best_qmi   = 0;

    const int32_t *st = kmer_score_table;

    /*
     * mtx_col_len : size of score table
     * when k = 1, mtx_col_len should be 5
     * when k = 2, mtx_col_len should be 25
     * when k = 3, mtx_col_len should be 125
     */
    int mtx_col_len = 1;
    for(size_t i = 0; i < kmer_size; i++){
        mtx_col_len *= 5;
    }
    #define table_index(rmi, qmi) ((qmi) * mtx_col_len + (rmi))
    #define tb(i, j) table_index(i, j)

    for(size_t rmi = 1; rmi <= rl; rmi++){
		for(size_t qmi = 1; qmi <= ql; qmi++){
            //fprintf(stderr, "proccecing (rmi, qmi) = (%d, %d)\n", rmi, qmi);
            //showMatrix(mtx);
            const size_t rsi = rmi - 1; // reference index on reference_sequence, i.e.   rs[rp]
            const size_t qsi = qmi - 1; // query index on query_sequence,         i.e.   qs[qp]
            Matrix upper_left_matrix = mtx[rmi - 1][qmi - 1]; upper_left_matrix.context.slide(rs[rsi], qs[qsi], kmer_size);
            Matrix left_matrix       = mtx[rmi - 1][qmi    ]; left_matrix.context.slide(      rs[rsi],     '-', kmer_size);
            Matrix upper_matrix      = mtx[rmi    ][qmi - 1]; upper_matrix.context.slide(         '-', qs[qsi], kmer_size);
            //fprintf(stderr, "(%d, %d)finish sliding context\n", rmi, qmi);

            const int32_t     m_score = upper_left_matrix.score + st[tb(upper_left_matrix.context.ref_seq, upper_left_matrix.context.query_seq)];
            const int32_t  left_score =       left_matrix.score + st[tb(      left_matrix.context.ref_seq,       left_matrix.context.query_seq)];
            const int32_t upper_score =      upper_matrix.score + st[tb(     upper_matrix.context.ref_seq,      upper_matrix.context.query_seq)];
            //fprintf(stderr, "(%d, %d)finish calculating scores\n", rmi, qmi);
            TraceBackDirection current_best_traceback_direction = TB_START;
            TraceBackDirection& cbt = current_best_traceback_direction;
			int32_t score_max = 0;
            //fprintf(stderr, "(m, left, up) = (%d, %d, %d)\n", m_score, left_score, upper_score);
			if(score_max < m_score){
				score_max = m_score;
                cbt = TB_MATCH;
			}
            if(score_max < left_score){
				score_max = left_score;
                cbt = TB_LEFT;
            }
            if(score_max < upper_score){
				score_max = upper_score;
                cbt = TB_UP;
            }
            if(score_max == 0){
                cbt = TB_START;
            }

            Matrix& current_cell = mtx[rmi][qmi];
            Matrix& cc = current_cell;
            cc.score = score_max;
            cc.tb_dir = cbt;
            if(cbt == TB_MATCH){
                cc.context = upper_left_matrix.context;
            } else if(cbt == TB_LEFT){
                cc.context = left_matrix.context;
            } else if(cbt == TB_UP){
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
                dummy_context.ref_seq   = all_gap_index;
                dummy_context.query_seq = all_gap_index;
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
			if(score_max >= current_best_score){
                current_best_score = score_max;
                current_best_rmi = rmi;
                current_best_qmi = qmi;
			}
		}
	}
    const int32_t best_score = current_best_score;
    const int32_t best_rmi   = current_best_rmi;
    const int32_t best_qmi   = current_best_qmi;

    ASSERT(best_score == mtx[best_rmi][best_qmi].score, "Logic error");
    //showMatrix(mtx);
    Matrix smith_waterman_result;
    Matrix sw = smith_waterman_result;
    if(best_score > 0){
        sw.score = best_score;
        sw.ipos  = best_rmi;
        sw.jpos  = best_qmi;
    }else{
        sw.score = 0;
        sw.ipos  = 0;
        sw.jpos  = 0;
    }
    return sw;
}

#ifdef MAIN
int main(int argc, char *argv[]){
	if(argc < 3){
		printf("few args\n");
		return 1;
	}
    //fprintf(stderr, "begin to set values\n");
    const char* const reference_sequence = argv[1];
    const char* const query_sequence = argv[2];
    Matrix result = smith_waterman(reference_sequence, query_sequence, 1, score_1mer_mod);
    printf("%d\n", result.score);

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
