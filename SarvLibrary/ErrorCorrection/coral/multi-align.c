/**
 *  Coral: short reads error correction with multiple alignments
 *  Copyright (C) 2011 Leena Salmela <leena.salmela@cs.helsinki.fi>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "define.h"
#include "multi-align.h"
#include "reverse.h"

int quick = 0;
int banded = 0;
int full = 0;

/* Compare two read_pos structures according to read ids */
int compare_read(const void *s1, const void *s2) {
    read_pos *r1, *r2;

    r1 = (read_pos *)s1;
    r2 = (read_pos *)s2;

    if (r1->read - r2->read == 0) {
        return r1->occ - r2->occ;
    } else {
        return r1->read - r2->read;
    }
}

/* Compare two read_pos structures according to positions */
int compare_pos(const void *s1, const void *s2) {
    read_pos *r1, *r2;

    r1 = (read_pos *)s1;
    r2 = (read_pos *)s2;

    return r1->pos - r2->pos;
}

/* Compare two read_pos structures according to counts */
int compare_count(const void *s1, const void *s2) {
    read_pos *r1, *r2;

    r1 = (read_pos *)s1;
    r2 = (read_pos *)s2;

    return r2->count - r1->count;
}

/* Update the consensus sequence of an alignment starting at position
   start and ending at position end */
void update_consensus(align *alignment, int start, int end) {
    int i;

    if (end >= MAX_CONTIG_LEN) {
        printf("multi-align: max contig len too short.\n");
        exit(1);
    }

#ifdef DEBUG_ALIGN
    printf("Updating consensus: %d -> %d\n", start, end);
#endif

    for(i = start; i < end; i++) {
        if ((alignment->contigA[i] == 0 && alignment->contigC[i] == 0 && 
             alignment->contigG[i] == 0 && alignment->contigT[i] == 0) ||
            (alignment->contigN[i] > alignment->contigA[i] && 
             alignment->contigN[i] > alignment->contigC[i] &&
             alignment->contigN[i] > alignment->contigG[i] &&
             alignment->contigN[i] > alignment->contigT[i])) {
            alignment->consensus[i] = '-';
        } else if (alignment->contigA[i] > alignment->contigC[i] &&
                   alignment->contigA[i] > alignment->contigG[i] &&
                   alignment->contigA[i] > alignment->contigT[i]) {
            alignment->consensus[i] = 'A';
        } else if (alignment->contigC[i] > alignment->contigG[i] &&
                   alignment->contigC[i] > alignment->contigT[i]) {
            alignment->consensus[i] = 'C';
        } else if (alignment->contigG[i] > alignment->contigT[i]) {
            alignment->consensus[i] = 'G';
        } else {
            alignment->consensus[i] = 'T';
        }
    }
}

/**
 * Align one more read against the alignment.
 * - alignment: The alignment computed so far
 * - reads: The set of reads
 * - pos: The set of reads for which the alignment is computed
 * - end: Length of alignment
 * - p: The id to pos of the new read to align (the reads with smaller ids have already been aligned)
 * Returns the length of the alignment after the new read has been added.
 */
int align_read(align *alignment, char **reads, char **quals, read_pos *pos, 
               int end,  int p, double max_error_rate,
               int match_reward, int mm_penalty, int gap_penalty) {
    int i,j;
    int m;

    char *read = reads[pos[p].read];
    char *qual = (quals == NULL) ? NULL: quals[pos[p].read];

    int len;

    char r[MAX_READ_LENGTH];
    char rq[MAX_READ_LENGTH];
    char er[MAX_READ_LENGTH];
    char erq[MAX_READ_LENGTH];

    int k;
    int kk;
    int m2;

    int edits;
    int old_end;

#pragma omp critical
    {
    m = strlen(read);
    if (pos[p].ori == 'C') {
        // Form the reverse complement before aligning
        for(i = 0; i < (signed int)m; i++) {
            switch(read[i]) {
            case 'a':
            case'A':
                r[m-i-1] = 'T';
                break;
            case 'c':
            case'C':
                r[m-i-1] = 'G';
                break;
            case 'g':
            case'G':
                r[m-i-1] = 'C';
                break;
            case 't':
            case'T':
                r[m-i-1] = 'A';
                break;
            case 'N':
            case 'n':
                r[m-i-1] = 'N';
                break;

            }
            if (quals != NULL)
                rq[m-i-1] = qual[i];
        }
        r[m] = '\0';
        read = r;
        if (quals != NULL) {
            rq[m] = '\0';
            qual = rq;
        }
    } else {
        // Make a local copy
        for(i = 0; i < (signed int)m; i++) {
            r[i] = read[i];
            if (quals != NULL)
                rq[i] = qual[i];
        }
        r[m] = '\0';
        read = r;
        if (quals != NULL) {
            rq[m] = '\0';
            qual = rq;
        }
    }
    }

    if (end <= 0) {
        // length of consensus is 0
        // update base counts
        for(i = 0; i < (signed int)m; i++) {
            alignment->contigA[i] = 0;
            alignment->contigC[i] = 0;
            alignment->contigG[i] = 0;
            alignment->contigT[i] = 0;
            alignment->contigN[i] = 0;
            switch(read[i]) {
            case 'a':
            case 'A':
                alignment->contigA[i]++;
                break;
            case 'c':
            case 'C':
                alignment->contigC[i]++;
                break;
            case 'g':
            case 'G':
                alignment->contigG[i]++;
                break;
            case 't':
            case 'T':
                alignment->contigT[i]++;
                break;
            case '-':
            case'N':
            case 'n':
                alignment->contigN[i]++;
                break;
            }
        }
        alignment->offset = pos[p].pos;
        pos[p].pos = 0;
        strcpy(pos[p].edited_read, reads[pos[p].read]);
        if (quals != NULL) {
            memcpy(pos[p].edited_qual, quals[pos[p].read], m);
            pos[p].edited_qual[m] = '\0';
        }
        pos[p].aligned = 1;
        return m;
    }

    update_consensus(alignment, 0, end);
    alignment->consensus[end] = '\0';

    pos[p].pos = pos[p].pos - alignment->offset;


#ifdef DEBUG_ALIGN
    printf("Aligning read %s (%d, %c, %d, %d)\n", read, pos[p].read, pos[p].ori, pos[p].pos, pos[p].count);
    printf("Against conse %s\n", alignment->consensus);
#endif

    if (gap_penalty > 100*mm_penalty) {
        /* There won't be gaps... */

        for(i = 0, j = pos[p].pos, edits = 0;
            i < m && j < end; i++, j++) {
            if (j >= 0) {
                if (alignment->consensus[j] != toupper(read[i]))
                    edits++;
            }
        }

        if (edits > m*max_error_rate) {
#ifdef DEBUG_ALIGN
            printf("Too many errors in read %d %c\n", pos[p].read, pos[p].ori);
#endif
            /* alignment->ok = 0; */
            return end;
        }

    } else {

        /* Try a quick gapless alignment positioned according to the k-mer*/

        if (pos[p].pos > 0) {
            for(k = 0, j = 0; k < pos[p].pos; j++) {
                if (alignment->consensus[j] != '-')
                    k++;
            }
	    while(alignment->consensus[j] == '-')j++;

            pos[p].pos = j;
        }

        for(i = 0, j = pos[p].pos, edits = 0;
            i < m && j < end && edits < 1; i++, j++) {
            if (j >= 0) {
                if (alignment->consensus[j] == '-') {
                    i--;
                } else {
                    if (alignment->consensus[j] != toupper(read[i]))
                        edits++;
                }
            }
        }
    }

    if (edits < 1 || gap_penalty > 100*mm_penalty) {
        quick++;
        /* Quick alignment found! */
#ifdef DEBUG_ALIGN
        printf("Quick alignment found for read %d (%c): %d\n", pos[p].read, 
               pos[p].ori, pos[p].pos);
#endif
        if (pos[p].pos < 0) {
            /* Insert columns to the beginning of the alignment */
            for (k = end; k >= 0; k--) {
                alignment->contigA[k-pos[p].pos] = alignment->contigA[k];
                alignment->contigC[k-pos[p].pos] = alignment->contigC[k];
                alignment->contigG[k-pos[p].pos] = alignment->contigG[k];
                alignment->contigT[k-pos[p].pos] = alignment->contigT[k];
                alignment->contigN[k-pos[p].pos] = alignment->contigN[k];
                alignment->consensus[k-pos[p].pos] = alignment->consensus[k];
            }
            for(k = 0; k < -pos[p].pos; k++) {
                alignment->contigA[k] = 0;
                alignment->contigC[k] = 0;
                alignment->contigG[k] = 0;
                alignment->contigT[k] = 0;
                alignment->contigN[k] = 0;
                alignment->consensus[k] = 'A';
            }
            end += -pos[p].pos;
            alignment->offset += pos[p].pos;
            for(k = 0; k < p; k++) {
                if (pos[k].aligned)
                    pos[k].pos -= pos[p].pos;
            }
            pos[p].pos -= pos[p].pos;
        }

        for(i = 0, j = pos[p].pos, k=0; i < m && j < end; i++, j++, k++) {
            if (alignment->consensus[j] == '-' && 
                gap_penalty <= 100*mm_penalty) {
                er[k] = '-';
                if (quals != NULL) {
                    if (i > 0 ) {
                        erq[k] = (qual[i-1] + qual[i])/2;
                    } else {
                        erq[k] = qual[i];
                    }
                }
                i--;
                alignment->contigN[j]++;
            } else {
                er[k] = read[i];
                if (quals != NULL)
                    erq[k] = qual[i];
                switch(read[i]) {
                case 'A':
                case 'a':
                    alignment->contigA[j]++;
                break;
                case 'C':
                case 'c':
                    alignment->contigC[j]++;
                break;
                case 'G':
                case 'g':
                    alignment->contigG[j]++;
                break;
                case 'T':
                case 't':
                    alignment->contigT[j]++;
                break;
                case 'N':
                case 'n':
                case '-':
                    alignment->contigN[j]++;
                break;
                }
            }
        }
        for(; i < m; i++, j++, k++) {
            end++;
            alignment->contigA[j] = 0;
            alignment->contigC[j] = 0;
            alignment->contigG[j] = 0;
            alignment->contigT[j] = 0;
            alignment->contigN[j] = 0;
            er[k] = read[i];
            if (quals != NULL)
                erq[k] = qual[i];
            switch(read[i]) {
            case 'A':
            case 'a':
                alignment->contigA[j]++;
            break;
            case 'C':
            case 'c':
                alignment->contigC[j]++;
            break;
            case 'G':
            case 'g':
                alignment->contigG[j]++;
            break;
            case 'T':
            case 't':
                alignment->contigT[j]++;
            break;
            case 'N':
            case 'n':
            case '-':
                alignment->contigN[j]++;
            break;
            }
        }
        er[k] = '\0';

#ifdef DEBUG_ALIGN
        printf("ER: %s\n", er);
#endif

        if (quals != NULL)
            erq[k] = '\0';
        if (pos[p].ori == 'U') {
            for(i = 0; i < k; i++){
                pos[p].edited_read[i] = er[i];
                if (quals != NULL)
                    pos[p].edited_qual[i] = erq[i];
            }
            pos[p].edited_read[i] = '\0';
            if (quals != NULL)
                pos[p].edited_qual[i] = '\0';
        } else {
            for(i = 0; i < k; i++) {
                if (quals != NULL)
                    pos[p].edited_qual[k-i-1] = erq[i];
                switch(er[i]) {
                case 'A':
                case 'a':
                    pos[p].edited_read[k-i-1] = 'T';
                break;
                case 'C':
                case 'c':
                    pos[p].edited_read[k-i-1] = 'G';
                break;
                case 'G':
                case 'g':
                    pos[p].edited_read[k-i-1] = 'C';
                break;
                case 'T':
                case 't':
                    pos[p].edited_read[k-i-1] = 'A';
                break;
                default:
                    pos[p].edited_read[k-i-1] = read[i];
                }
            }
            pos[p].edited_read[i] = '\0';
            if (quals != NULL)
                pos[p].edited_qual[i] = '\0';
        }

        pos[p].aligned = 1;
        return end;
    }

    if (pos[p].occ == 1) {
        banded++;
        
        /* The k-mer occurs only once -> use banded alignment */
        int s, e;

        s = pos[p].pos - (max_error_rate * m + 1);
        e = pos[p].pos + (max_error_rate * m + 1);

#ifdef DEBUG_ALIGN
        printf("Band: %d -> %d\n", s, e);
#endif

        for(i = s; i <= e; i++) {
            if (i >= 0) {
                if (alignment->consensus[i] == '-') {
                    e++;
                }
                if (e >= end) {
                    break;
                }
            }
        }

        for(j = 0; j < m; j++) {

            if (s < end) {
                int gapi;
                int gapj;

                if (s >= 0) {
                    gapi = (s == end-1)? 0 : gap_penalty;
                    if (alignment->consensus[s] == toupper(read[j])) {
                        if (alignment->dp[s][j] + match_reward >= alignment->dp[s+1][j] - gapi) {
                            alignment->dp[s+1][j+1] = alignment->dp[s][j] + match_reward;
                            alignment->dp_trace[s+1][j+1] = 'M';
                        } else {
                            alignment->dp[s+1][j+1] = alignment->dp[s+1][j] - gapi;
                            alignment->dp_trace[s+1][j+1] = 'I';
                        }
                    } else {
                        if (alignment->dp[s][j] - mm_penalty >= alignment->dp[s+1][j] - gapi) {
                            alignment->dp[s+1][j+1] = alignment->dp[s][j] - mm_penalty;
                            alignment->dp_trace[s+1][j+1] = 'D';
                        } else {
                            alignment->dp[s+1][j+1] = alignment->dp[s+1][j] - gapi;
                            alignment->dp_trace[s+1][j+1] = 'I';
                        }
                    }
                }

                for(i = (s < 0) ? 0 : s+1; i < ((e > end) ? end : e); i++) {
                    gapi = (i == end-1)? 0 : gap_penalty;
                    gapj = (j == m-1 || alignment->consensus[i] == '-')? 0 : gap_penalty;
                    if (alignment->consensus[i] == toupper(read[j])) {
                        if (alignment->dp[i][j] + match_reward >= alignment->dp[i][j+1] - gapj  &&
                            alignment->dp[i][j] + match_reward >= alignment->dp[i+1][j] - gapi)  {
                            alignment->dp[i+1][j+1] = alignment->dp[i][j] + match_reward;
                            alignment->dp_trace[i+1][j+1] = 'M';
                        } else if (alignment->dp[i][j+1] - gapj >= alignment->dp[i+1][j] - gapi) {
                            alignment->dp[i+1][j+1] = alignment->dp[i][j+1] - gapj;
                            alignment->dp_trace[i+1][j+1] = 'J';
                        } else {
                            alignment->dp[i+1][j+1] = alignment->dp[i+1][j] - gapi;
                            alignment->dp_trace[i+1][j+1] = 'I';
                        }
                    } else {
                        if (alignment->dp[i][j] - mm_penalty >= alignment->dp[i][j+1] - gapj  &&
                            alignment->dp[i][j] - mm_penalty >= alignment->dp[i+1][j] - gapi)  {
                            alignment->dp[i+1][j+1] = alignment->dp[i][j] - mm_penalty;
                            alignment->dp_trace[i+1][j+1] = 'D';
                        } else if (alignment->dp[i][j+1] - gapj >= alignment->dp[i+1][j] - gapi) {
                            alignment->dp[i+1][j+1] = alignment->dp[i][j+1] - gapj;
                            alignment->dp_trace[i+1][j+1] = 'J';
                        } else {
                            alignment->dp[i+1][j+1] = alignment->dp[i+1][j] - gapi;
                            alignment->dp_trace[i+1][j+1] = 'I';
                        }
                    }
                }

                if (e < end && e >= 0) {
                    gapj = (j == m-1 || alignment->consensus[e] == '-')? 0 : gap_penalty;
                    if (alignment->consensus[e] == toupper(read[j])) {
                        if (alignment->dp[e][j] + match_reward >= alignment->dp[e][j+1] - gapj) {
                            alignment->dp[e+1][j+1] = alignment->dp[e][j] + match_reward;
                            alignment->dp_trace[e+1][j+1] = 'M';
                        } else {
                            alignment->dp[e+1][j+1] = alignment->dp[e][j+1] - gapj;
                            alignment->dp_trace[e+1][j+1] = 'J';
                        }
                    } else {
                        if (alignment->dp[e][j] - mm_penalty >= alignment->dp[e][j+1] - gapj) {
                            alignment->dp[e+1][j+1] = alignment->dp[e][j] - mm_penalty;
                            alignment->dp_trace[e+1][j+1] = 'D';
                        } else {
                            alignment->dp[e+1][j+1] = alignment->dp[e][j+1] - gapj;
                            alignment->dp_trace[e+1][j+1] = 'J';
                        }
                    }
                }
            } else {
                alignment->dp[end][j+1] = alignment->dp[end][j];
                alignment->dp_trace[end][j+1] = 'I';
            }

            s++;
            e++;
            while(s >= 0 && s < end-1 && alignment->consensus[s] == '-') {
                s++;
		/* alignment->dp[s][j+1] = alignment->dp[s-1][j+1]; */
		/* alignment->dp_trace[s][j+1] = 'J'; */
            }
            while(e >= 0 && e < end-1 && alignment->consensus[e] == '-') {
                e++;
		alignment->dp[e][j+1] = alignment->dp[e-1][j+1];
		alignment->dp_trace[e][j+1] = 'J';
            }
        }

        for(i = e-1; i < end; i++) {
            alignment->dp[i+1][m] = alignment->dp[i][m];
            alignment->dp_trace[i+1][m] = 'J';
        }

    } else {
        full++;

        /* The k-mer occurs several times -> compute full alignment */

        /* Align the read against the consensus */
        /*
        for(i = 0; i <= end; i++) {
            alignment->dp[i][0] = 0;
            alignment->dp_trace[i][0] = 'J';
        }
        for(j = 0; j <= m; j++) {
            alignment->dp[0][j] = 0;
            alignment->dp_trace[0][j] = 'I';
        }
        */

        for(i = 0; i < end; i++) {
            int gapi = (i == end-1)? 0 : gap_penalty;
            for(j = 0; j < m; j++) {
                int gapj = (j == m-1 || alignment->consensus[i] == '-')? 0 : gap_penalty;
                if (alignment->consensus[i] == toupper(read[j])) {
                    if (alignment->dp[i][j] + match_reward >= alignment->dp[i][j+1] - gapj &&
                        alignment->dp[i][j] + match_reward >= alignment->dp[i+1][j] - gapi) {
                        alignment->dp[i+1][j+1] = alignment->dp[i][j] + match_reward;
                        alignment->dp_trace[i+1][j+1] = 'M';
                    } else if (alignment->dp[i][j+1] - gapj >= alignment->dp[i+1][j] - gapi) {
                        alignment->dp[i+1][j+1] = alignment->dp[i][j+1] - gapj;
                        alignment->dp_trace[i+1][j+1] = 'J';
                    } else {
                        alignment->dp[i+1][j+1] = alignment->dp[i+1][j] - gapi;
                        alignment->dp_trace[i+1][j+1] = 'I';
                    }
                } else {
                    if (alignment->dp[i][j] - mm_penalty >= alignment->dp[i][j+1] - gapj &&
                        alignment->dp[i][j] - mm_penalty >= alignment->dp[i+1][j] - gapi) {
                        alignment->dp[i+1][j+1] = alignment->dp[i][j] - mm_penalty;
                        alignment->dp_trace[i+1][j+1] = 'D';
                    } else if (alignment->dp[i][j+1] - gapj >= alignment->dp[i+1][j] - gapi) {
                        alignment->dp[i+1][j+1] = alignment->dp[i][j+1] - gapj;
                        alignment->dp_trace[i+1][j+1] = 'J';
                    } else {
                        alignment->dp[i+1][j+1] = alignment->dp[i+1][j] - gapi;
                        alignment->dp_trace[i+1][j+1] = 'I';
                    }
                }
            }
        }
    }

    /* Traceback in the dp array - first count edits */
    i = end;
    j = m;
    len = 0;
    edits = 0;

    while(i >= 0 && alignment->dp_trace[i][j] == 'J') {
#ifdef DEBUG_ALIGN
        printf("%c(%d)", alignment->dp_trace[i][j], alignment->dp[i][j]);
#endif
        i--;
    }
    
    while(j > 0) {
        if (alignment->dp_trace[i][j] == 'D' || 
            alignment->dp_trace[i][j] == 'M') {
#ifdef DEBUG_ALIGN
            printf("%c(%d)", alignment->dp_trace[i][j], alignment->dp[i][j]);
#endif
            len++;
            if (alignment->dp_trace[i][j] == 'D')
                edits++;
            i--;
            j--;
        } else if (alignment->dp_trace[i][j] == 'J') {
#ifdef DEBUG_ALIGN
            printf("%c(%d)", alignment->dp_trace[i][j], alignment->dp[i][j]);
#endif
            if (j > 0 && j < m && alignment->consensus[i-1] != '-') {
                len++;
                edits++;
            }
            i--;
        } else if (alignment->dp_trace[i][j] == 'I') {
#ifdef DEBUG_ALIGN
            printf("%c(%d)", alignment->dp_trace[i][j], alignment->dp[i][j]);
#endif
            if (i > 0 && i < end) {
                len++;
                edits++;
            }
            j--;
        } else {
            printf("Undefined dp_trace %d %d\n", i, j);
            alignment->ok = 0;
            return end;
        }
    }

#ifdef DEBUG_ALIGN
    printf("\nTracing done. %d %d\n", edits, len);
#endif

    if (edits > len*max_error_rate) {
#ifdef DEBUG_ALIGN
        printf("Too many errors in read %d %c\n", pos[p].read, pos[p].ori);
#endif
        alignment->ok = 0;
        return end;
    }


    /* Traceback in the dp array */
    i = end;
    j = m;
    len = 0;
    edits = 0;
    old_end = end;

    while(i >= 0 && alignment->dp_trace[i][j] == 'J') {
#ifdef DEBUG_ALIGN
        printf("%c", alignment->dp_trace[i][j]);
#endif
        i--;
    }
    
    while(j > 0) {
        if (alignment->dp_trace[i][j] == 'D' || 
            alignment->dp_trace[i][j] == 'M') {
#ifdef DEBUG_ALIGN
            printf("%c", alignment->dp_trace[i][j]);
#endif
            er[len] = read[j-1];
            if (read[j-1] == 0) {
                printf("Whooot! %d %d\n", j-1, m);
            }
            if (quals != NULL)
                erq[len] = qual[j-1];
            len++;
            if (alignment->dp_trace[i][j] == 'D')
                edits++;
            i--;
            j--;
        } else if (alignment->dp_trace[i][j] == 'J') {
#ifdef DEBUG_ALIGN
            printf("%c", alignment->dp_trace[i][j]);
#endif
            er[len] = '-';
            if (quals != NULL) {
                erq[len] = (qual[j] + qual[j-1])/2;
            }
            len++;
            if (j > 0 && j < m && alignment->consensus[i-1] != '-') {
                edits++;
            }
            i--;
        } else if (alignment->dp_trace[i][j] == 'I') {
#ifdef DEBUG_ALIGN
            printf("%c", alignment->dp_trace[i][j]);
#endif
            er[len] = read[j-1];
            if (read[j-1] == 0) {
                printf("Whooot? %d %d\n", j-1, m);
            }
            if (quals != NULL)
                erq[len] = qual[j-1];
            len++;
            if (i > 0 && i < old_end) {
                edits++;
            }
            if (i == 0) {
                alignment->offset--;
            }
            j--;

            /* Insertion in the read. Make space for the insertion in consensus*/
            for(k = end; k >= i; k--) {
                alignment->contigA[k+1] = alignment->contigA[k];
                alignment->contigC[k+1] = alignment->contigC[k];
                alignment->contigG[k+1] = alignment->contigG[k];
                alignment->contigT[k+1] = alignment->contigT[k];
                alignment->contigN[k+1] = alignment->contigN[k];
            }
            alignment->contigA[i] = 0;
            alignment->contigC[i] = 0;
            alignment->contigG[i] = 0;
            alignment->contigT[i] = 0;
            alignment->contigN[i] = 0;
            end++;
            /* add the insertion to all aligned reads */
            for(k = 0; k < p; k++) {
                if (pos[k].aligned) {
                    if (pos[k].pos >= i) {
                        pos[k].pos++;
                    } else if (pos[k].pos < i && pos[k].pos + (signed) strlen(pos[k].edited_read) > i) {
                        alignment->contigN[i]++;
                        m2 = strlen(pos[k].edited_read);
                        if (pos[k].ori == 'C') {
                            for(kk = m2+1; kk >= 0; kk--) {
                                if (pos[k].pos + m2 - kk < i) {
                                    pos[k].edited_read[kk] = 
                                        pos[k].edited_read[kk-1];
                                    if (quals != NULL)
                                        pos[k].edited_qual[kk] = 
                                            pos[k].edited_qual[kk-1];
                                } else if (pos[k].pos + m2 - kk == i) {
                                    pos[k].edited_read[kk] = '-';
                                    if (quals != NULL) {
                                        pos[k].edited_qual[kk] = 
                                            (pos[k].edited_qual[kk-1] + pos[k].edited_qual[kk+1])/2;
                                }
                                }
                            }
                        } else {
                            for(kk = m2+1; kk >= 0; kk--) {
                                if (pos[k].pos + kk > i) {
                                    pos[k].edited_read[kk] = 
                                        pos[k].edited_read[kk-1];
                                    if (quals != NULL)
                                        pos[k].edited_qual[kk] = 
                                            pos[k].edited_qual[kk-1];
                                } else if (pos[k].pos + kk == i) {
                                    pos[k].edited_read[kk] = '-';
                                    if (quals != NULL) {
                                        pos[k].edited_qual[kk] = 
                                            (pos[k].edited_qual[kk-1] + pos[k].edited_qual[kk+1])/2;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } else {
            printf("Undefined dp_trace %d %d\n", i, j);
            alignment->ok = 0;
            return end;
        }
    }

#ifdef DEBUG_ALIGN
    printf("\nTracing done.\n");
#endif

    pos[p].pos = i;
    i--;

#ifdef DEBUG_ALIGN
    printf("Read position: %d -> %d\n", i, i+len);
    er[len] = '\0';
    printf("ER: %s\n", er);
#endif

    /* Update alignment counts */
    for(j = 0; j < len; j++) {
        switch(er[j]) {
        case 'A':
        case 'a':
            alignment->contigA[i+len-j]++;
        break;
        case 'C':
        case 'c':
            alignment->contigC[i+len-j]++;
        break;
        case 'G':
        case 'g':
            alignment->contigG[i+len-j]++;
        break;
        case 'T':
        case 't':
            alignment->contigT[i+len-j]++;
        break;
        case 'N':
        case 'n':
        case '-':
             alignment->contigN[i+len-j]++;
        break;
        default:
            printf("Trash in er %d %d/%d\n", er[j], j, len);
        }
    }

    if (pos[p].ori == 'C') {
        for(j = 0; j < len; j++) {
            if (quals != NULL)
                pos[p].edited_qual[j] = erq[j];
            switch(er[j]) {
            case 'A':
            case 'a':
                pos[p].edited_read[j] = 'T';
            break;
            case 'C':
            case 'c':
                pos[p].edited_read[j] = 'G';
            break;
            case 'G':
            case 'g':
                pos[p].edited_read[j] = 'C';
            break;
            case 'T':
            case 't':
                pos[p].edited_read[j] = 'A';
            break;
            default:
                pos[p].edited_read[j] = er[j];
                break;
            }
        }
        pos[p].edited_read[len] = '\0';
        if (quals != NULL)
            pos[p].edited_qual[len] = '\0';
    } else {
        for(j = 0; j < len; j++) {
            pos[p].edited_read[j] = er[len-j-1];
            if (quals != NULL)
                pos[p].edited_qual[j] = erq[len-j-1];
        }
        pos[p].edited_read[len] = '\0';
        if (quals != NULL)
            pos[p].edited_qual[len] = '\0';
    }

    pos[p].aligned = 1;
    return end;
}

/**
 * Compute a multiple alignment.
 * - alignment: The alignment data structure.
 * - reads: Array of reads
 * - subset: input: read ids, orientations and approximate positions for alignment
             output: read ids, orientations, positions, and read with inserted gaps
 * - size: the number of reads to align
 * Returns the number of reads aligned. This may be smaller than size if reads are specified
 * multiple times in the input data.
 */
int multi_align(align *alignment, char **reads, char **qual,
                read_pos *subset, int size,
                double max_error_rate, int max_aligned_reads,
                int match_reward, int mm_penalty, int gap_penalty) {
    int i,j;
    int end;
    int offset;

    /* Initialize alignment */
    alignment->len = 0;
    alignment->offset = 0;
    alignment->ok = 1;

    /* Sort the reads according to read ids */
    qsort(subset, size, sizeof(read_pos), compare_read);

    /* Remove any read that appears twice */
    i = 0;
    subset[0].count = 1;
    subset[0].aligned = 0;
    for(j = 1; j < size; j++) {
        if (subset[i].read != subset[j].read) {
            i++;
            if (i != j) {
                subset[i].read = subset[j].read;
                subset[i].pos = subset[j].pos;
                subset[i].ori = subset[j].ori;
                subset[i].occ = subset[j].occ;
                subset[i].count = 1;
            }
            subset[i].aligned = 0;
        } else {
            subset[i].count++;
        }
    }
    
    size = i+1;
    
    if (size > max_aligned_reads) {
#ifdef DEBUG_ALIGN
        printf("Too many reads to align %d\n", size);
#endif

        alignment->ok = 0;
        return 0;
    }

    /* Sort the reads according to counts */
    qsort(subset, size, sizeof(read_pos), compare_count);

    offset = -subset[0].pos;
    for(i = 0; i < size; i++) {
        subset[i].pos = subset[i].pos+offset;
    }

    /* Form the alignment */
    end = 0;
    for(i = 0; i < size && alignment->ok; i++) {
        j = align_read(alignment, reads, qual, subset, end, i, max_error_rate,
                       match_reward, mm_penalty, gap_penalty);
        if (j > 0)
            end = j;

        if (!alignment->ok)
            return 0;

#ifdef DEBUG_ALIGN
        printf("Consensus ends at %d\n", end);
        printf("Consensus offset is %d\n", alignment->offset);
#endif
    }

    alignment->len = end;

    update_consensus(alignment, 0, end);

#ifdef DEBUG_ALIGN
    printf("Consensus ends at %d\n", alignment->len);
#endif

    i = 0;
    for( j = 0; j < size; j++) {
        if (subset[j].aligned) {
            if (i != j) {
                subset[i].read = subset[j].read;
                subset[i].pos = subset[j].pos;
                subset[i].ori = subset[j].ori;
                subset[i].occ = subset[j].occ;
                subset[i].count = subset[j].count;
                subset[i].aligned = subset[j].aligned;

                strcpy(subset[i].edited_read, subset[j].edited_read);
                if (qual != NULL) {
                    memcpy(subset[i].edited_qual, subset[j].edited_qual, 
                           strlen(subset[i].edited_read));
                }

                subset[i].edits = subset[j].edits;                
            }
            i++;
        }
    }

    return i;
}

/* Remove such reads from read_align that do not share s k-length
   alignment with read r. Return the number of such reads */
int kalign_share(align *alignment, read_pos *read_align, int size, int k, int r, char **reads) {
    int s, e;
    int i;
    int s2, e2;
    int sc, ec;

    int count = 0;

    if (!alignment->ok)
        return size;

    s = read_align[r].pos;
    e = s + strlen(read_align[r].edited_read);

    for(i = 0; i < size; i++) {
        s2 = read_align[i].pos;
        e2 = s2 + strlen(read_align[i].edited_read);
        sc = s < s2 ? s2:s;
        ec = e < e2 ? e:e2;
        if (ec - sc >= k) {
            if (count != i) {
                read_align[count].read = read_align[i].read;
                read_align[count].pos = read_align[i].pos;
                read_align[count].ori = read_align[i].ori;
                read_align[count].edits = read_align[i].edits;
                strcpy(read_align[count].edited_read, 
                       read_align[i].edited_read);
            }
            count++;
        }
    }
    return count;
}

/* Compute the width of the part of the alignment that is shared by
   all reads. A negative number if returned if no such area is
   found. */
int common_width(align *alignment, read_pos *read_align, int size) {
    int s,e;
    int i;

    int s2, e2; 

    s = 0;
    e = MAX_READ_LENGTH;

    for(i = 0; i < size; i++) {
        s2 = read_align[i].pos;
        e2 = s2 + strlen(read_align[i].edited_read);
        if (s2 > s)
            s = s2;
        if (e2 < e)
            e = e2;
    }

    return e-s;
}

/* Compute the quality of an alignment after the alignment has been calculated. */
double align_quality(align *alignment, read_pos *read_align, int size) {
    int i,j;
    int len;
    int c;

    uchar buf[MAX_READ_LENGTH];
    double min = 1.0;

    if (!alignment->ok)
        return 0.0;

    for(i = 0; i < size; i++) {
        if (read_align[i].ori == 'U') {
            strcpy((char *)buf, read_align[i].edited_read);
        } else {
            reverse(read_align[i].edited_read, (char *)buf);
        }

        len = 0;
        c = 0;
        for(j = 0; j < (signed int)strlen((char *)buf); j++) {
            if (alignment->consensus[read_align[i].pos+j] != '-' ||
                buf[j] != '-') {
                len++;
                if (buf[j] == alignment->consensus[read_align[i].pos+j]) {
                    c++;
                }
            }
        }
        read_align[i].edits = len-c;
        if ((double)c / (double)len < min)
            min = (double)c / (double)len;
    }

    return min;
}

/* Get the consensus sequence of an alignment */
int get_consensus(align *alignment, char **cons) {
    *cons = alignment->consensus;
    return alignment->len;
}

void print_stats() {
    printf("Quick alignments: %d, Banded alignments: %d, Full alignments: %d\n",
           quick, banded, full);
}
