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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <omp.h>

#include "reverse.h"
#include "multi-align.h"
#include "q-gram-index-ht.h"


void usage(char *prog_name) {
    printf("Usage: %s -f[q,s] <input file> -o <output file> [options]\n", prog_name);

    printf("Required parameters:\n");
    printf("-f or -fq or -fs <file> Use -f for a fasta file, -fq for standard fastq file\n");
    printf("                        and -fs for Solexa fasta file\n");
    printf("-o <file>               Output file for corrected reads. Format is fasta if\n");
    printf("                        the input file is fasta and fastq if the input file\n");
    printf("                        is fastq. The output file cannot be the same file as\n");
    printf("                        the input file.\n\n");
    /*    printf("-n <int>                The number of reads in the input file.\n\n");*/

    printf("Quick options for sequencing technologies:\n");
    printf("-454                    Equal to: -mr %d -mm %d -g %d\n",
           MATCH_REWARD_454, MM_PENALTY_454, GAP_PENALTY_454);
    printf("-illumina               Equal to: -mr %d -mm %d -g %d\n\n",
           MATCH_REWARD_ILLUMINA, MM_PENALTY_ILLUMINA, GAP_PENALTY_ILLUMINA);

    printf("Other options:\n");
    printf("-k <int>                The length of a k-mer for indexing the reads. [%d]\n", DEFAULT_Q);
    /* printf("-m <int>                Maximum number of different k-mers indexed (for memory\n"); */
    /* printf("                        allocation) [%d]\n", MAX_NUM_QGRAMS); */
    printf("-e <float>              Maximum allowed error rate of a read in the multiple\n");
    printf("                        alignment. [%.2f]\n", MAX_ERROR_RATE);
    printf("-t <float>              The minimum proportion of agreeing reads to consider\n");
    printf("                        the consensus of a multiple alignment trustworthy.\n");
    printf("                        [%.2f]\n", THRESHOLD);
    printf("-q <float>              A threshold for proportion of quality values for bases\n");
    printf("                        agreeing with the consensus to consider the consensus\n");
    printf("                        trustworthy. Not used if input is fasta. [%.2f]\n", QUALITY_THRESHOLD);
    printf("-cq [<int>]             If the integer parameter is not present, old quality\n");
    printf("                        scores are retained for corrected mismatches. If the\n");
    printf("                        integer parameter is present, the quality scores of\n");
    printf("                        all corrected bases are set to that value. By default\n");
    printf("                        quality scores of corrected bases are computed as\n");
    printf("                        explained in the paper.\n");
    printf("-a <int>                Maximum number of reads to compute multiple alignments\n");
    printf("                        for. [%d]\n", MAX_ALIGNED_READS);
    printf("-p <int>                Number of threads to use. [%d]\n", NUM_THREADS);
    printf("-r <int>                The number of times the k-mer index is built during a\n");
    printf("                        correction run.\n");
    printf("-i <int>                Index only the first <int> reads. By default all\n");
    printf("                        reads are indexed.\n");
    printf("-j <int>                Do not correct the first <int> reads. By default all\n");
    printf("                        reads are corrected.\n");
    printf("-s <file>               Write statistics of computed alignments to a file. By\n");
    printf("                        default statistics are not written.\n");
    printf("-c <file>               Write the consensus sequences of computed alignments\n");
    printf("                        to a file. By default the consensus sequences are not\n");
    printf("                        written.\n");
    printf("-mr <int>               Reward for matching bases when computing alignments.\n");
    printf("                        [%d]\n", MATCH_REWARD);
    printf("-mm <int>               Penalty for mismatches when computing alignments. [%d]\n", MM_PENALTY);
    printf("-g <int>                Gap penalty when computing alignments. If gap penalty\n");
    printf("                        is higher than 100 times mismatch penalty, it is\n");
    printf("                        assumed that no gaps will occur. [%d]\n", GAP_PENALTY);
}

#define BUF_SIZE 65536
#define MAX_DUPL_READS 10000

void correct_errors(uchar **reads, int *read_length, char **qual, 
                    int64_t num_reads, int q, int max_grams, int p, uchar *num_edits, 
                    double max_error_rate, double threshold, double quality_threshold, 
                    int max_aligned_reads, int match_reward, int mm_penalty, int gap_penalty,
                    char *stat_file, char *consensus_file, 
                    int corrected_quality_scheme, int corrected_quality, int rebuilds,
                    int64_t indexed_reads, int64_t min_correct_id) {
    char *line, *line2, *lineq, *lineq2;

    int64_t i, j, k;
    int64_t ii, jj;
    FILE *f;
    FILE *sf;

    int tid; /* Thread id of the current thread */
    int dist; /* Edit distance between the original and corrected read */

    read_pos *align_pos; /* Arguments to alignment calculation */
    char *pos;
    char *consensus;     /* Consensus sequence returned by alignment calculation */
    int clen;            /* Length of consensus */

    align *alignment;    /* Data structures for alignment calculation */

    double quality;      /* Quality of an alignment */
    int share;           /* Number of reads that share a k-mer with the base read */

    int64_t num_good = 0;    /* Number of alignments that were used for correction */
    int64_t num_bad = 0;     /* Number of alignments that were not good enough for correction */
    int64_t num_perfect = 0; /* Number of alignments with no disagreeing columns */

    char *used_reads;    /* For recording which reads have been used to build a consensus */

    int bits;
    uint64_t mask;  /* Bit mask for masking the succinct representations */

    uint64_t gram;           /* Succinct representation of a q-gram */
    uint64_t gram_reverse;   /* Succinct representation of the reverse of a q-gram */
    uint64_t current;
    int lastN;

    index_t *index;
    int iter;

    char *aligned_reads;
    char *corrected_reads;

    int part;

    for(i = 0; i < num_reads; i++) {
        num_edits[i] = 0;
    }

    aligned_reads = (char *)malloc(num_reads*sizeof(char));
    corrected_reads = (char *)malloc(num_reads*sizeof(char));
    if (aligned_reads == NULL || corrected_reads == NULL) {
        printf("Malloc failed\n");
        exit(1);
    }
    for(i = 0; i < num_reads; i++) {
        aligned_reads[i] = 0;
        corrected_reads[i] = 0;
    }

    if (consensus_file != NULL) {
        f = fopen(consensus_file, "w");
        if (!f) {
            printf("Could not open consensus file %s\n", consensus_file);
            exit(1);
        }
        used_reads = (char *)malloc(num_reads*sizeof(char));
        if (used_reads == NULL) {
            printf("Malloc failed\n");
            exit(1);
        }
        for(i = 0; i < num_reads; i++)
            used_reads[i] = 0;
    }

    if (stat_file != NULL) {
        sf = fopen(stat_file, "w");
        if (!sf) {
            printf("Could not open statistics fille %s\n", stat_file);
            exit(1);
        }
    }

    omp_set_num_threads(p);

    for(iter = 0; iter < rebuilds; iter++) {
    if (iter > 0)
        free_index(index);

    i = build_index(reads, indexed_reads, q, (void **)&index, (uchar *)PREFIX, max_grams, 
                    1/*(int)(1.0/(1.0-threshold))*/, max_aligned_reads);

    if (i != 0) {
        printf("Building q-gram index failed: %s\n", error_index(i));
        exit(1);
    }

    printf("Q-gram index built.\n");

    printf("Correcting reads\n");


#pragma omp parallel private(line, line2, lineq, lineq2, i, j, k, align_pos, pos, consensus, clen, ii, jj, dist, alignment, tid, quality, gram, gram_reverse, lastN, current, bits, mask, share, part)

    {

    bits = 2;
    mask = ((uint64_t)1 << (2*index->q))-1;

#pragma omp critical
    {
        /* Is malloc thread safe? */
        align_pos = (read_pos *)malloc(MAX_DUPL_READS*sizeof(read_pos));
        alignment = (align *)malloc(sizeof(align));
        line = (char *)malloc(BUF_SIZE*sizeof(char));
        line2 = (char *)malloc(BUF_SIZE*sizeof(char));
        if (qual != NULL) {
            lineq = (char *)malloc(BUF_SIZE*sizeof(char));
            lineq2 = (char *)malloc(BUF_SIZE*sizeof(char));
        } else {
            lineq = NULL;
            lineq2 = NULL;
        }
        for(i = 0; i < MAX_DUPL_READS; i++) {
            align_pos[i].edited_read = 
                (char *)malloc(MAX_READ_LENGTH*sizeof(char));
            if (qual != NULL) {
                align_pos[i].edited_qual = 
                    (char *)malloc(MAX_READ_LENGTH*sizeof(char));
            }
        }

        for(i = 0; i < MAX_CONTIG_LEN+10; i++) {
            alignment->dp[i][0] = 0;
            alignment->dp_trace[i][0] = 'J';
        }
        for(j = 0; j < MAX_CONTIG_LEN+10; j++) {
            alignment->dp[0][j] = 0;
            alignment->dp_trace[0][j] = 'I';
        }

    }
    tid = omp_get_thread_num();
    p = omp_get_num_threads();
    printf("Thread %d/%d starting\n", tid, p);

    /* How to split the work between the threads:
     * - Threads should not work at the same time with ids close to
     *   each other because if so the likelihood of two threads
     *   working on the same read will be higher. 
     */
    for(i = iter*(int64_t)num_reads/rebuilds+tid*10000; i < (iter+1)*(int64_t)num_reads/rebuilds; i++) {
        if (i >= min_correct_id && num_edits[i] <= max_error_rate*strlen((char *)reads[i])) {

            gram = 0;
            gram_reverse = 0;
            lastN = -1; /* Initialize so that we will index the first q-gram only after reading q chars */

            for (part = 0; part < (signed int)strlen((char *)reads[i]); part +=50) {
                jj = 0;

                for(j = part; j < part+50 && j < (signed int)strlen((char *)reads[i]); j++) {
                    if (index->reads[i][j] == 'N' || index->reads[i][j] == 'n') {
                        lastN = j;
                    }

                    switch(index->reads[i][j]) {
                    case 'A':
                    case 'a':
                        gram = gram << bits | CODE_A;
                    gram_reverse = gram_reverse >> bits | ((uint64_t)CODE_T << ((index->q-1)*bits));
                    break;
                    case 'C':
                    case 'c':
                        gram = gram << bits | CODE_C;
                    gram_reverse = gram_reverse >> bits | ((uint64_t)CODE_G << ((index->q-1)*bits));
                    break;
                    case 'G':
                    case 'g':
                        gram = gram << bits | CODE_G;
                    gram_reverse = gram_reverse >> bits | ((uint64_t)CODE_C << ((index->q-1)*bits));
                    break;
                    case 'T':
                    case 't':
                        gram = gram << bits | CODE_T;
                    gram_reverse = gram_reverse >> bits | ((uint64_t)CODE_A << ((index->q-1)*bits));
                    break;
                    default:
                        gram = gram << bits;
                        gram_reverse = gram_reverse >> bits;
                        break;
                    }

                    gram = gram & mask;
                    gram_reverse = gram_reverse & mask;

                    if (lastN > j-index->q || gram_reverse == gram)
                        continue;

                    if (gram_reverse < gram) {
                        current = gram_reverse;
                    } else {
                        current = gram;
                    }

                    /* decipher the q-gram (forward and reverse representations, line and line2) */
                    for(k = 0; k < index->q; k++) {
                        switch((gram >> 2*(index->q-k-1)) & 0x03) {
                        case 0:
                            line[k] = 'A';
                            line2[index->q-k-1] = 'T';
                            break;
                        case 1:
                            line[k] = 'C';
                            line2[index->q-k-1] = 'G';
                            break;
                        case 2:
                            line[k] = 'G';
                            line2[index->q-k-1] = 'C';
                            break;
                        case 3:
                            line[k] = 'T';
                            line2[index->q-k-1] = 'A';
                            break;
                        }
                    }
                    line[index->q] = '\0';
                    line2[index->q] = '\0';

#pragma omp critical
                    {
                        k = (*index->gram2id)[current]-1;
                    }

                    if (k < 0)
                        continue;

                    if (i >= indexed_reads && jj < MAX_DUPL_READS) {
                        /* This read is not in the index so we add it here. */
                        align_pos[jj].read = i;
                        /* q-gram is always in forward orientation (however, because of using 
                           parts of the read as base read, the read might already be corrected
                           and then we cannot find the q-gram anymore) */
                        pos = strstr((char *)reads[i], line);
                        if (pos != NULL) {
                            align_pos[jj].ori = 'U';
                            align_pos[jj].pos = (j-index->q+1)-
                                ((long long int)pos -(long long int)reads[i]);
                            align_pos[jj].occ = 1;
                            pos = strstr(pos+1, line);
                            while(pos != NULL) {
                                align_pos[jj].occ++;
                                pos = strstr(pos+1, line);
                            }
                            jj++;
                        }
                    }

                    /* initiliaze input for alignment calculation */
                    for(ii = 0; ii < (signed int)index->counts[k] && jj < MAX_DUPL_READS; ii++) {
                        if (num_edits[index->gram_reads[k][ii]] > 
                            max_error_rate*strlen((char *)reads[index->gram_reads[k][ii]]))
                            continue;
                        align_pos[jj].read = index->gram_reads[k][ii];
                        /* Search for the q-gram in forward orientation */
                        pos = strstr((char *)reads[index->gram_reads[k][ii]], line);
                        if (pos != NULL) {
                            align_pos[jj].ori = 'U';
                            align_pos[jj].pos = (j-index->q+1)-
                                ((long long int)pos -(long long int)reads[index->gram_reads[k][ii]]);
                            align_pos[jj].occ = 1;
                            while(ii + 1 < (signed int)index->counts[k] &&
                                  index->gram_reads[k][ii+1] == index->gram_reads[k][ii]) {
                                ii++;
                                align_pos[jj].occ++;
                            }
                            jj++;
                        } else {
                            /* Forward q-gram not found, search for the reverse complement */
                            pos = strstr((char *)reads[index->gram_reads[k][ii]], line2);
                            if (pos != NULL) {
                                align_pos[jj].ori = 'C';
                                align_pos[jj].pos =  (j-index->q+1) 
                                    - (strlen((const char *)reads[index->gram_reads[k][ii]]) -
                                       ((int64_t)pos - (int64_t)reads[index->gram_reads[k][ii]]) -
                                       index->q);
                                align_pos[jj].occ = 1;
                                while(ii + 1 < (signed int)index->counts[k] &&
                                      index->gram_reads[k][ii+1] == index->gram_reads[k][ii]) {
                                    ii++;
                                    align_pos[jj].occ++;
                                }
                                jj++;
                            }
                        }
                    }
                }
        
                j = jj;

                if (j >= 1.0/(1.0-threshold) && j < MAX_DUPL_READS) {

#pragma omp critical
                    {
                        for(ii = 0; ii < j; ii++) {
                            aligned_reads[align_pos[ii].read] = 1;
                        }
                    }

                    /* Compute multiple alignment */
                    j = multi_align(alignment, (char **)reads, qual, align_pos, j, 
                                    max_error_rate, max_aligned_reads,
                                    match_reward, mm_penalty, gap_penalty);

                    for(k = 0; k < j; k++)
                        if (align_pos[k].read == (ulong) i)
                            break;

                    clen = get_consensus(alignment, &consensus);
                    consensus[clen] = '\0';
            
                    if (k >= j) {
                        quality = 0.0;
                        share = 0;
                    } else {
                        share = kalign_share(alignment, align_pos, j, index->q, k, 
                                             (char **)reads);
                    }

                    if (share != j && share > 0) {
                        /* Recompute multiple alignment */
                        j = multi_align(alignment, (char **)reads, qual, align_pos, 
                                        share, max_error_rate, max_aligned_reads,
                                        match_reward, mm_penalty, gap_penalty);

                        for(k = 0; k < j; k++)
                            if (align_pos[k].read == (ulong) i)
                                break;

                        clen = get_consensus(alignment, &consensus);
                        consensus[clen] = '\0';
                
                        if (k >= j) {
                            quality = 0.0;
                            share = 0;
                        } else {
                            share = kalign_share(alignment, align_pos, j, index->q, k,
                                                 (char **)reads);
                        }
                    }

                    if (share == j) {
                        quality = align_quality(alignment, align_pos, j);
                    } else {
                        quality = 0.0;
                    }
#pragma omp critical
                    {
                        if (stat_file != NULL) {
                            snprintf(line2, BUF_SIZE, "%ld\t%f\t%d\t%ld\n", i,
                                     1.0-quality, share, j);
                            if (!fwrite(line2, sizeof(char), strlen(line2), sf)) {
                                printf("Could not write to statistics file\n");
                            }
                        }


                        /* update statistics */
                        if (quality >= 1.0) {
                            num_perfect++;
                        } else if (quality > 1.0-max_error_rate) {
                            num_good++;
                        } else {
                            num_bad++;
                        }
                    }


#ifdef DEBUG
#pragma omp critical
                    {
                        /* if (quality > 1.0-max_error_rate && quality < 1.0 && share == j) { */
                        if (quality < 1.0) {
                            printf("******************\n");
                            printf("%d: %s\n", i, reads[i]);
                    
                            printf("  %s\n", consensus);
                    
                            for(k = 0; k < j; k++) {
                                printf("%c ", align_pos[k].ori);
                                for(ii = 0; ii < align_pos[k].pos; ii++) {
                                    printf(" ");
                                }
                                if (align_pos[k].ori == 'C') {
                                    reverse(align_pos[k].edited_read, line2);
                                    printf("%s\n", line2);
                                } else {
                                    printf("%s\n", align_pos[k].edited_read);
                                }
                            }

                            if (qual != NULL) {
                                for(k = 0; k < j; k++) {
                                    printf("%c ", align_pos[k].ori);
                                    for(ii = 0; ii < align_pos[k].pos; ii++) {
                                        printf("   ");
                                    }
                                    if (align_pos[k].ori == 'C') {
                                        for(ii = strlen(align_pos[k].edited_read)-1; ii >= 0; ii--) {
                                            printf("%d ", align_pos[k].edited_qual[ii]);
                                        }
                                        printf("\n");
                                    } else {
                                        for(ii = 0; ii < (signed)strlen(align_pos[k].edited_read); ii++) {
                                            printf("%d ", align_pos[k].edited_qual[ii]);
                                        }
                                        printf("\n");
                                    }
                                }
                            }
                        }
                        printf("OK: %d, Quality: %f Share: %d/%d\n", alignment->ok, quality, share, j);
                    }
#endif

                    if (quality > 1.0-max_error_rate && share == j) {           
#pragma omp critical
                        {
                            for(ii = 0; ii < j; ii++) {
                                corrected_reads[align_pos[ii].read] = 1;
                            }
                        }
                    }


                    if (quality > 1.0-max_error_rate && quality < 1.0 && share == j) {
                        /* Alignment has an appropriate quality for correction */

                        clen = get_consensus(alignment, &consensus);

                        /* Quality values for the consensus */
                        int cons_qual[2*MAX_READ_LENGTH+1];
                        int tot_qual[2*MAX_READ_LENGTH+1];

                        if (qual != NULL) {
                            for(k = 0; k < clen; k++) {
                                cons_qual[k] = 0;
                                tot_qual[k] = 0;
                            }

                            for(k = 0; k < j; k++) {
                                /* Copy to line2 a representation of the
                                   edited read in the same orientation as the consensus */
                                if (align_pos[k].ori == 'U') {
                                    strcpy(line2, align_pos[k].edited_read);
                                    if (qual != NULL)
                                        memcpy(lineq2, align_pos[k].edited_qual, 
                                               strlen(line2));
                                } else {
                                    reverse(align_pos[k].edited_read, line2);
                                    if (qual != NULL) {
                                        for(jj = 0; jj < (signed)strlen(line2); jj++)
                                            lineq2[strlen(line2)-jj-1] = 
                                                align_pos[k].edited_qual[jj];
                                        lineq2[strlen(line2)] = '\0';
                                    }
                                }

                                for(ii = 0; ii < (signed int)strlen(line2); ii++) {
                                    if (consensus[align_pos[k].pos + ii] == line2[ii]) {
                                        cons_qual[align_pos[k].pos + ii] += lineq2[ii];
                                    }
                                    tot_qual[align_pos[k].pos + ii] += lineq2[ii];
                                }
                            }
                        }
                
                        for(k = 0; k < j; k++) {
                            if ((signed)align_pos[k].read <= min_correct_id) {
                                continue;
                            }

                            if (align_pos[k].edits > 0) {
                        
                                /* Copy to line2 a representation of the
                                   edited read in the same orientation as the consensus */
                                if (align_pos[k].ori == 'U') {
                                    strcpy(line2, align_pos[k].edited_read);
                                    if (qual != NULL)
                                        memcpy(lineq2, align_pos[k].edited_qual, 
                                               strlen(line2));
                                } else {
                                    reverse(align_pos[k].edited_read, line2);
                                    if (qual != NULL) {
                                        for(jj = 0; jj < (signed)strlen(line2); jj++)
                                            lineq2[strlen(line2)-jj-1] = 
                                                align_pos[k].edited_qual[jj];
                                        lineq2[strlen(line2)] = '\0';
                                    }
                                }
                    
#ifdef DEBUG
                                printf("Original  read %s\n", reads[align_pos[k].read]);
#endif
                                jj = 0;
                                dist = 0;
                                for(ii = 0; ii < (signed int)strlen(line2); ii++) {
                                    int consC; /* Number of reads supporting consensus at this position */
                                    //int editC; /* Number of reads supporting this read at this position */
                        
                                    switch(consensus[align_pos[k].pos+ii]) {
                                    case 'A':
                                    case 'a':
                                        consC = alignment->contigA[align_pos[k].pos+ii];
                                    break;
                                    case 'C':
                                    case 'c':
                                        consC = alignment->contigC[align_pos[k].pos+ii];
                                    break;
                                    case 'G':
                                    case 'g':
                                        consC = alignment->contigG[align_pos[k].pos+ii];
                                    break;
                                    case 'T':
                                    case 't':
                                        consC = alignment->contigT[align_pos[k].pos+ii];
                                    break;
                                    case 'N':
                                    case 'n':
                                    case '-':
                                        consC = alignment->contigN[align_pos[k].pos+ii];
                                    break;
                                    default:
                                        consC = 0;
                                        break;
                                    }
                        
                                    /* switch(line2[ii]) { */
                                    /* case 'A': */
                                    /* case 'a': */
                                    /*     editC = alignment->contigA[align_pos[k].pos+ii]; */
                                    /* break; */
                                    /* case 'C': */
                                    /* case 'c': */
                                    /*     editC = alignment->contigC[align_pos[k].pos+ii]; */
                                    /* break; */
                                    /* case 'G': */
                                    /* case 'g': */
                                    /*     editC = alignment->contigG[align_pos[k].pos+ii]; */
                                    /* break; */
                                    /* case 'T': */
                                    /* case 't': */
                                    /*     editC = alignment->contigT[align_pos[k].pos+ii]; */
                                    /* break; */
                                    /* case 'N': */
                                    /* case 'n': */
                                    /* case '-': */
                                    /*     editC = alignment->contigN[align_pos[k].pos+ii]; */
                                    /* break; */
                                    /* default: */
                                    /*     editC = 0; */
                                    /*     break; */
                                    /* } */
                        
                                    if ((double)consC / (double)(alignment->contigA[align_pos[k].pos+ii]+
                                                                 alignment->contigC[align_pos[k].pos+ii]+
                                                                 alignment->contigG[align_pos[k].pos+ii]+
                                                                 alignment->contigT[align_pos[k].pos+ii]+
                                                                 alignment->contigN[align_pos[k].pos+ii]) >= threshold
                                        && 
                                        (qual == NULL? 1 : 
                                         (double)cons_qual[align_pos[k].pos+ii]/(double)tot_qual[align_pos[k].pos+ii]
                                         >= quality_threshold)
                                        /*cons_qual[align_pos[k].pos+ii] - lineq2[ii] > quality_threshold)*/) {
                                        if (line2[ii] != consensus[align_pos[k].pos+ii])
                                            dist++;
                                        if (consensus[align_pos[k].pos+ii] != '-') {
                                            line[jj] = consensus[align_pos[k].pos+ii];
                                            if (qual != NULL) {
                                                if (line2[ii] == 
                                                    consensus[align_pos[k].pos+ii]) {
                                                    lineq[jj] = lineq2[ii];
                                                } else {
                                                    if (corrected_quality_scheme == FLAT) {
                                                        lineq[jj] = corrected_quality;
                                                    } else if (corrected_quality_scheme == ESTIMATE ||
                                                               line2[ii] == '-') {
                                                        lineq[jj] = cons_qual[align_pos[k].pos+ii]/consC;
                                                    } else {
                                                        lineq[jj] = lineq2[ii];
                                                    }
                                                }
                                            }
                                            jj++;
                                        }
                                    } else {
                                        if (line2[ii] != '-') {
                                            line[jj] = line2[ii];
                                            if (qual != NULL)
                                                lineq[jj] = lineq2[ii];
                                            jj++;
                                        }
                                    }
                                }
                                line[jj] = '\0';
                                if (qual != NULL)
                                    lineq[jj] = '\0';
                        
                                // Check that the corrected read fits into the array
                                if ((int)strlen(line) <= read_length[align_pos[k].read] + MAX_INSERTIONS) {
                                    /* Copy the original read in the original orientation to the reads array */
#pragma omp critical
                                    {
                                    if (align_pos[k].ori == 'U') {
                                        strcpy((char *)reads[align_pos[k].read], line);
                                        if (qual != NULL)
                                            memcpy(qual[align_pos[k].read], lineq, 
                                                   strlen(line));
                                    } else {
                                        reverse(line, (char *)reads[align_pos[k].read]);
                                        if (qual != NULL) {
                                            for(jj = 0; jj < (signed)strlen(line); jj++) {
                                                qual[align_pos[k].read][strlen(line)-jj-1] = lineq[jj];
                                            }
                                            qual[align_pos[k].read][strlen(line)] = '\0';
                                        }
                                    }
                        
#ifdef DEBUG
                                    printf("Corrected read %s\n", reads[align_pos[k].read]);
#endif
                                    num_edits[align_pos[k].read] += dist;
                                    }
                                }
                            }
                        }
                    } else {
                        j = share;
                    }
            
                    if (quality > 1.0-max_error_rate) {
#pragma omp critical
                        {
                            /* Write consensus sequence to file */
                            if (consensus_file != NULL) {
                                for(ii = 0, jj = 0; jj < (int64_t)strlen(consensus); jj++) {
                                    if (consensus[jj] != '-') {
                                        consensus[ii] = consensus[jj];
                                        ii++;
                                    }
                                }
                                consensus[ii] = '\0';
                        
                                snprintf(line, BUF_SIZE, ">consensus_%ld (%f)\n%s\n", i, quality, consensus);
                                if (!fwrite(line, sizeof(char), strlen(line), f)) {
                                    printf("Writing to consensus file failed\n");
                                }
                        
                                for(k = 0; k < j; k++) {
                                    used_reads[align_pos[k].read] = 1;
                                }
                            }
                        }
                    }
                }
            }
        }

        /* Print some hint of progress */
        if (i % 10000 == 9999) {
            printf("Thread %d: %ld/%ld Perfect alignments: %ld, Good alignments: %ld, Bad alignments: %ld\n", 
                   tid, i, num_reads, num_perfect, num_good, num_bad);
            /* printf("Thread %d: %d/%d\n", tid, i, num_reads); */
            i += (p-1)*10000;
        }
    }

#pragma omp critical
    {
        for(i = 0; i < MAX_DUPL_READS; i++) {
            free(align_pos[i].edited_read);
            if (qual != NULL)
                free(align_pos[i].edited_qual);
        }
        free(align_pos);
        free(alignment);
        free(line);
        free(line2);
        if (qual != NULL) {
            free(lineq);
            free(lineq2);
        }
        printf("Thread %d/%d finishing\n", tid, p);
    }

    } // end omp parallel

    }

    printf("Perfect alignments: %ld, Good alignments: %ld, Bad alignments: %ld\n", 
           num_perfect, num_good, num_bad);

    int num_aligned = 0;
    int num_corrected = 0;
    for(i = 0; i < num_reads; i++) {
        if (aligned_reads[i] > 0)
            num_aligned++;
        if (corrected_reads[i] > 0)
            num_corrected++;
    }

    printf("%d reads were aligned, %d reads were aligned in a good alignment\n",
           num_aligned, num_corrected);

    if (stat_file != NULL)
        fclose(sf);

    line = (char *)malloc(BUF_SIZE*sizeof(char));

    if (consensus_file != NULL) {
        for(i = 0; i < num_reads; i++) {
            if (!used_reads[i]) {
                snprintf(line, BUF_SIZE, ">unused_read_%ld\n%s\n", i, reads[i]);
                if (!fwrite(line, sizeof(char), strlen(line), f)) {
                    printf("Could not write to consensus file\n");
                }
            }
        }
        fclose(f);
    }
}


/**
 * Main routine of error correction.
 * Command line options:
 * -f[q] <input file>
 * -o <output file>
 * -n <number of reads>
 * -k <k-mer length>
 * -m <maximum number of different k-mers>
 * -e <max error rate for alignment>
 * -t <minimum proportion of agreeing reads to consider the consensus trustworthy>
 * -q <threshold for the sum of quality values of bases agreeing with consensus minus 
 *    the quality of the base to correct>
 * -a <maximum number of reads to compute multiple alignments for>
 * -p <number of threads to use>
 * -pr <prefix for considered q-grams> (flag not implemented)
 * -s <statistics file>
 * -c <consensus file>
 * -r <number of rebuilds>
 */
int main(int argc, char *argv[]) {
    char *fasta_file = NULL;
    char *fastq_file = NULL;
    int64_t num_reads = 0;
    char *out_file = NULL;
    char *stat_file = NULL;
    char *consensus_file = NULL;

    uchar **reads;
    char **qual;
    int *read_length;

    char *line = (char *)malloc(BUF_SIZE*sizeof(char));
    char *line2 = (char *)malloc(BUF_SIZE*sizeof(char));

    int64_t i,k,j,jj;
    int64_t ii;

    FILE *f;
    FILE *of;
    uchar *buf4reads = NULL;
    char *buf4qual = NULL;

    int c;

    int p = NUM_THREADS; /* Number of threads to use */

    /* Number of edit operations performed on each read */
    uchar *num_edits;

    int q = DEFAULT_Q;   /* The length of a q-gram */

    char *prog_name = argv[0];
    uint64_t max_grams = MAX_NUM_QGRAMS;
    double max_error_rate = MAX_ERROR_RATE;
    double threshold = THRESHOLD;
    double quality_threshold = QUALITY_THRESHOLD;
    int max_aligned_reads = MAX_ALIGNED_READS;

    int match_reward = MATCH_REWARD;
    int mm_penalty = MM_PENALTY;
    int gap_penalty = GAP_PENALTY;

    int corrected_quality_scheme = ESTIMATE;
    int corrected_quality = 0;

    int quality_offset = 64;

    int rebuilds = 1;
    int64_t indexed_reads = -1;
    int64_t min_correct_id = -1;

    while(argc > 0) {
        if (!strcmp(argv[0], "-f")) {
            if (argc > 1) {
                fasta_file = argv[1];
                fastq_file = NULL;
            } else {
                usage(prog_name);
                return 1;
            }
            argc--;
            argv++;
        } else if (!strcmp(argv[0], "-fq")) {
            if (argc > 1) {
                fastq_file = argv[1];
                fasta_file = NULL;
                quality_offset = 33;
            } else {
                usage(prog_name);
                return 1;
            }
            argc--;
            argv++;
        } else if (!strcmp(argv[0], "-fs")) {
            if (argc > 1) {
                fastq_file = argv[1];
                fasta_file = NULL;
                quality_offset = 64;
            } else {
                usage(prog_name);
                return 1;
            }
            argc--;
            argv++;
        } else if (!strcmp(argv[0], "-o")) {
            if (argc > 1) {
                out_file = argv[1];
            } else {
                usage(prog_name);
                return 1;
            }
            argc--;
            argv++;
        /* } else if (!strcmp(argv[0], "-n")) { */
        /*     if (argc > 1) { */
        /*         num_reads = atoi(argv[1]); */
        /*     } else { */
        /*         usage(prog_name); */
        /*         return 1; */
        /*     } */
        /*     argc--; */
        /*     argv++; */
        } else if (!strcmp(argv[0], "-k")) {
            if (argc > 1) {
                q = atoi(argv[1]);
            } else {
                usage(prog_name);
                return 1;
            }
            argc--;
            argv++;
        } else if (!strcmp(argv[0], "-p")) {
            if (argc > 1) {
                p = atoi(argv[1]);
            } else {
                usage(prog_name);
                return 1;
            }
            argc--;
            argv++;
        /* } else if (!strcmp(argv[0], "-m")) { */
        /*     if (argc > 1) { */
        /*         max_grams = atol(argv[1]); */
        /*     } else { */
        /*         usage(prog_name); */
        /*         return 1; */
        /*     } */
        /*     argc--; */
        /*     argv++; */
        } else if (!strcmp(argv[0], "-e")) {
            if (argc > 1) {
                max_error_rate = atof(argv[1]);
            } else {
                usage(prog_name);
                return 1;
            }
            argc--;
            argv++;
        } else if (!strcmp(argv[0], "-t")) {
            if (argc > 1) {
                threshold = atof(argv[1]);
            } else {
                usage(prog_name);
                return 1;
            }
            argc--;
            argv++;
        } else if (!strcmp(argv[0], "-q")) {
            if (argc > 1) {
                quality_threshold = atof(argv[1]);
            } else {
                usage(prog_name);
                return 1;
            }
            argc--;
            argv++;
        } else if (!strcmp(argv[0], "-cq")) {
            if (argc > 1 && argv[1][0] != '-') {
                corrected_quality_scheme = FLAT;
                corrected_quality = atoi(argv[1]);
                argc--;
                argv++;
            } else {
                corrected_quality_scheme = KEEP;
            }
        } else if (!strcmp(argv[0], "-a")) {
            if (argc > 1) {
                max_aligned_reads = atoi(argv[1]);
            } else {
                usage(prog_name);
                return 1;
            }
            argc--;
            argv++;
        } else if (!strcmp(argv[0], "-c")) {
            if (argc > 1) {
                consensus_file = argv[1];
            } else {
                usage(prog_name);
                return 1;
            }
            argc--;
            argv++;
        } else if (!strcmp(argv[0], "-s")) {
            if (argc > 1) {
                stat_file = argv[1];
            } else {
                usage(prog_name);
                return 1;
            }
            argc--;
            argv++;
        } else if (!strcmp(argv[0], "-mr")) {
            if (argc > 1) {
                match_reward = atoi(argv[1]);
            } else {
                usage(prog_name);
                return 1;
            }
            argc--;
            argv++;
        } else if (!strcmp(argv[0], "-mm")) {
            if (argc > 1) {
                mm_penalty = atoi(argv[1]);
            } else {
                usage(prog_name);
                return 1;
            }
            argc--;
            argv++;
        } else if (!strcmp(argv[0], "-g")) {
            if (argc > 1) {
                gap_penalty = atoi(argv[1]);
            } else {
                usage(prog_name);
                return 1;
            }
            argc--;
            argv++;
        } else if (!strcmp(argv[0], "-r")) {
            if (argc > 1) {
                rebuilds = atoi(argv[1]);
            } else {
                usage(prog_name);
                return 1;
            }
            argc--;
            argv++;
        } else if (!strcmp(argv[0], "-i")) {
            if (argc > 1) {
                indexed_reads = atoi(argv[1]);
            } else {
                usage(prog_name);
                return 1;
            }
            argc--;
            argv++;
        } else if (!strcmp(argv[0], "-j")) {
            if (argc > 1) {
                min_correct_id = atoi(argv[1]);
            } else {
                usage(prog_name);
                return 1;
            }
            argc--;
            argv++;
        } else if (!strcmp(argv[0], "-454")) {
            match_reward = MATCH_REWARD_454;
            mm_penalty = MM_PENALTY_454;
            gap_penalty = GAP_PENALTY_454;
        } else if (!strcmp(argv[0], "-illumina")) {
            match_reward = MATCH_REWARD_ILLUMINA;
            mm_penalty = MM_PENALTY_ILLUMINA;
            gap_penalty = GAP_PENALTY_ILLUMINA;
        }
        argc--;
        argv++;
    }
    
    if (out_file == NULL/* || num_reads == 0*/) {
        usage(prog_name);
        return 1;
    }

    if ((fasta_file != NULL && strcmp(out_file, fasta_file) == 0)||
	(fastq_file != NULL && strcmp(out_file, fastq_file) == 0)) {
	printf("The input and output files must not be the same file!\n");
	usage(prog_name);
	return 1;
    }
    i = 0;
    k = BUF_SIZE;

    if (fasta_file != NULL) {
        qual = NULL;
        printf("Counting reads in file %s\n", fasta_file);
        f = fopen(fasta_file, "r");
        if (f == NULL) {
            printf("Could not open file %s\n", fasta_file);
            abort();
        }
        c = fgetc(f);
        while(c != EOF) {
            if (c != '>') {
                printf("Fasta file %s has invalid format\n", fasta_file);
                abort();
            }

            num_reads++;

            while((c = fgetc(f)) != '\n' && c != EOF);
            while((c = fgetc(f)) != '>' && c != EOF);
        }

        printf("Found %ld reads\n", num_reads);
        if (indexed_reads < 0 || indexed_reads > num_reads)
            indexed_reads = num_reads;
        reads = (uchar **)malloc(num_reads * sizeof(uchar *));
        read_length = (int *)malloc(num_reads * sizeof(int));
        num_edits = (uchar *)malloc(num_reads * sizeof(uchar));

        printf("Reading reads from file %s\n", fasta_file);

        /* Read the reads from the file */
        f = fopen(fasta_file, "r");
        if (f == NULL) {
            printf("Could not open file %s\n", fasta_file);
            abort();
        }
        c = fgetc(f);
        while(c != EOF) {
            j = 0;
            line[j++] = c;
            while((c = fgetc(f)) != '\n' && c != EOF && j < BUF_SIZE) {
                line[j++] = c;
            }
            line[j] = '\0';

            if (j >= BUF_SIZE || c == EOF || line[0] != '>') {
                printf("Fasta file %s has invalid format\n", fasta_file);
                abort();
            }

            j = 0;
            while((c = fgetc(f)) != '>' && c != EOF && j < BUF_SIZE) {
                if (c != '\n')
                    line[j++] = c;
            }
            line[j] = '\0';

            if (j >= BUF_SIZE) {
                printf("Read length exceeds the maximum.\n");
                abort();
            }

            if (j+MAX_INSERTIONS+1+k > BUF_SIZE) {
                buf4reads = (uchar *)malloc(BUF_SIZE*sizeof(uchar));
                if (buf4reads == NULL) {
                    printf("Could not malloc space for reads\n");
                    abort();
                }
                k = 0;
            }
            
            strncpy((char *)&buf4reads[k], line, j+1);
            reads[i] = &buf4reads[k];
            read_length[i] = j;
            
            i++;
            k += j+1;
            for(ii = 0; ii < MAX_INSERTIONS; ii++) {
                buf4reads[k++] = '\0';
            }

#ifdef DEBUG_Q
            printf("Original: %s\n", reads[i-1]);
#endif

        }
        fclose(f);

        printf("Found %ld reads\n", i);
    } else if (fastq_file != NULL) {
        printf("Counting reads in file %s\n", fastq_file);
        f = fopen(fastq_file, "r");
        if (f == NULL) {
            printf("Could not open file %s\n", fastq_file);
            abort();
        }
        c = fgetc(f);
        while(c != EOF) {
            if (c != '@') {
                printf("Fastq file %s has invalid format\n", fastq_file);
                abort();
            }

            num_reads++;

            /*comment*/
            while((c = fgetc(f)) != '\n' && c != EOF);
            /*read*/
            while((c = fgetc(f)) != '\n' && c != EOF);
            /*comment*/
            while((c = fgetc(f)) != '\n' && c != EOF);
            /*qualities*/
            while((c = fgetc(f)) != '\n' && c != EOF);
            c = fgetc(f);
        }

        printf("Found %ld reads\n", num_reads);
        if (indexed_reads < 0 || indexed_reads > num_reads)
            indexed_reads = num_reads;
        reads = (uchar **)malloc(num_reads * sizeof(uchar *));
        read_length = (int *)malloc(num_reads * sizeof(int));
        num_edits = (uchar *)malloc(num_reads * sizeof(uchar));
        qual = (char **)malloc(num_reads * sizeof(char *));

        printf("Reading reads from file %s\n", fastq_file);

        /* Read the reads from the file */
        f = fopen(fastq_file, "r");
        if (f == NULL) {
            printf("Could not open file %s\n", fastq_file);
            return 1;
        }
        while((c = fgetc(f)) != EOF) {
            if (c != '@') {
                printf("Fastq file %s has invalid format\n", fastq_file);
                exit(1);
            }
            /* Skip read id */
            while(c != '\n' && c != EOF) c = fgetc(f); 
            
            /* Bases */
            j = 0;
            while((c = fgetc(f)) != '\n') {if (c == EOF) break; line[j++] = c;}
            line[j] = '\0';

            /* Skip the separating line */
            if (fgetc(f) != '+') {
                printf("Fastq file %s has invalid format\n", fastq_file);
                exit(1);
            }
            while(fgetc(f) != '\n');

            /* Qualities */
            jj = 0;
            while((c = fgetc(f)) != '\n') {
                if (c == EOF) break; 
                line2[jj++] = c-quality_offset;
            }
            line2[jj] = '\0';
            if (j != jj) {
                printf("Length of bases and qualities not the same: %ld, %ld\n",
                       j, jj);
                exit(1);
            }

            
            if (j+MAX_INSERTIONS+1+k > BUF_SIZE) {
                buf4reads = (uchar *)malloc(BUF_SIZE*sizeof(uchar));
                buf4qual = (char *)malloc(BUF_SIZE*sizeof(char));
                if (buf4reads == NULL || buf4qual == NULL) {
                    printf("Could not malloc space for reads\n");
                    abort();
                }
                k = 0;
            }
            
            strncpy((char *)&buf4reads[k], line, j+1);
	    memcpy((char *)&buf4qual[k], line2, j+1);
            //strncpy((char *)&buf4qual[k], line2, j+1);
            reads[i] = &buf4reads[k];
            qual[i] = &buf4qual[k];
            read_length[i] = j;
            
            i++;
            k += j+1;
            for(ii = 0; ii < MAX_INSERTIONS; ii++) {
                buf4reads[k] = '\0';
                buf4qual[k++] = '\0';
            }

#ifdef DEBUG_Q
            printf("Original: %s\n", reads[i-1]);
            printf("          ");
            for(jj= 0; jj < j; jj++)
                printf("%d ", qual[i-1][jj]);
            printf("\n");
#endif

        }
        fclose(f);

        printf("Found %ld reads\n", i);

    } else {
        usage(prog_name);
        abort();
    }

    correct_errors(reads, read_length, qual, num_reads, q, max_grams, p, num_edits, 
                   max_error_rate, threshold, quality_threshold,
                   max_aligned_reads, match_reward, mm_penalty, gap_penalty,
                   stat_file, consensus_file, 
                   corrected_quality_scheme, corrected_quality, rebuilds, indexed_reads,
                   min_correct_id);

    print_stats();

    printf("Correction finished. Outputting corrected reads.\n");

    if (fasta_file) {

        f = fopen(out_file, "w");

        /*  Open the original file to read the read names from there */
        of = fopen(fasta_file, "r");
        if (of == NULL) {
            printf("Could not open file %s\n", fasta_file);
            abort();
        }
        c = fgetc(of);

        for(i = 0; i < num_reads; i++) {

            /* Read the comment line */
            j = 0;
            line2[j++] = c;
            while((c = fgetc(of)) != '\n' && c != EOF && j < BUF_SIZE) {
                line2[j++] = c;
            }
            line2[j] = '\0';

            if (j >= BUF_SIZE || c == EOF || line2[0] != '>') {
                printf("Fasta file %s has invalid format\n", fasta_file);
                abort();
            }

            /* Skip the read */
            while((c = fgetc(of)) != '>' && c != EOF);

            snprintf(line, BUF_SIZE, "%s ( %d edit operations)\n", line2, num_edits[i]);
            if (!fwrite(line, sizeof(char), strlen(line), f)) {
                printf("Could not write to output file\n");
            }
            snprintf(line, BUF_SIZE, "%s\n", reads[i]);
            if (!fwrite(line, sizeof(char), strlen(line), f)) {
                printf("Could not write to output file\n");
            }
        }
        fclose(f);
    } else if (fastq_file) {
        f = fopen(out_file, "w");

        /*  Open the original file to read the comments from there */
        of = fopen(fastq_file, "r");
        if (of == NULL) {
            printf("Could not open file %s\n", fastq_file);
            abort();
        }
        c = fgetc(of);

        for(i = 0; i < num_reads; i++) {

            if (c != '@') {
                printf("Fastq file %s has invalid format\n", fastq_file);
                abort();
            }

            /*comment*/
            j = 0;
            line2[j++] = c;
            while((c = fgetc(of)) != '\n' && c != EOF) {
                line2[j++] = c;
            }
            line2[j] = '\0';

            /*skip the read*/
            while((c = fgetc(of)) != '\n' && c != EOF);

            snprintf(line, BUF_SIZE, "%s ( %d edit operations)\n", line2, num_edits[i]);
            if (!fwrite(line, sizeof(char), strlen(line), f)) {
                printf("Could not write to output file\n");
            }
            snprintf(line, BUF_SIZE, "%s\n", reads[i]);
            if (!fwrite(line, sizeof(char), strlen(line), f)) {
                printf("Could not write to output file\n");
            }

            /*comment*/
            j = 0;
            while((c = fgetc(of)) != '\n' && c != EOF) {
                line2[j++] = c;
            }
            line2[j] = '\0';

            /*skip qualities*/
            while((c = fgetc(of)) != '\n' && c != EOF);
            c = fgetc(of);

            snprintf(line, BUF_SIZE, "%s\n", line2);
            if (!fwrite(line, sizeof(char), strlen(line), f)) {
                printf("Could not write to output file\n");
            }
            for(j = 0; j < (signed)strlen((char *)reads[i]); j++) {
                line[j] = qual[i][j]+quality_offset;
            }
            line[strlen((char *)reads[i])] = '\n';
            line[1+strlen((char *)reads[i])] = '\0';

            if (!fwrite(line, sizeof(char), strlen(line), f)) {
                printf("Could not write to output file\n");
            }
        }
        fclose(f);
        fclose(of);
    }


    return 0;
}
