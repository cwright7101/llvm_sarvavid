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
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <limits.h>

#include "q-gram-index-ht.h"

#define INVALID ((uint64_t)0xffffffffffffffff)

/* Error descriptions for the construction of the q-gram index */
const char *error_strings[] = {
    "",
    "Malloc failed.", 
    "Hash table too small",
    "Pattern was too short.",
    "Filename too long.",
    "File open failed.",
    "File write failed.",
    "Stating a file failed.",
    "File read failed.",
    "File has invalid format.",
    "Position list too long.",
    "Invalid build options.",
    "Too many matches.",
    "Unexpected error when building the index."
};

#ifdef STATS
// Old code that does not work with the current version of the index
void calculate_stats(index_t *index, int max_pos_list_len) {
    int *pos_list_len;
    int i;

    pos_list_len = malloc(sizeof(int)*(max_pos_list_len+1));

    for(i = 0; i <= max_pos_list_len; i++) {
        pos_list_len[i] = 0;
    }

    for(i = 0; i < index->size; i++) {
        if (i == 0) {
            pos_list_len[index->index[i].plen]++;
        } else {
            pos_list_len[index->index[i].plen- index->index[i-1].plen]++;
        }
    }

    printf("\nposition list lengths:\n");
    for(i = 0; i <= max_pos_list_len; i++) {
        if (pos_list_len[i]>0)
            printf("\t%d:\t%d\n", i, pos_list_len[i]);
    }
    free(pos_list_len);
}
#endif

/**
 * Construct the q-gram index for a set of reads
 * - reads: the set of reads
 * - num_reads: number of reads in the set
 * - q: the length of a q-gram'
 * - index: a pointer to the constructed index will be returned here
 * - prefix: the index contains only q-grams starting with this prefix. Not 
 *   compatible with the new speedup techniques (can be used but some speedup 
 *   will be lost)
 * - size: The maximum number of different q-grams that can be inserted to this index
 * - min_reads_per_gram: Minimum number of reads per q-gram to keep the read list
 * - max_reads_per_gram: Maximum number of reads per q-gram to keep the read list
 * Returns 0 on success, error code otherwise.
 */
int build_index(uchar **reads, ulong num_reads, int q, void **index, uchar *prefix, uint64_t size,
                int min_reads_per_gram, int max_reads_per_gram) {
    index_t *ind = (index_t *)malloc(sizeof(index_t)); /* The q-gram index */
    long long int i, k, j;
    long long int count = 0; /* The total number of q-grams to be inserted to the occurrence lists */
    uint64_t gram_count = 0; /* The number of different q-grams. */
    uint64_t gram;           /* Succinct representation of a q-gram */
    uint64_t gram_reverse;   /* Succinct representation of the reverse of a q-gram */
    uint64_t mask;           /* Bit mask for masking the succinct representations */

    int bits;                /* Number of bits per base, i.e. 2. */

    int lastN;               /* Index to the last indeterminate nucletide encountered. */

    /*uint64_t *gram_array;     Initial q-gram id to q-gram mapping */
    uint64_t *count_array;   /* Initial q-gram id to number of occurrences mapping */

    uint64_t sum;

    uint64_t pre_gram;       /* Succinct presentation of the prefix */
    uint64_t pre_mask;       /* Mask for the prefix */

    *index = ind;

    ind->reads = reads;
    ind->num_reads = num_reads;

    if (q > 0) {
        ind->q = q;
        if (ind->q < 3 || ind->q > 31) {
            return 11;
        }
    } else {
        ind->q = DEFAULT_Q;
    }

    printf("Length of k-mers: %d\n", ind->q);

    if (strlen((char *)prefix) > 0)
        printf("Using prefix %s\n", prefix);

    /* Initialize the succinct presentation of the prefix and a bit
       mask for masking the prefix of a q-gram. */
    pre_gram = 0;
    pre_mask = 0;
    for(i = 0; i < (signed int)strlen((char *)prefix); i++) {
        pre_gram = (pre_gram << 2) | ((prefix[i] >> 1) & 0x03);
        pre_mask = (pre_mask << 2) | 0x3;
        
    }

    /* Count the number of q-grams */
    bits = 2;
    mask = ((uint64_t)1 << (2*ind->q))-1;

    /* gram_array = (uint64_t *)malloc(size*sizeof(uint64_t)); */
    ind->gram2id = new __gnu_cxx::hash_map<uint64_t,uint64_t>;
    count_array = (uint64_t *)malloc(size*sizeof(uint64_t));

    if (ind->gram2id == NULL || count_array == NULL)
        return 1;

    for(i = 0; i < (int64_t)size; i++) {
        /* gram_array[i] = INVALID; */
        count_array[i] = 0;
    }

    printf("Counting k-mers\n");

    gram_count = 0;
    for(k = 0; k < (int64_t)num_reads; k++) {
        /* Print some hint of progress */
        if (k % 100000 == 0)
            printf("%Ld reads, %lu different k-mers\n", k, gram_count);

        gram = 0;
        gram_reverse = 0;
        lastN = -1; /* Initialize so that we will index the first q-gram only after reading q chars */

        for(i = 0; i < (signed int)strlen((char *)reads[k]); i++) {
            if (gram_count >= size) {
                size = 2*size;
                count_array = (uint64_t *)realloc(count_array, 
                                                  size*sizeof(uint64_t));

                if (count_array == NULL)
                    return 1;

                for(j = size/2; j < (int64_t)size; j++) {
                    count_array[j] = 0;
                }
            }

            if (ind->reads[k][i] == 'N' || ind->reads[k][i] == 'n') {
                lastN = i;
            }

            switch(ind->reads[k][i]) {
            case 'A':
            case 'a':
                gram = gram << bits | CODE_A;
                gram_reverse = gram_reverse >> bits | ((uint64_t)CODE_T << ((ind->q-1)*bits));
                break;
            case 'C':
            case 'c':
                gram = gram << bits | CODE_C;
                gram_reverse = gram_reverse >> bits | ((uint64_t)CODE_G << ((ind->q-1)*bits));
                break;
            case 'G':
            case 'g':
                gram = gram << bits | CODE_G;
                gram_reverse = gram_reverse >> bits | ((uint64_t)CODE_C << ((ind->q-1)*bits));
                break;
            case 'T':
            case 't':
                gram = gram << bits | CODE_T;
                gram_reverse = gram_reverse >> bits | ((uint64_t)CODE_A << ((ind->q-1)*bits));
                break;
            default:
                gram = gram << bits;
                gram_reverse = gram_reverse >> bits;
                break;
            }

            gram = gram & mask;
            gram_reverse = gram_reverse & mask;

            if (gram < gram_reverse) {
                if (lastN <= i-ind->q && ((gram & pre_mask) == pre_gram)) {
                    if ((*ind->gram2id)[gram] == 0) {
                        /* A new q-gram */
                        gram_count++;
                        (*ind->gram2id)[gram] = gram_count;
                        /* gram_array[(*ind->gram2id)[gram]-1] = gram; */
                    }

                    count_array[(*ind->gram2id)[gram]-1]++;
                }
            } else if (gram_reverse < gram) {
                if (lastN <= i-ind->q && ((gram_reverse & pre_mask) == pre_gram)) {
                    if ((*ind->gram2id)[gram_reverse] == 0) {
                        /* A new q-gram */
                        gram_count++;
                        (*ind->gram2id)[gram_reverse] = gram_count;
                        /* gram_array[(*ind->gram2id)[gram_reverse]-1] = gram_reverse; */
                    }

                    count_array[(*ind->gram2id)[gram_reverse]-1]++;
                }
            }
        }
    }

    ind->size = gram_count;

    count = 0;
    for(i = 0; i < (int64_t) gram_count; i++) {
        if (count_array[i] >= (uint64_t)min_reads_per_gram &&
            count_array[i] <= (uint64_t)max_reads_per_gram) {
            count += count_array[i];
        }
    }

    printf("%lu different k-mers, %Lu total k-mers\n", gram_count, count);


    /* use arrays of correct size for the real index */
    ind->read_id_array = (ulong *)malloc(count*sizeof(ulong));
    /* ind->grams = (uint64_t *)malloc(gram_count*sizeof(uint64_t)); */
    ind->counts = (uint64_t *)malloc(gram_count*sizeof(uint64_t));
    ind->gram_reads = (ulong **)malloc(gram_count*sizeof(ulong *));

    if (ind->read_id_array == NULL || /* ind->grams == NULL || */
        ind->counts == NULL || ind->gram_reads == NULL) {
        return 1;
    }

    for(i = 0; i < (int64_t)gram_count; i++) {
        /* ind->grams[i] = gram_array[i]; */
        ind->counts[i] = 0;
    }
    /* free(gram_array); */

    /* Initialize the gram_reads array */
    sum = 0;
    for(i = 0; i < (int64_t)gram_count; i++) {
        if (count_array[i] >= (uint64_t) min_reads_per_gram && 
            count_array[i] <= (uint64_t)max_reads_per_gram) {
            ind->gram_reads[i] = &ind->read_id_array[sum];
            sum += count_array[i];
        } else {
            ind->gram_reads[i] = NULL;
        }
    }

    free(count_array);

    printf("Constructing read lists\n");

    for(k = 0; k < (int64_t)num_reads; k++) {
        /* Print a hint of progress */
        if (k % 100000 == 0)
            printf("%Ld reads\n", k);

        gram = 0;
        gram_reverse = 0;
        lastN = -1; /* Initialize so that we will index the first q-gram only after reading q chars */
        for(i = 0; i < (signed int)strlen((char *)reads[k]); i++) {

            if (ind->reads[k][i] == 'N' || ind->reads[k][i] == 'n') {
                lastN = i;
            }

            switch(ind->reads[k][i]) {
            case 'A':
            case 'a':
                gram = gram << bits | CODE_A;
                gram_reverse = gram_reverse >> bits | ((uint64_t)CODE_T << ((ind->q-1)*bits));
                break;
            case 'C':
            case 'c':
                gram = gram << bits | CODE_C;
                gram_reverse = gram_reverse >> bits | ((uint64_t)CODE_G << ((ind->q-1)*bits));
                break;
            case 'G':
            case 'g':
                gram = gram << bits | CODE_G;
                gram_reverse = gram_reverse >> bits | ((uint64_t)CODE_C << ((ind->q-1)*bits));
                break;
            case 'T':
            case 't':
                gram = gram << bits | CODE_T;
                gram_reverse = gram_reverse >> bits | ((uint64_t)CODE_A << ((ind->q-1)*bits));
                break;
            default:
                gram = gram << bits;
                gram_reverse = gram_reverse >> bits;
                break;
            }

            gram = gram & mask;
            gram_reverse = gram_reverse & mask;

            if (gram < gram_reverse) {
                if (lastN <= i-ind->q  && ((gram & pre_mask) == pre_gram)) {
                    if (ind->gram_reads[(*ind->gram2id)[gram]-1] != NULL) {
                        ind->gram_reads[(*ind->gram2id)[gram]-1][ind->counts[(*ind->gram2id)[gram]-1]] = k;
                        ind->counts[(*ind->gram2id)[gram]-1]++;
                    }
                }
            } else if (gram_reverse < gram) {
                if (lastN <= i-ind->q  && ((gram_reverse & pre_mask) == pre_gram)) {
                    if (ind->gram_reads[(*ind->gram2id)[gram_reverse]-1] != NULL) {
                        ind->gram_reads[(*ind->gram2id)[gram_reverse]-1]
                            [ind->counts[(*ind->gram2id)[gram_reverse]-1]] = k;
                        ind->counts[(*ind->gram2id)[gram_reverse]-1]++;
                    }
                }
            }
        }
    }

    return 0;
}


/**
 * Free the space used by the index.
 */
int free_index(void *index) {
    index_t *ind = (index_t *)index;

    free(ind->read_id_array);
    /* free(ind->grams); */
    free(ind->counts);
    free(ind->gram_reads);
    (*ind->gram2id).clear();
    delete ind->gram2id;
    free(ind);

    return 0;
}

/**
 * Return an error description for an error code.
 */
const char *error_index(int e) {
    return error_strings[e];
}

