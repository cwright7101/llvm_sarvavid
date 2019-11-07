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
#ifndef Q_GRAM_INDEX_HT_H
#define Q_GRAM_INDEX_HT_H

#include "define.h"

#include <ext/hash_map>

/* Codes for the bases */
#define CODE_A 0x00
#define CODE_C 0x01
#define CODE_G 0x02
#define CODE_T 0x03

/*The q-gram index*/
typedef struct {
    ulong size;          /* The number of q-grams in the index */
    uchar q;             /* The length of a q-gram */  
    uchar **reads;       /* The reads for which the index is built */
    ulong num_reads;     /* Number of reads */

    ulong *read_id_array;  /* Space for holding the arrays of read ids for the q-grams*/
    /*uint64_t *grams;     Compressed representations of q-gram (id -> q-gram map) */
    uint64_t *counts;    /* The number of occurrences for each q-gram (id -> count) */
    ulong **gram_reads;    /* Pointers to the beginning to read id arrays for each q-gram */
    __gnu_cxx::hash_map<uint64_t,uint64_t> *gram2id;  /* Mapping from q-gram to its id */
} index_t;

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
                int min_reads_per_gram, int max_reads_per_gram);

/**
 * Free the space used by the index.
 */
int free_index(void *index);

/**
 * Return an error description for an error code.
 */
const char *error_index(int e);

#endif
