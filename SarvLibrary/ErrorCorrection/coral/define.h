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
#ifndef DEFINE_H
#define DEFINE_H

#include <stdint.h>

typedef unsigned char uchar;
typedef unsigned long ulong;

#define DEFAULT_Q 21  /* The default length of a q-gram */

#define MAX_NUM_QGRAMS 1048576 /* The maximum number of q-grams. Used
                                  for allocating some initial memory
                                  structures. */

#define PREFIX "" /* Index only q-grams starting with this prefix.
                     Does not work well with the new speedup
                     techniques. */

#define MAX_ERROR_RATE 0.07 /* Maximum allowed proportion of
                               non-agreeing columns in an alignment */

#define THRESHOLD 0.75 /* The minimum proportion of agreeing reads to
                          consider the consensus trustworthy */

#define QUALITY_THRESHOLD 0.75 /* The sum of qualities for agreeing bases minus the 
                                  quality of the base to be corrected must exceed this 
                                  for correction to occur */

#define MAX_ALIGNED_READS 1000 /* The maximum number of reads alignments
                                  will be computed for. */


#define MAX_READ_LENGTH 10000 /* Maximum length of a read */

#define MAX_INSERTIONS 50 /* Maximum number of positions that can be
                             added to a read because of insertions */

#define NUM_THREADS 8 /* Number of threads to use */

/* Scoring scheme for alignments (gaps in the beginning are free) */
#define GAP_PENALTY 3
#define MM_PENALTY 2
#define MATCH_REWARD 2

/* Default parameters for Illumina reads */
#define GAP_PENALTY_ILLUMINA 1000
#define MM_PENALTY_ILLUMINA 1
#define MATCH_REWARD_ILLUMINA 1

/* Default parameters for 454 reads */
#define GAP_PENALTY_454 3
#define MM_PENALTY_454 2
#define MATCH_REWARD_454 2

/* Schemes for quality scores of corrected bases */
#define ESTIMATE 0
#define KEEP 1
#define FLAT 2

#endif
