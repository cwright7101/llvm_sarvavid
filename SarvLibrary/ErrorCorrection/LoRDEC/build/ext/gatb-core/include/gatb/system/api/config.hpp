/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *   Authors: R.Chikhi, G.Rizk, E.Drezen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

#define INT128_FOUND            1

#define STR_LIBRARY_VERSION     "1.0.7"
#define STR_COMPILATION_DATE    "2019-01-23 12:48:14"
#define STR_COMPILATION_FLAGS   " -O3  -DINT128_FOUND  -DWITH_LAMBDA_EXPRESSION -std=c++0x    -DWITH_MPHF  -DDONT_USE_TR1  -Wno-invalid-offsetof"
#define STR_COMPILER            "/usr/bin/gcc  (6.5.0)"

#define STR_OPERATING_SYSTEM    "Linux-4.15.0-43-generic"

#define KSIZE_1  32
#define KSIZE_2  64
#define KSIZE_3  96   
#define KSIZE_4  128

#define PREC_1  ((KSIZE_1+31)/32)
#define PREC_2  ((KSIZE_2+31)/32)
#define PREC_3  ((KSIZE_3+31)/32)
#define PREC_4  ((KSIZE_4+31)/32)

#ifdef GATB_USE_CUSTOM_ALLOCATOR
    #define CUSTOM_MEM_ALLOC  1
#else
    #define CUSTOM_MEM_ALLOC  0
#endif

#define GATB_HDF5_NB_ITEMS_PER_BLOCK (4*1024)
#define GATB_HDF5_CLEANUP_WORKAROUND 4
