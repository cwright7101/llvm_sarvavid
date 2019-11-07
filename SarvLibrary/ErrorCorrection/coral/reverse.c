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
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "reverse.h"

/* Compute the reverse complement of orig and save to reversed */
void reverse(char *orig, char *reversed) {
    int i;
    int m;

    m = strlen(orig);
    for(i = 0; i < (signed int)strlen(orig); i++) {
        switch(orig[i]) {
        case 'A':
            reversed[m-i-1] = 'T';
            break;
        case 'C':
            reversed[m-i-1] = 'G';
            break;
        case 'G':
            reversed[m-i-1] = 'C';
            break;
        case 'T':
            reversed[m-i-1] = 'A';
            break;
        default:
            reversed[m-i-1] = orig[i];
            break;
        }
    }
    reversed[m] = '\0';
}

#ifdef MAIN
int read_line(FILE *f, char *buf, int len) {
    int i;
    int c;
    
    c = fgetc(f);
    i = 0;
    while(c != EOF && c != '\n') {
        if (i < len-1) {
            buf[i] = c;
            i++;
        }
        c = fgetc(f);
    }
    
    buf[i] = '\0';

    if (c == EOF)
        return 0;
    else
        return 1;
}


int main(int argc, char **argv) {
    char buf[1024];
    char buf2[1024];

    FILE *f = fopen(argv[1], "r");

    while(read_line(f, buf, 1024)) {
        printf("%s\n", buf);

        read_line(f, buf, 1024);
        reverse_base(buf, buf2);
        printf("%s\n", buf2);
    }

    fclose(f);

    return 0;
}
#endif
