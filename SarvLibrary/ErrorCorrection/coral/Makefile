CFLAGS	= -O3 -Wall -g -fopenmp -Wno-deprecated -m64 # -DDEBUG -DDEBUG_ALIGN
LDFLAGS	= -fopenmp -lm -m64
CC	= g++

CFILES	= multi-align.c q-gram-index-ht.c reverse.c error-correct.c
HFILES  = define.h multi-align.h q-gram-index-ht.h reverse.h
EXTRA = Makefile COPYING

OBJS	= $(CFILES:.c=.o)

PROG	= coral
VERSION = 1.4.1

all:		$(PROG)

$(PROG):	$(OBJS)
		$(CC) $(LDFLAGS) -o $@ $(OBJS)

$(OBJS):	%.o: %.c
		$(CC) $(CFLAGS) -c $<

clean:
		rm $(OBJS)

dist:
		mkdir $(PROG)-$(VERSION)
		cp $(CFILES) $(PROG)-$(VERSION)/
		cp $(HFILES) $(PROG)-$(VERSION)/
		cp $(EXTRA) $(PROG)-$(VERSION)/
		sed 's/VERSION/$(VERSION)/g' < README > $(PROG)-$(VERSION)/README
		tar zcvf $(PROG)-$(VERSION).tar.gz $(PROG)-$(VERSION)

dist-clean:
		rm -r $(PROG)-$(VERSION)/
