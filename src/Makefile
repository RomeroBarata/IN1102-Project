CC=gcc
CFLAGS=-Wall -c -O3
EXEDIR=.

all: $(EXEDIR)/mvfcmv

$(EXEDIR)/mvfcmv: stex.o util.o mvfcmv.o
	$(CC) stex.o util.o mvfcmv.o -o $@ -lm

mvfcmv.o: mvfcmv.c
	$(CC) $(CFLAGS) mvfcmv.c

util.o: util.c
	$(CC) $(CFLAGS) util.c

stex.o: stex.c
	$(CC) $(CFLAGS) stex.c

clean:
	rm -rf *o $(EXEDIR)/mvfcmv
