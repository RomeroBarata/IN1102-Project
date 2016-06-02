CC=gcc
CFLAGS=-Wall -c -O3
EXEDIR=.

all: $(EXEDIR)/mvfcmv

$(EXEDIR)/mvfcmv: util.o mvfcmv.o
	$(CC) util.o mvfcmv.o -o $@ -lm

mvfcmv.o: mvfcmv.c
	$(CC) $(CFLAGS) mvfcmv.c

util.o: util.c
	$(CC) $(CFLAGS) util.c

clean:
	rm -rf *o $(EXEDIR)/mvfcmv
