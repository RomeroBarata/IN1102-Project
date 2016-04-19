CC=gcc
CFLAGS=-Wall -c -O3
EXEDIR=.

all: $(EXEDIR)/mvfcmv

$(EXEDIR)/mvfcmv: mvfcmv.o
	$(CC) mvfcmv.o -o $@ -lm

mvfcmv.o: mvfcmv.c
	$(CC) $(CFLAGS) mvfcmv.c

clean:
	rm -rf *o $(EXEDIR)/mvfcmv
