
CC = icc
#CHECK = -DCHECK
CFLAGS = -O3 $(CHECK) -xHost -ansi-alias -ipo -fp-model precise -openmp

all: exe_1d 
exe_1d: exe_1d3p exe_1d5p

.PHONY: clean 

exe_1d3p: 1d.np.c
	$(CC) $(CFLAGS) -Dpoint=3 $^ -o $@ 
exe_1d5p: 1d.np.c
	$(CC) $(CFLAGS) -Dpoint=5 $^ -o $@ 

clean:
	rm -f   exe_* 

