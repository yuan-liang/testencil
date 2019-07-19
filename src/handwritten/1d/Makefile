
CC = icc
CHECK = -DCHECK
#DEBUG = -Ddebug


ifeq ($(CC), icc)
        CFLAGS = -O3 $(CHECK) $(DEBUG) -xHost -ansi-alias -ipo -fp-model precise -parallel -qopenmp
else
        # for gcc
       CFLAGS   = -O3 $(check) -march=native -mtune=native -ftree-vectorize   -fopenmp
endif


all: exe_1dp exe_1d 
exe_1d: exe_1d3p exe_1d5p
exe_1dp: exe_1d3pp exe_1d5pp

.PHONY: clean 

exe_1d3pp: 1d.p.c
	$(CC) $(CFLAGS) -Dpoint=3 $^ -o $@ 
exe_1d5pp: 1d.p.c
	$(CC) $(CFLAGS) -Dpoint=5 $^ -o $@ 
exe_1d3p: 1d.np.c
	$(CC) $(CFLAGS) -Dpoint=3 $^ -o $@ 
exe_1d5p: 1d.np.c
	$(CC) $(CFLAGS) -Dpoint=5 $^ -o $@ 

clean:
	rm -f   exe_* 

