
CC = icc
#CHECK = -DCHECK
CFLAGS = -O3 $(CHECK) -xHost -ansi-alias -ipo -fp-model precise -openmp

all: exe_2d 
exe_2d: exe_2d5p exe_2d9p exe_2dgol

.PHONY: clean 


exe_2d5p: 2d.np.c
	$(CC) $(CFLAGS) -Dpoint=5 $^ -o $@

exe_2d9p: 2d.np.c
	$(CC) $(CFLAGS) -Dpoint=9 $^ -o $@

exe_2dgol: 2d.np.c
	$(CC) $(CFLAGS) -Dpoint=0 $^ -o $@

clean:
	rm -f   exe_*

