
CC = icc
#CHECK = -DCHECK
CFLAGS = -O3 $(CHECK) -xHost -ansi-alias -ipo -fp-model precise -qopenmp

all: exe_3d 
exe_3d: exe_3d7p exe_3d27p

.PHONY: clean 


exe_3d7p: 3d.np.c
	$(CC) $(CFLAGS) -Dpoint=7 $^ -o $@

exe_3d27p: 3d.np.c
	$(CC) $(CFLAGS) -Dpoint=27 $^ -o $@

clean:
	rm -f   exe_*

