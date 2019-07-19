#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>

#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))
#define myabs(x,y)  ((x) > (y)? ((x)-(y)) : ((y)-(x))) 
#define myceil(x,y)  (int)ceil(((double)x)/((double)y)) // if x and y are integers, myceil(x,y) = (x-1)/y + 1
#define myfloor(x,y)  (int)floor(((double)x)/((double)y)) // if x and y are integers, myceil(x,y) = (x-1)/y + 1


#if !defined(point)
#define point 3
#endif


#if defined(debug)

#if point == 3
#define  kernel(A) if((t == tt+tb) && (x == xmin)) printf("%d\t%d\t%d\t%d\t%d\n",tt,xx,t,xmin,xmax)//A[(t+1)%2][x] = 0.25 * ((A[t%2][x+1] + 2.0 * A[t%2][x]) + A[t%2][x-1])
#define  kernel_b1(A) if((t == tt) && (x == xmin)) printf("%d\t%d\t%d\t%d\t%d\n",tt,xx,t,xmin,xmax)//A[(t+1)%2][x] = 0.25 * ((A[t%2][x+1] + 2.0 * A[t%2][x]) + A[t%2][x-1])
#define  kernel_boundary(A) A[(t+1)%2][x] = 0.25 * ((A[t%2][x==N-1?0:(x+1)] + 2.0 * A[t%2][x]) + A[t%2][x==0?(N-1):(x-1)])
#define XSLOPE  1
#elif point == 5
#define  kernel(A) if(t == tt+tb && x == xmin) printf("%d\t%d\t%d\t%d\t%d\n",tt,xx,t,xmin,xmax)//A[(t+1)%2][x] = 0.25 * ((A[t%2][x+1] + 2.0 * A[t%2][x]) + A[t%2][x-1])
#define  kernel_b1(A) if(t == tt && x == xmin) printf("%d\t%d\t%d\t%d\t%d\n",tt,xx,t,xmin,xmax)//A[(t+1)%2][x] = 0.25 * ((A[t%2][x+1] + 2.0 * A[t%2][x]) + A[t%2][x-1])
#define  kernel_boundary(A)  A[(t+1)%2][x] = 0.125 * (1.4*A[t%2][x<2?(N-2+x):(x-2)] + 1.6*A[t%2][x==0?(N-1):(x-1)] + 2.0 * A[t%2][x] + 1.9*A[t%2][x==N-1?0:(x+1)] + 1.1*A[t%2][(x>N-3)?(x-N+2):(x+2)]);
#define XSLOPE  2
#endif

#else

#if point == 3
#define  kernel(A) A[(t+1)%2][x] = 0.25 * ((A[t%2][x+1] + 2.0 * A[t%2][x]) + A[t%2][x-1])
#define  kernel_boundary(A) A[(t+1)%2][x] = 0.25 * ((A[t%2][x==N-1?0:(x+1)] + 2.0 * A[t%2][x]) + A[t%2][x==0?(N-1):(x-1)])
#define XSLOPE  1
#elif point == 5
#define  kernel(A)  A[(t+1)%2][x] = 0.125 * (1.4*A[t%2][x-2] + 1.6*A[t%2][x-1] + 2.0 * A[t%2][x] + 1.9*A[t%2][x+1] + 1.1*A[t%2][x+2]);
#define  kernel_boundary(A)  A[(t+1)%2][x] = 0.125 * (1.4*A[t%2][x<2?(N-2+x):(x-2)] + 1.6*A[t%2][x==0?(N-1):(x-1)] + 2.0 * A[t%2][x] + 1.9*A[t%2][x==N-1?0:(x+1)] + 1.1*A[t%2][(x>N-3)?(x-N+2):(x+2)]);
#define XSLOPE  2
#endif

#endif

#ifdef CHECK
#define TOLERANCE 0
#endif

int main(int argc, char * argv[]) {

	struct timeval start, end;

	long int  i;
	int N  = atoi(argv[1]);
	int T  = atoi(argv[2]);
	int Bx = atoi(argv[3]);
	int tb = atoi(argv[4]);

	if(Bx<(2*XSLOPE+1) || Bx>N || tb>(((Bx-1)/2)/XSLOPE)){ 
		return 0;
	}

	double (*A)[N] = (double (*)[N+2*XSLOPE])malloc(sizeof(double)*(N+2*XSLOPE)*2);
#ifdef CHECK
	double (*B)[N] = (double (*)[N+2*XSLOPE])malloc(sizeof(double)*(N+2*XSLOPE)*2);
#endif

	srand(100); 
	for (i = 0; i < N+2*XSLOPE; i++) {
		A[0][i] = 1.0 * (rand() % 1024);
		A[1][i] = 0;
#ifdef CHECK
		B[0][i] = A[0][i];
		B[1][i] = 0;
#endif
	}


	int bx = Bx - 2 * tb * XSLOPE;
	if(bx < 2 * XSLOPE) return 0; // bx_first_B1 + bx_last_B1 must be larger than 2 * XSLOPE.

	int ix = Bx + bx;

	int nblock = myfloor(N, ix);
	int nrestpoints = N % ix;

	int bx_first_B1 = (bx + nrestpoints)/2;
	int bx_last_B1  = (bx + nrestpoints) - bx_first_B1;

	int xright[2] = {Bx + bx_first_B1,  Bx + bx_first_B1 + (Bx + bx)/2};


#if defined(debug)
	printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\n",Bx,bx,nblock,bx_first_B1,bx_last_B1,xright[0],xright[1]);
#endif

	int level = 0;
	int x, xx, t, tt;
	register int xmin, xmax;

	gettimeofday(&start, 0);

	for (tt = -tb; tt < T ;  tt += tb ){
#pragma omp parallel for private(xmin,xmax,t,x)
		for(xx = 0; xx <nblock; xx++) {
			if ( (level == 0) || (xx != nblock - 1) ){
				for(t= max(tt, 0) ; t <min( tt + 2*tb,  T); t++){
					xmin = max( 0, xright[level] - Bx + xx*ix + myabs((tt+tb),(t+1))*XSLOPE);
					xmax = min( N, xright[level]      + xx*ix - myabs((tt+tb),(t+1))*XSLOPE);
#pragma ivdep
#pragma vector always
					for(x = xmin; x < xmax; x++){
						kernel(A);
					}
				}
			}
			else{
				for(t=tt ; t <min( tt + 2*tb,  T); t++){
					xmax =     bx_first_B1 + tb * XSLOPE - myabs((tt+tb),(t+1))*XSLOPE;
					xmin = N - bx_last_B1  - tb * XSLOPE + myabs((tt+tb),(t+1))*XSLOPE;
					for(x = 0; x < XSLOPE; x++){
						kernel_boundary(A);
					}
#pragma ivdep
#pragma vector always
					for(; x < xmax; x++){
						kernel(A);
					}
#pragma ivdep
#pragma vector always
					for(x = xmin; x < N-XSLOPE; x++){
						kernel(A);
					}
					for(; x < N; x++){
						kernel_boundary(A);
					}
				}
			}
		}
		level = 1 - level;
	}

	gettimeofday(&end, 0);
	printf("GStencil/s = %f\n",((double)N * T) / (double)(end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec) * 1.0e-6) / 1000000000L);

#ifdef CHECK
	for (t = 0; t < T; t++) {
		for (x = 0; x < N; x++) {
			kernel_boundary(B);
		}
	}
	for (i = 0; i < N; i++) {
		if(myabs(A[T%2][i], B[T%2][i]) > TOLERANCE)
			printf("Naive[%d] = %f, Check = %f: FAILED!\n", i, B[T%2][i], A[T%2][i]);
	}
#endif
}
