#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>

#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))
#define myabs(x,y)  (((x) > (y))? ((x)-(y)) : ((y)-(x)))
#define myceil(x,y)  (int)ceil(((double)x)/((double)y)) // if x and y are integers, myceil(x,y) = (x-1)/y + 1
#define myfloor(x,y)  (int)floor(((double)x)/((double)y)) // if x and y are integers, myceil(x,y) = (x-1)/y + 1

#if !defined(point) 
#define point 5
#endif


int updatedpionts = 0;

#if point == 5

#ifdef CHECK

#define  kernel(A) updatedpionts++;A[(t+1)%2][x][y] = 0.125 * (A[t%2][x+1][y] - 2.0 * A[t%2][x][y] + A[t%2][x-1][y]) + \
									  0.125 * (A[t%2][x][y+1] - 2.0 * A[t%2][x][y] + A[t%2][x][y-1]) + \
A[t%2][x][y];
#define  kernel_boundary_left(A) updatedpionts++;A[(t+1)%2][x][y] = 0.125 * (A[t%2][x==NX-1?0:(x+1)][y] - 2.0 * A[t%2][x][y] + A[t%2][x==0?(NX-1):(x-1)][y]) + \
													0.125 * (A[t%2][x][y==NY-1?0:(y+1)] - 2.0 * A[t%2][x][y] + A[t%2][x][y==0?(NY-1):(y-1)]) + \
A[t%2][x][y];
#define  kernel_boundary_right(A) updatedpionts++;A[(t+1)%2][x][y] = 0.125 * (A[t%2][x==NX-1?0:(x+1)][y] - 2.0 * A[t%2][x][y] + A[t%2][x==0?(NX-1):(x-1)][y]) + \
													 0.125 * (A[t%2][x][y==NY-1?0:(y+1)] - 2.0 * A[t%2][x][y] + A[t%2][x][y==0?(NY-1):(y-1)]) + \
A[t%2][x][y];
#define  kernel_boundary_top(A) updatedpionts++;A[(t+1)%2][x][y] = 0.125 * (A[t%2][x==NX-1?0:(x+1)][y] - 2.0 * A[t%2][x][y] + A[t%2][x==0?(NX-1):(x-1)][y]) + \
												   0.125 * (A[t%2][x][y==NY-1?0:(y+1)] - 2.0 * A[t%2][x][y] + A[t%2][x][y==0?(NY-1):(y-1)]) + \
A[t%2][x][y];
#define  kernel_boundary_bottom(A) updatedpionts++;A[(t+1)%2][x][y] = 0.125 * (A[t%2][x==NX-1?0:(x+1)][y] - 2.0 * A[t%2][x][y] + A[t%2][x==0?(NX-1):(x-1)][y]) + \
													  0.125 * (A[t%2][x][y==NY-1?0:(y+1)] - 2.0 * A[t%2][x][y] + A[t%2][x][y==0?(NY-1):(y-1)]) + \
A[t%2][x][y];
#define  kernel_boundary_left_bottom(A) updatedpionts++;A[(t+1)%2][x][y] = 0.125 * (A[t%2][x==NX-1?0:(x+1)][y] - 2.0 * A[t%2][x][y] + A[t%2][x==0?(NX-1):(x-1)][y]) + \
														   0.125 * (A[t%2][x][y==NY-1?0:(y+1)] - 2.0 * A[t%2][x][y] + A[t%2][x][y==0?(NY-1):(y-1)]) + \
A[t%2][x][y];
#define  kernel_boundary_left_top(A) updatedpionts++;A[(t+1)%2][x][y] = 0.125 * (A[t%2][x==NX-1?0:(x+1)][y] - 2.0 * A[t%2][x][y] + A[t%2][x==0?(NX-1):(x-1)][y]) + \
														0.125 * (A[t%2][x][y==NY-1?0:(y+1)] - 2.0 * A[t%2][x][y] + A[t%2][x][y==0?(NY-1):(y-1)]) + \
A[t%2][x][y];
#define  kernel_boundary_right_bottom(A) updatedpionts++;A[(t+1)%2][x][y] = 0.125 * (A[t%2][x==NX-1?0:(x+1)][y] - 2.0 * A[t%2][x][y] + A[t%2][x==0?(NX-1):(x-1)][y]) + \
															0.125 * (A[t%2][x][y==NY-1?0:(y+1)] - 2.0 * A[t%2][x][y] + A[t%2][x][y==0?(NY-1):(y-1)]) + \
A[t%2][x][y];
#define  kernel_boundary_right_top(A) updatedpionts++;A[(t+1)%2][x][y] = 0.125 * (A[t%2][x==NX-1?0:(x+1)][y] - 2.0 * A[t%2][x][y] + A[t%2][x==0?(NX-1):(x-1)][y]) + \
														 0.125 * (A[t%2][x][y==NY-1?0:(y+1)] - 2.0 * A[t%2][x][y] + A[t%2][x][y==0?(NY-1):(y-1)]) + \
A[t%2][x][y];
#define  kernel_boundary(A) updatedpionts++;A[(t+1)%2][x][y] = 0.125 * (A[t%2][x==NX-1?0:(x+1)][y] - 2.0 * A[t%2][x][y] + A[t%2][x==0?(NX-1):(x-1)][y]) + \
											   0.125 * (A[t%2][x][y==NY-1?0:(y+1)] - 2.0 * A[t%2][x][y] + A[t%2][x][y==0?(NY-1):(y-1)]) + \
A[t%2][x][y];
#else

#define  kernel(A) A[(t+1)%2][x][y] = 0.125 * (A[t%2][x+1][y] - 2.0 * A[t%2][x][y] + A[t%2][x-1][y]) + \
																		  0.125 * (A[t%2][x][y+1] - 2.0 * A[t%2][x][y] + A[t%2][x][y-1]) + \
A[t%2][x][y];
#define  kernel_boundary_left(A) A[(t+1)%2][x][y] = 0.125 * (A[t%2][x==NX-1?0:(x+1)][y] - 2.0 * A[t%2][x][y] + A[t%2][x==0?(NX-1):(x-1)][y]) + \
																										0.125 * (A[t%2][x][y==NY-1?0:(y+1)] - 2.0 * A[t%2][x][y] + A[t%2][x][y==0?(NY-1):(y-1)]) + \
A[t%2][x][y];
#define  kernel_boundary_right(A) A[(t+1)%2][x][y] = 0.125 * (A[t%2][x==NX-1?0:(x+1)][y] - 2.0 * A[t%2][x][y] + A[t%2][x==0?(NX-1):(x-1)][y]) + \
																										 0.125 * (A[t%2][x][y==NY-1?0:(y+1)] - 2.0 * A[t%2][x][y] + A[t%2][x][y==0?(NY-1):(y-1)]) + \
A[t%2][x][y];
#define  kernel_boundary_top(A) A[(t+1)%2][x][y] = 0.125 * (A[t%2][x==NX-1?0:(x+1)][y] - 2.0 * A[t%2][x][y] + A[t%2][x==0?(NX-1):(x-1)][y]) + \
																								   0.125 * (A[t%2][x][y==NY-1?0:(y+1)] - 2.0 * A[t%2][x][y] + A[t%2][x][y==0?(NY-1):(y-1)]) + \
A[t%2][x][y];
#define  kernel_boundary_bottom(A) A[(t+1)%2][x][y] = 0.125 * (A[t%2][x==NX-1?0:(x+1)][y] - 2.0 * A[t%2][x][y] + A[t%2][x==0?(NX-1):(x-1)][y]) + \
																										  0.125 * (A[t%2][x][y==NY-1?0:(y+1)] - 2.0 * A[t%2][x][y] + A[t%2][x][y==0?(NY-1):(y-1)]) + \
A[t%2][x][y];
#define  kernel_boundary_left_bottom(A) A[(t+1)%2][x][y] = 0.125 * (A[t%2][x==NX-1?0:(x+1)][y] - 2.0 * A[t%2][x][y] + A[t%2][x==0?(NX-1):(x-1)][y]) + \
																												   0.125 * (A[t%2][x][y==NY-1?0:(y+1)] - 2.0 * A[t%2][x][y] + A[t%2][x][y==0?(NY-1):(y-1)]) + \
A[t%2][x][y];
#define  kernel_boundary_left_top(A) A[(t+1)%2][x][y] = 0.125 * (A[t%2][x==NX-1?0:(x+1)][y] - 2.0 * A[t%2][x][y] + A[t%2][x==0?(NX-1):(x-1)][y]) + \
																												0.125 * (A[t%2][x][y==NY-1?0:(y+1)] - 2.0 * A[t%2][x][y] + A[t%2][x][y==0?(NY-1):(y-1)]) + \
A[t%2][x][y];
#define  kernel_boundary_right_bottom(A) A[(t+1)%2][x][y] = 0.125 * (A[t%2][x==NX-1?0:(x+1)][y] - 2.0 * A[t%2][x][y] + A[t%2][x==0?(NX-1):(x-1)][y]) + \
																														0.125 * (A[t%2][x][y==NY-1?0:(y+1)] - 2.0 * A[t%2][x][y] + A[t%2][x][y==0?(NY-1):(y-1)]) + \
A[t%2][x][y];
#define  kernel_boundary_right_top(A) A[(t+1)%2][x][y] = 0.125 * (A[t%2][x==NX-1?0:(x+1)][y] - 2.0 * A[t%2][x][y] + A[t%2][x==0?(NX-1):(x-1)][y]) + \
																												 0.125 * (A[t%2][x][y==NY-1?0:(y+1)] - 2.0 * A[t%2][x][y] + A[t%2][x][y==0?(NY-1):(y-1)]) + \
A[t%2][x][y];
#define  kernel_boundary(A) A[(t+1)%2][x][y] = 0.125 * (A[t%2][x==NX-1?0:(x+1)][y] - 2.0 * A[t%2][x][y] + A[t%2][x==0?(NX-1):(x-1)][y]) + \
																						   0.125 * (A[t%2][x][y==NY-1?0:(y+1)] - 2.0 * A[t%2][x][y] + A[t%2][x][y==0?(NY-1):(y-1)]) + \
A[t%2][x][y];


#endif
#define XSLOPE 1
#define YSLOPE 1
#define DATA_TYPE double
#elif point == 9
#define  kernel(A) A[(t+1)%2][x][y] =  0.96 * A[t%2][x][y] + \
									   0.0051 * (A[t%2][x+1][y] +  A[t%2][x-1][y] + A[t%2][x][y+1]+A[t%2][x][y-1]) + \
0.0049 * (A[t%2][x+1][y-1] + A[t%2][x-1][y+1] + A[t%2][x-1][y-1] + A[t%2][x+1][y+1]); 
#define XSLOPE 1
#define YSLOPE 1
#define DATA_TYPE double
#elif point == 0
#define kernel(A)  A[(t+1)%2][x][y] = b2s23(A[t%2][x][y], A[t%2][x-1][y+1] + A[t%2][x-1][y] + \
		A[t%2][x-1][y-1] + A[t%2][x][y+1] + \
		A[t%2][x][y-1]   + A[t%2][x+1][y+1] + \
		A[t%2][x+1][y]   + A[t%2][x+1][y-1]);
#define XSLOPE 1
#define YSLOPE 1
#define DATA_TYPE int
int b2s23(int cell, int neighbors) {
	if((cell == 1 && ((neighbors < 2) || (neighbors > 3)))) {
		return 0;
	}
	if((cell == 1 && (neighbors == 2 || neighbors == 3))) {
		return 1;
	}
	if((cell == 0 && neighbors == 3)) {
		return 1;
	}
	return cell;
}
#endif


#ifdef CHECK
#define TOLERANCE  0
#endif


int main(int argc, char * argv[]) {

	struct timeval start, end;

	long int t, i, j;
	int NX = atoi(argv[1]);
	int NY = atoi(argv[2]);
	int T  = atoi(argv[3]);
	int Bx = atoi(argv[4]);
	int By = atoi(argv[5]);
	int tb = atoi(argv[6]);


	if(Bx<(2*XSLOPE+1) || By<(2*YSLOPE+1) || Bx>NX || By>NY || tb>min(((Bx-1)/2)/XSLOPE,((By-1)/2)/YSLOPE)){
		return 0;
	}

	DATA_TYPE (*A)[NX][NY] = (DATA_TYPE (*)[NX][NY])malloc(sizeof(DATA_TYPE)*(NX)*(NY)*2);
	if(NULL == A) return 0;

#ifdef CHECK
	DATA_TYPE (*B)[NX][NY] = (DATA_TYPE (*)[NX][NY])malloc(sizeof(DATA_TYPE)*(NX)*(NY)*2);
	if(NULL == B) return 0;
#endif

	srand(100);

	for (i = 0; i < NX; i++) {
		for (j = 0; j < NY; j++) {
			A[0][i][j] = 0;//(DATA_TYPE) (1.0 * (rand() % 1024));
#if  point == 0
			A[0][i][j] = ((int)A[0][i][j])%2;
#endif
			A[1][i][j] = 0;
#ifdef CHECK
			B[0][i][j] = A[0][i][j];
			B[1][i][j] = 0;
#endif
		}
	}

	int bx = Bx-2*(tb*XSLOPE);
	int by = By-2*(tb*YSLOPE);

	if(bx < 2 * XSLOPE) return 0;
	if(by < 2 * YSLOPE) return 0;

	int ix = Bx+bx;
	int iy = By+by; // ix and iy are even.

	int xnb0  = NX / ix;
	int ynb0  = NY / iy;

	int xnrestpoints = NX % ix;
	int ynrestpoints = NY % iy;


	int bx_first_B11 = (bx + xnrestpoints)/2;
	int bx_last_B11  = (bx + xnrestpoints) - bx_first_B11;

	int by_first_B12 = (by + ynrestpoints)/2;
	int by_last_B12  = (by + ynrestpoints) - by_first_B12;

	int xnb11 = xnb0;
	int ynb11 = ynb0;
	int xnb12 = xnb0;
	int ynb12 = ynb0;  // for periodic boundary stencils xnb11 and ynb12 must be larger than 2.
	int xnb2  = xnb0;
	int ynb2  = ynb0;

	int nb02[2]  = {xnb0 * ynb0, xnb2 * ynb2};
	int nb1[2]   = {(xnb11-1) * ynb11, xnb12 * (ynb12-1)};
	int xnb02[2] = {xnb0, xnb2-1};
	int xnb1[2]  = {xnb11-1, xnb12};
	/*
	   int xleft02[2]   = {XSLOPE + bx, XSLOPE - (Bx-bx)/2};
	   int ybottom02[2] = {YSLOPE + by, YSLOPE - (By-by)/2};
	   int xleft11[2]   = {XSLOPE,		 XSLOPE + ix/2};
	   int ybottom11[2] = {YSLOPE + by, YSLOPE - (By-by)/2};
	   int xleft12[2]   = {XSLOPE + bx, XSLOPE - (Bx-bx)/2};
	   int ybottom12[2] = {YSLOPE,		 YSLOPE + iy/2};
	 */

	int xleft02[2]   = {bx_first_B11,      bx_first_B11 + ix/2};
	int ybottom02[2] = {by_first_B12,      by_first_B12 + iy/2};
	int xleft11[2]   = {bx_first_B11 + Bx, bx_first_B11 + (Bx-bx)/2};
	int ybottom11[2] = {by_first_B12,      by_first_B12 + (By+by)/2};
	int xleft12[2]   = {bx_first_B11,      bx_first_B11 + (Bx+bx)/2};
	int ybottom12[2] = {by_first_B12 + By, by_first_B12 + (By-by)/2};

	int level = 0;
	int tt,n;
	int x, y;
	register int ymin, ymax;
	int xmin,xmax;

	gettimeofday(&start,0);

	for(tt = -tb; tt < T; tt += tb){
#pragma omp parallel for schedule(dynamic)  private(xmin,xmax,ymin,ymax,t,x,y)
		for(n = 0; n < nb02[level]; n++){
			if(level == 0 || n < nb02[level] - (xnb2 + ynb2 -1)){ // there is no boundary block of B0
				for(t = max(tt,0); t < min(tt + 2*tb, T); t++) { 
					xmin =   xleft02[level] + (n%xnb02[level]) * ix      + myabs(t+1,tt+tb) * XSLOPE;
					xmax =   xleft02[level] + (n%xnb02[level]) * ix + Bx - myabs(t+1,tt+tb) * XSLOPE;
					ymin = ybottom02[level] + (n/xnb02[level]) * iy      + myabs(t+1,tt+tb) * YSLOPE;
					ymax = ybottom02[level] + (n/xnb02[level]) * iy + By - myabs(t+1,tt+tb) * YSLOPE;
					for(x = xmin; x < xmax; x++) {
#pragma ivdep
#pragma vector always
						for(y = ymin; y < ymax; y++) {
							kernel(A);
						}
					}
				}
			}
			else if(n < nb02[level] - xnb2){  // processing ynb2 - 1 left  boundary blocks of B2
				for(t = tt; t < min(tt + 2*tb, T); t++) { 
					xmin = NX - bx_last_B11  - (Bx-bx)/2 + myabs(t+1,tt+tb) * XSLOPE;
					xmax =      bx_first_B11 + (Bx-bx)/2 - myabs(t+1,tt+tb) * XSLOPE;
					ymin = by_first_B12 + iy/2 +      (nb02[level] - xnb2 - n - 1) * iy + myabs(t+1,tt+tb) * YSLOPE;
					ymax = by_first_B12 + iy/2 + By + (nb02[level] - xnb2 - n - 1) * iy - myabs(t+1,tt+tb) * YSLOPE;

					for(x = 0; x < XSLOPE; x++){
						for(y = ymin; y < ymax; y++) {
							kernel_boundary_left(A);
						}
					}
					for(; x < xmax; x++) {
#pragma ivdep
#pragma vector always
						for(y = ymin; y < ymax; y++) {
							kernel(A);
						}
					}
					for(x = xmin; x < NX-XSLOPE; x++) {
#pragma ivdep
#pragma vector always
						for(y = ymin; y < ymax; y++) {
							kernel(A);
						}
					}
					for(; x < NX; x++){
						for(y = ymin; y < ymax; y++) {
							kernel_boundary_right(A); 
						}
					}
				}
			}
			else if(n < nb02[level] - 1){  // processing xnb2 - 1 bottom boundary blocks of B2
				for(t = tt; t < min(tt + 2*tb, T); t++) { 
					ymax = by_first_B12 + (By-by)/2 - myabs(t+1,tt+tb) * YSLOPE;
					ymin = NY - by_last_B12 - (By-by)/2 + myabs(t+1,tt+tb) * YSLOPE;
					xmin = bx_first_B11 + ix/2 +      (nb02[level] - 1 - n - 1) * ix + myabs(t+1,tt+tb) * XSLOPE; 
					xmax = bx_first_B11 + ix/2 + Bx + (nb02[level] - 1 - n - 1) * ix - myabs(t+1,tt+tb) * XSLOPE; 

					for(x = xmin; x < xmax; x++) {
						for(y = 0; y < YSLOPE; y++) {
							kernel_boundary_bottom(A);
						}
					}
					for(x = xmin; x < xmax; x++) {
#pragma ivdep
#pragma vector always
						for(y = YSLOPE; y < ymax; y++) {
							kernel(A);
						}
					}
					for(x = xmin; x < xmax; x++) {
#pragma ivdep
#pragma vector always
						for(y = ymin; y < NY - YSLOPE; y++) {
							kernel(A);
						}
					}
					for(x = xmin; x < xmax; x++) {
						for(y = NY - YSLOPE; y < NY; y++) {
							kernel_boundary_top(A); 
						}
					}
				}
			}
			else{ // processing corner block of B2
				for(t = tt; t < min(tt + 2*tb, T); t++) { 
					xmin = NX - bx_last_B11  - (Bx-bx)/2 + myabs(t+1,tt+tb) * XSLOPE;
					xmax =      bx_first_B11 + (Bx-bx)/2 - myabs(t+1,tt+tb) * XSLOPE;
					ymax =      by_first_B12 + (By-by)/2 - myabs(t+1,tt+tb) * YSLOPE;
					ymin = NY - by_last_B12  - (By-by)/2 + myabs(t+1,tt+tb) * YSLOPE;


					for(x = 0; x < XSLOPE; x++){
						for(y = 0; y < YSLOPE; y++) {
							kernel_boundary_left_bottom(A);
						}
					}
					for(x = 0; x < XSLOPE; x++){
						for(y = YSLOPE; y < ymax; y++) {
							kernel_boundary_left(A);
						}
					}
					for(x = XSLOPE; x < xmax; x++){
						for(y = 0; y < YSLOPE; y++) {
							kernel_boundary_bottom(A);
						}
					}
					for(x = XSLOPE; x < xmax; x++){
#pragma ivdep
#pragma vector always
						for(y = YSLOPE; y < ymax; y++) {
							kernel(A);
						}
					}


					for(x = 0; x < XSLOPE; x++){
						for(y = NY - YSLOPE; y < NY; y++) {
							kernel_boundary_left_top(A);
						}
					}
					for(x = 0; x < XSLOPE; x++){
						for(y = ymin; y < NY - YSLOPE; y++) {
							kernel_boundary_left(A);
						}
					}
					for(x = XSLOPE; x < xmax; x++){
						for(y = NY - YSLOPE; y < NY; y++) {
							kernel_boundary_top(A);
						}
					}
					for(x = XSLOPE; x < xmax; x++){
#pragma ivdep
#pragma vector always
						for(y = ymin; y < NY - YSLOPE; y++) {
							kernel(A);
						}
					}


					for(x = NX - XSLOPE; x < NX; x++){
						for(y = 0; y < YSLOPE; y++) {
							kernel_boundary_right_bottom(A);
						}
					}
					for(x = NX - XSLOPE; x < NX; x++){
						for(y = YSLOPE; y < ymax; y++) {
							kernel_boundary_right(A);
						}
					}
					for(x = xmin; x < NX - XSLOPE; x++){
						for(y = 0; y < YSLOPE; y++) {
							kernel_boundary_bottom(A);
						}
					}
					for(x = xmin; x < NX - XSLOPE; x++){
#pragma ivdep
#pragma vector always
						for(y = YSLOPE; y < ymax; y++) {
							kernel(A);
						}
					}


					for(x = NX - XSLOPE; x < NX; x++){
						for(y = NY - YSLOPE; y < NY; y++) {
							kernel_boundary_right_top(A);
						}
					}
					for(x = NX - XSLOPE; x < NX; x++){
						for(y = ymin; y < NY - YSLOPE; y++) {
							kernel_boundary_right(A);
						}
					}
					for(x = xmin; x < NX - XSLOPE; x++){
						for(y = NY - YSLOPE; y < NY; y++) {
							kernel_boundary_top(A);
						}
					}
					for(x = xmin; x < NX - XSLOPE; x++){
#pragma ivdep
#pragma vector always
						for(y = ymin; y < NY - YSLOPE; y++) {
							kernel(A);
						}
					}
				}
			}
		}

#pragma omp parallel for schedule(dynamic)  private(xmin,xmax,ymin,ymax,t,x,y)
		for(n = 0; n < xnb11 * ynb11 + xnb12 * ynb12; n++) {
			if(n < (xnb11-1) * ynb11 + xnb12 * (ynb12 - 1)){
				for(t = tt + tb; t < min(tt + 2*tb, T); t++) {
					if(n < nb1[level]) {
						xmin = xleft11[level]   + (n%xnb1[level]) * ix       - (t+1-tt-tb) * XSLOPE;
						xmax = xleft11[level]   + (n%xnb1[level]) * ix  + bx + (t+1-tt-tb) * XSLOPE;
						ymin = ybottom11[level] + (n/xnb1[level]) * iy       + (t+1-tt-tb) * YSLOPE;
						ymax = ybottom11[level] + (n/xnb1[level]) * iy  + By - (t+1-tt-tb) * YSLOPE;
					}
					else {
						xmin = xleft12[level]   + ((n-nb1[level])%xnb1[1-level]) * ix      + (t+1-tt-tb) * XSLOPE;
						xmax = xleft12[level]   + ((n-nb1[level])%xnb1[1-level]) * ix + Bx - (t+1-tt-tb) * XSLOPE;
						ymin = ybottom12[level] + ((n-nb1[level])/xnb1[1-level]) * iy      - (t+1-tt-tb) * YSLOPE;
						ymax = ybottom12[level] + ((n-nb1[level])/xnb1[1-level]) * iy + by + (t+1-tt-tb) * YSLOPE;
					}
					for(x = xmin; x < xmax; x++) {
#pragma ivdep
#pragma vector always
						for(y = ymin; y < ymax; y++) {
							kernel(A);
						}
					}
				}
			}
			else if(n < xnb11 * ynb11 + xnb12 * (ynb12 - 1)){ // processing left column B1 blocks
				for(t = tt + tb; t < min(tt + 2*tb, T); t++) {
					if(level == 0) {
						xmin = NX - bx_last_B11 - (t+1-tt-tb) * XSLOPE;
						xmax =     bx_first_B11 + (t+1-tt-tb) * XSLOPE;
						ymin = by_first_B12 + (n - (xnb11-1) * ynb11 - xnb12 * (ynb12 - 1)) * iy       + (t+1-tt-tb) * YSLOPE;
						ymax = by_first_B12 + (n - (xnb11-1) * ynb11 - xnb12 * (ynb12 - 1)) * iy  + By - (t+1-tt-tb) * YSLOPE;
					}
					else {
						xmin = NX - bx_last_B11 - (Bx - bx)/2 + (t+1-tt-tb) * XSLOPE;
						xmax =     bx_first_B11 + (Bx - bx)/2 - (t+1-tt-tb) * XSLOPE;
						ymin = by_first_B12 + (By - by)/2 + (n - (xnb11-1) * ynb11 - xnb12 * (ynb12 - 1)) * iy      - (t+1-tt-tb) * YSLOPE;
						ymax = by_first_B12 + (By - by)/2 + (n - (xnb11-1) * ynb11 - xnb12 * (ynb12 - 1)) * iy + by + (t+1-tt-tb) * YSLOPE;
					}

					for(x = 0; x < XSLOPE; x++){
						for(y = ymin; y < ymax; y++) {
							kernel_boundary_left(A);
						}
					}
					for(x = XSLOPE; x < xmax; x++) {
#pragma ivdep
#pragma vector always
						for(y = ymin; y < ymax; y++) {
							kernel(A);
						}
					}
					for(x = xmin; x < NX-XSLOPE; x++) {
#pragma ivdep
#pragma vector always
						for(y = ymin; y < ymax; y++) {
							kernel(A);
						}
					}
					for(x = NX - XSLOPE; x < NX; x++){
						for(y = ymin; y < ymax; y++) {
							kernel_boundary_right(A); 
						}
					}
				}
			}
			else{ // processing bottom row B1 blocks
				for(t = tt + tb; t < min(tt + 2*tb, T); t++) {
					if(level == 0) {
						xmin = bx_first_B11 + (n - xnb11 * ynb11 - xnb12 * (ynb12 - 1)) * ix      + (t+1-tt-tb) * XSLOPE;
						xmax = bx_first_B11 + (n - xnb11 * ynb11 - xnb12 * (ynb12 - 1)) * ix + Bx - (t+1-tt-tb) * XSLOPE;
						ymin = NY - by_last_B12 - (t+1-tt-tb) * YSLOPE;
						ymax =     by_first_B12 + (t+1-tt-tb) * YSLOPE;
					}
					else {
						xmin = bx_first_B11 + (Bx - bx)/2 + (n - xnb11 * ynb11 - xnb12 * (ynb12 - 1)) * ix      - (t+1-tt-tb) * XSLOPE;
						xmax = bx_first_B11 + (Bx - bx)/2 + (n - xnb11 * ynb11 - xnb12 * (ynb12 - 1)) * ix + bx + (t+1-tt-tb) * XSLOPE;
						ymin = NY - by_last_B12 + (t+1-tt-tb) * YSLOPE - (By - by)/2;
						ymax =     by_first_B12 - (t+1-tt-tb) * YSLOPE + (By - by)/2;
					}

					for(x = xmin; x < xmax; x++) {
						for(y = 0; y < YSLOPE; y++) {
							kernel_boundary_bottom(A);
						}
					}
					for(x = xmin; x < xmax; x++) {
#pragma ivdep
#pragma vector always
						for(y = YSLOPE; y < ymax; y++) {
							kernel(A);
						}
					}
					for(x = xmin; x < xmax; x++) {
#pragma ivdep
#pragma vector always
						for(y = ymin; y < NY - YSLOPE; y++) {
							kernel(A);
						}
					}
					for(x = xmin; x < xmax; x++) {
						for(y = NY-YSLOPE; y < NY; y++) {
							kernel_boundary_top(A); 
						}
					}
				}
			}
		}
		level = 1 - level;
	}

	gettimeofday(&end,0);

#ifdef CHECK	
	printf("%d\t%d\n", updatedpionts, NX * NY * T);
#endif

	printf("GStencil/s = %f\n", ((double)NX * NY * T) / (double)(end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec) * 1.0e-6) / 1000000000L);

#ifdef CHECK
	for (t = 0; t < T; t++) {
		for (x = 0; x < NX; x++) {
			for (y = 0; y < NY; y++) {
				kernel_boundary(B);
			}
		}
	}
	for (i = 0; i < NX; i++) {
		for (j = 0; j < NY; j++) {
			if(myabs(A[T%2][i][j],B[T%2][i][j]) > TOLERANCE)
				printf("Naive[%d][%d] = %f, Check = %f: FAILED!\n", i, j, B[T%2][i][j], A[T%2][i][j]);
		}
	}
#endif
}
