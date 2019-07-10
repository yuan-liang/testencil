#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>

#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))
#define myabs(x,y)  (((x) > (y))? ((x)-(y)) : ((y)-(x)))

#if !defined(point) 
#define point 7
#endif

#if point == 27
#define kernel(A) A[(t+1)%2][x][y][z] = 0.54 * (A[(t)%2][x][y][z]) + \
										0.03 * (A[(t)%2][x][y][z-1] + A[(t)%2][x][y-1][z] + \
												A[(t)%2][x][y+1][z] + A[(t)%2][x][y][z+1] + \
												A[(t)%2][x-1][y][z] + A[(t)%2][x+1][y][z])+ \
										0.01 * (A[(t)%2][x-1][y][z-1] + A[(t)%2][x-1][y-1][z] + \
												A[(t)%2][x-1][y+1][z] + A[(t)%2][x-1][y][z+1] + \
												A[(t)%2][x][y-1][z-1] + A[(t)%2][x][y+1][z-1] + \
												A[(t)%2][x][y-1][z+1] + A[(t)%2][x][y+1][z+1] + \
												A[(t)%2][x+1][y][z-1] + A[(t)%2][x+1][y-1][z] + \
												A[(t)%2][x+1][y+1][z] + A[(t)%2][x+1][y][z+1])+ \
										0.02 * (A[(t)%2][x-1][y-1][z-1] + A[(t)%2][x-1][y+1][z-1] + \
												A[(t)%2][x-1][y-1][z+1] + A[(t)%2][x-1][y+1][z+1] + \
												A[(t)%2][x+1][y-1][z-1] + A[(t)%2][x+1][y+1][z-1] + \
												A[(t)%2][x+1][y-1][z+1] + A[(t)%2][x+1][y+1][z+1]);
#define XSLOPE 1
#define YSLOPE 1
#define ZSLOPE 1
#elif point == 7
#define  kernel(A) A[(t+1)%2][x][y][z] = 0.64 * (A[t%2][x][y][z]) + \
										 0.06 * (A[t%2][x - 1][y][z] + A[t%2][x][y - 1][z] + \
												 A[t%2][x][y][z - 1] + A[t%2][x + 1][y][z] + \
												 A[t%2][x][y + 1][z] + A[t%2][x][y][z + 1]);
#define XSLOPE 1
#define YSLOPE 1
#define ZSLOPE 1
#endif

#ifdef CHECK
#define TOLERANCE  0
#endif

int main(int argc, char * argv[]) {

	struct timeval start, end;

	long int t, i, j, k;
	int NX = atoi(argv[1]);
	int NY = atoi(argv[2]);
	int NZ = atoi(argv[3]);
	int T  = atoi(argv[4]);
	int Bx = atoi(argv[5]);
	int By = atoi(argv[6]);
	int tb = atoi(argv[7]);
	
	if(Bx<(2*XSLOPE+1) || By<(2*YSLOPE+1) || Bx>NX || By>NY || tb>min(((Bx-1)/2)/XSLOPE,((By-1)/2)/YSLOPE)){
		return 0;
	}
	
	double (*A)[NX+2*XSLOPE][NY+2*YSLOPE][NZ+2*ZSLOPE] = (double (*)[NX+2*XSLOPE][NY+2*YSLOPE][NZ+2*ZSLOPE])malloc(sizeof(double)*(NX+2*XSLOPE)*(NY+2*YSLOPE)*(NZ+2*ZSLOPE)*2);
	if(NULL == A) return 0;

#ifdef CHECK
    double (*B)[NX+2*XSLOPE][NY+2*YSLOPE][NZ+2*ZSLOPE] = (double (*)[NX+2*XSLOPE][NY+2*YSLOPE][NZ+2*ZSLOPE])malloc(sizeof(double)*(NX+2*XSLOPE)*(NY+2*YSLOPE)*(NZ+2*ZSLOPE)*2);
	if(NULL == B) return 0;
#endif

	srand(100);

	for (i = 0; i < NX+2*XSLOPE; i++) {
		for (j = 0; j < NY+2*YSLOPE; j++) {
			for (k = 0; k < NZ+2*ZSLOPE; k++) {
				A[0][i][j][k] = (double) (1.0 * (rand() % 1024));
				A[1][i][j][k] = 0;
#ifdef CHECK
				B[0][i][j][k] = A[0][i][j][k];
				B[1][i][j][k] = 0;
#endif
			}
		}
	}



	int bx = Bx-2*(tb*XSLOPE);
	int by = By-2*(tb*YSLOPE);

	int ix = Bx+bx;
	int iy = By+by;

	int xnb0 = ceild(NX,ix);
	int ynb0 = ceild(NY,iy);
	int xnb11 = ceild(NX-ix/2+1,ix) + 1;
	int ynb12 = ceild(NY-iy/2+1,iy) + 1;
	int xnb12 = xnb0;
	int ynb11 = ynb0;
	int xnb2 = max(xnb11,xnb0);
	int ynb2 = max(ynb12,ynb0);

	int nb1[2] = {xnb12 * ynb12, xnb11 * ynb11};
	int nb02[2] = {xnb2 * ynb2, xnb0 * ynb0};
	int xnb1[2] = {xnb12, xnb11};
	int xnb02[2] = {xnb2, xnb0};

	int xleft02[2] = {XSLOPE-bx, XSLOPE+(Bx-bx)/2};
	int ybottom02[2] = {YSLOPE-by, YSLOPE+(By-by)/2};
	int xleft11[2] = {XSLOPE+(Bx-bx)/2, XSLOPE - bx};
	int ybottom11[2] = {YSLOPE-(By+by)/2, YSLOPE};
	int xleft12[2] = {XSLOPE-(Bx+bx)/2, XSLOPE};
	int ybottom12[2] = {YSLOPE+(By-by)/2, YSLOPE-by};

	int level = 1;
	int tt,n;
	int x, y, z;
	int ymin, ymax;
	int xmin, xmax;

	gettimeofday(&start,0);

	for(tt =- tb; tt < T; tt += tb){
#pragma omp parallel for schedule(dynamic)  private(xmin,xmax,ymin,ymax,t,x,y)
		for(n = 0; n < nb02[level]; n++){

			for(t = max(tt,0); t < min(tt + 2*tb, T); t++) { 

				xmin = max(   XSLOPE,   xleft02[level] + (n%xnb02[level]) * ix      - tb*XSLOPE + myabs(t+1,tt+tb) * XSLOPE);
				xmax = min(NX+XSLOPE,   xleft02[level] + (n%xnb02[level]) * ix + bx + tb*XSLOPE - myabs(t+1,tt+tb) * XSLOPE);
				ymin = max(   YSLOPE, ybottom02[level] + (n/xnb02[level]) * iy      - tb*YSLOPE + myabs(t+1,tt+tb) * YSLOPE);
				ymax = min(NY+YSLOPE, ybottom02[level] + (n/xnb02[level]) * iy + by + tb*YSLOPE - myabs(t+1,tt+tb) * YSLOPE);

				for(x = xmin; x < xmax; x++) {
					for(y = ymin; y < ymax; y++) {
#pragma ivdep
#pragma vector always
						for (z = ZSLOPE; z < NZ+ZSLOPE; z++) {
							kernel(A);
						}
					}
				}
			}
		}


#pragma omp parallel for schedule(dynamic)  private(xmin,xmax,ymin,ymax,t,x,y)
		for(n = 0; n < nb1[0] + nb1[1]; n++) {

			for(t = tt + tb; t < min(tt + 2*tb, T); t++) {
				if(n < nb1[level]) {
					xmin = max(     XSLOPE,   xleft11[level] + (n%xnb1[level]) * ix       - (t+1-tt-tb) * XSLOPE);
					xmax = min(NX + XSLOPE,   xleft11[level] + (n%xnb1[level]) * ix  + bx + (t+1-tt-tb) * XSLOPE);
					ymin = max(     YSLOPE, ybottom11[level] + (n/xnb1[level]) * iy       + (t+1-tt-tb) * YSLOPE);
					ymax = min(NY + YSLOPE, ybottom11[level] + (n/xnb1[level]) * iy  + By - (t+1-tt-tb) * YSLOPE);
				}
				else {

					xmin = max(     XSLOPE,   xleft12[level] + ((n-nb1[level])%xnb1[1-level]) * ix      + (t+1-tt-tb) * XSLOPE);
					xmax = min(NX + XSLOPE,   xleft12[level] + ((n-nb1[level])%xnb1[1-level]) * ix + Bx - (t+1-tt-tb) * XSLOPE);
					ymin = max(     YSLOPE, ybottom12[level] + ((n-nb1[level])/xnb1[1-level]) * iy      - (t+1-tt-tb) * YSLOPE);
					ymax = min(NY + YSLOPE, ybottom12[level] + ((n-nb1[level])/xnb1[1-level]) * iy + by + (t+1-tt-tb) * YSLOPE);

				}
				for(x = xmin; x < xmax; x++) {
					for(y = ymin; y < ymax; y++) {
#pragma ivdep
#pragma vector always
						for (z = ZSLOPE; z < NZ+ZSLOPE; z++) {
							kernel(A);
						}
					}
				}
			}
		}
		level = 1- level;
	}

	gettimeofday(&end,0);

	printf("MStencil/s = %f\n", ((double)NX * NY * NZ * T) / (double)(end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec) * 1.0e-6) / 1000000L);

#ifdef CHECK
	for (t = 0; t < T; t++) {
		for (x = XSLOPE; x < NX+XSLOPE; x++) {
			for (y = YSLOPE; y < NY+YSLOPE; y++) {
				for (z = ZSLOPE; z < NZ+ZSLOPE; z++) {
					kernel(B);
				}
			}
		}
	}
	for (i = XSLOPE; i < NX+XSLOPE; i++) {
		for (j = YSLOPE; j < NY+YSLOPE; j++) {
			for (k = ZSLOPE; k < NZ+ZSLOPE; k++) {
				if(myabs(A[T%2][i][j][k], B[T%2][i][j][k]) > TOLERANCE)
					printf("Naive[%d][%d][%d] = %f, Check() = %f: FAILED!\n", i, j, k, B[T%2][i][j][k], A[T%2][i][j][k]);
			}
		}
	}
#endif
}
