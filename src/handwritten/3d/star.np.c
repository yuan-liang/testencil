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
#define kernel(A,t) A[(t+1)%2][x][y][z] = 0.54 * (A[(t)%2][x][y][z]) + \
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
#define  kernel(A,t) A[(t+1)%2][x][y][z] = 0.64 * (A[(t)%2][x][y][z]) + \
										   0.06 * (A[(t)%2][x - 1][y][z] + A[(t)%2][x][y - 1][z] + \
												   A[(t)%2][x][y][z - 1] + A[(t)%2][x + 1][y][z] + \
												   A[(t)%2][x][y + 1][z] + A[(t)%2][x][y][z + 1]);
#define XSLOPE 1
#define YSLOPE 1
#define ZSLOPE 1
#endif

#ifdef CHECK
#define TOLERANCE  0
#endif

#if defined(NX) && defined(NY) && defined(NZ)
double A[2][NX+2*XSLOPE][NY+2*YSLOPE][NZ+2*ZSLOPE];
#endif

int main(int argc, char * argv[])
{

	struct timeval start, end;

	long int t, i, j, k;
	int NX = atoi(argv[1]);
	int NY = atoi(argv[2]);
	int NZ = atoi(argv[3]);
	int T  = atoi(argv[4]);
	int Bx = atoi(argv[5]);
	int tb = atoi(argv[6]);


	if(Bx<(2*XSLOPE+1) || Bx>NX || tb>((Bx-1)/2)/XSLOPE)
	{
		return 0;
	}

#if !defined(NX)
	double (*A)[NX+2*XSLOPE][NY+2*YSLOPE][NZ+2*ZSLOPE] = (double (*)[NX+2*XSLOPE][NY+2*YSLOPE][NZ+2*ZSLOPE])malloc(sizeof(double)*(NX+2*XSLOPE)*(NY+2*YSLOPE)*(NZ+2*ZSLOPE)*2);
	if(NULL == A) return 0;
#endif

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
	int ix=Bx+bx;//x方向跨度


	//奇数时间层,2个空间层,B0和B2,B11和B12对应xleft相反 ybottom相同
	/*   int xleft0[2] = {XSLOPE, XSLOPE-(ix/2)};
		 int ybottom0[2] = {XSLOPE-Bx/2, XSLOPE+(bx/2+1)};
		 int xleft11[2] = {XSLOPE+Bx/2, XSLOPE-(bx/2+1)};
		 int ybottom11[2] = {XSLOPE, XSLOPE+(ix/2)};
		 int xleft12[2] = {xleft11[1], xleft11[0]};
		 int ybottom12[2] = {ybottom11[0], ybottom11[1]};
		 int xleft2[2] = {xleft0[1], xleft0[0]};
		 int ybottom2[2] = {ybottom0[0], ybottom0[1]};
	 */
	int xleft0[2] = {XSLOPE+(bx/2+1), XSLOPE-(Bx/2)};
	int ybottom0[2] = {XSLOPE-(Bx/2), XSLOPE+(bx/2+1)};
	int xleft11[2] = {XSLOPE+(ix/2), XSLOPE};
	int ybottom11[2] = {XSLOPE, XSLOPE+(ix/2)};
	int xleft12[2] = {xleft11[1], xleft11[0]};
	int ybottom12[2] = {ybottom11[0], ybottom11[1]};
	int xleft2[2] = {xleft0[1], xleft0[0]};
	int ybottom2[2] = {ybottom0[0], ybottom0[1]};


	//偶数时间层的B0,B11,B12,B2分别在奇数层B2,B12,B11,B0的位置

	int xnb0[2] = {(NX+XSLOPE-1-xleft0[0])/ix+1, (NX+XSLOPE-1-xleft0[1])/ix+1};//B0在x方向的个数, 空间奇偶层
	int ynb0[2] = {(NY+XSLOPE-1-ybottom0[0])/ix+1, (NY+XSLOPE-1-ybottom0[1])/ix+1};//B0在y方向的个数
	int xnb11[2] = {(NX+XSLOPE-1-xleft11[0])/ix+1, (NX+XSLOPE-1-xleft11[1])/ix+1};//B11在x方向的个数
	int ynb11[2] = {(NY+XSLOPE-1-ybottom11[0])/ix+1, (NY+XSLOPE-1-ybottom11[1])/ix+1};//B11在y方向的个数
	int xnb12[2] = {xnb11[1], xnb11[0]};
	int ynb12[2] = {ynb11[0], ynb11[1]};
	int xnb2[2] = {xnb0[1], xnb0[0]};
	int ynb2[2] = {ynb0[0], ynb0[1]};

	int nb1[2] = {xnb11[0]*ynb11[0] + xnb11[1]*ynb11[1], xnb12[0]*ynb12[0] + xnb12[1]*ynb12[1]};//奇数时间层B11、B12总个数
	int nb02[2] = {xnb0[0]*ynb0[0] + xnb0[1]*ynb0[1], xnb2[0]*ynb2[0] + xnb2[1]*ynb2[1]};//奇、偶时间层B02总个数
	int xnb1[2] = {xnb11[0], xnb12[0]};
	int xnb02[2] = {xnb0[0], xnb2[0]};



	int level = 0;
	int tt,n;
	int x, y, z, xr;
	register int ymin, ymax;
	int xmin,xmax;
	int dy;

	gettimeofday(&start,0);


	for(tt=-tb; tt< T; tt+=tb)
	{
#pragma omp parallel for schedule(dynamic)  private(xmin,xmax,ymin,ymax,t,x,y,z)
		for(n=0; n < nb02[level]; n++)
		{
			for(t= max(tt,0); t <min( tt + 2*tb,  T); t++)
			{
				if(n<xnb0[level]*ynb0[0])
				{
					xmin = xleft0[level] + (Bx-bx)/2 + (n%xnb0[level])*ix      - tb*XSLOPE + myabs(tt+tb,t+1)*XSLOPE;
					xmax = xleft0[level] + (Bx-bx)/2 + (n%xnb0[level])*ix + bx + tb*XSLOPE - myabs(tt+tb,t+1)*XSLOPE;
					ymin =   ybottom0[0] + (Bx-bx)/2 + (n/xnb0[level])*ix      - tb*XSLOPE + myabs(tt+tb,t+1)*XSLOPE;
					ymax =   ybottom0[0] + (Bx-bx)/2 + (n/xnb0[level])*ix + bx + tb*XSLOPE - myabs(tt+tb,t+1)*XSLOPE;
				}
				else
				{
					xmin = xleft0[1-level] + (Bx-bx)/2 + ((n-xnb0[level]*ynb0[0])%xnb0[1-level])*ix      - tb*XSLOPE + myabs(tt+tb,t+1)*XSLOPE;
					xmax = xleft0[1-level] + (Bx-bx)/2 + ((n-xnb0[level]*ynb0[0])%xnb0[1-level])*ix + bx + tb*XSLOPE - myabs(tt+tb,t+1)*XSLOPE;
					ymin =     ybottom0[1] + (Bx-bx)/2 + ((n-xnb0[level]*ynb0[0])/xnb0[1-level])*ix      - tb*XSLOPE + myabs(tt+tb,t+1)*XSLOPE;
					ymax =     ybottom0[1] + (Bx-bx)/2 + ((n-xnb0[level]*ynb0[0])/xnb0[1-level])*ix + bx + tb*XSLOPE - myabs(tt+tb,t+1)*XSLOPE;
				}
				for(x=max(XSLOPE,xmin); x<min(NX+XSLOPE,xmax); x++)
				{
					for(y=max(ymin + myabs((xmin+xmax-1)/2,x), YSLOPE); y<min(ymax - myabs((xmin+xmax-1)/2,x), NY+YSLOPE); y++)
					{
#pragma simd
						for (z = ZSLOPE; z < NZ+ZSLOPE; z++)
						{
							kernel(A,t);
						}
					}
				}
			}
		}

#pragma omp parallel for  schedule(dynamic)  private(xmin,xmax,ymin,ymax,t,x,y,z,dy,xr,i)
		for(n=0; n <nb1[0] + nb1[1]; n++)
		{
			for(t= tt+tb ; t <min( tt + 2*tb + (bx>1),  T+1); t++)
			{
				if(n<nb1[level])  //B11
				{
					dy = -1;
					if(n<xnb11[level]*ynb11[0])  //空间奇数层B11
					{
						xmin = xleft11[level] + (n%xnb11[level]) * ix + XSLOPE;
						ymax = ybottom11[0] + (n/xnb11[level]) * ix + XSLOPE - 1;
					}
					else
					{
						xmin = xleft11[1-level] + ((n-xnb11[level]*ynb11[0])%xnb11[1-level]) * ix + XSLOPE;
						ymax = ybottom11[1] + ((n-xnb11[level]*ynb11[0])/xnb11[1-level]) * ix + XSLOPE - 1;
					}
				}
				else  //B12
				{
					dy = 1;
					if(n<nb1[level]+xnb12[level]*ynb11[0])
					{
						xmin = xleft12[level] + ((n-nb1[level])%xnb12[level]) * ix + XSLOPE;
						ymax = ybottom11[0] + ((n-nb1[level])/xnb12[level]) * ix + ix/2 - XSLOPE + 1;
					}
					else
					{
						xmin = xleft12[1-level] + ((n-nb1[level]-xnb12[level]*ynb11[0])%xnb12[1-level]) * ix + XSLOPE;
						ymax = ybottom11[1] + ((n-nb1[level]-xnb12[level]*ynb11[0])/xnb12[1-level]) * ix + ix/2 - XSLOPE + 1;
					}
				}
				xmax = xmin + ix/2 - bx/2 - (t - (tt+tb)) - 2*XSLOPE + 1;
				ymin = ymax + dy * (XSLOPE - 1 - ix/2 + bx/2 + XSLOPE + t - (tt+tb));

				for(i=0; i<bx/2+1+t-(tt+tb); i++)
				{
					x = xmin +      i +      max(ymin-dy*i-(NY+YSLOPE-1), 0);
					y = ymin - dy * i + dy * max(ymin-dy*i-(NY+YSLOPE-1), 0);
					
					//                    xr = x + min(NX+XSLOPE-x, 
					//                         (dy==-1) ? ymin - ymax - max(ymin-dy*i-(NY+YSLOPE-1), 0)
					//                                  : ymax - ymin - max(ymax-dy*i-(NY+YSLOPE), 0));
					//printf("%d\t%d\t%d\n",xr, max(xmax+i-(NX+XSLOPE),0), x + ((dy==-1) ? ymin - ymax - max(ymin-dy*i-(NY+YSLOPE-1), 0) - max(xmax+i-(NX+XSLOPE),0) : ymax - ymin - max(max(ymax-dy*i-(NY+YSLOPE), xmax+i-(NX+XSLOPE)), 0)));


					xr = x + ((dy==-1) ? ymin - ymax - max(ymin-dy*i-(NY+YSLOPE-1), 0) - max(xmax+i-(NX+XSLOPE),0)
						           : ymax - ymin - max(max(ymax-dy*i-(NY+YSLOPE), xmax+i-(NX+XSLOPE)), 0));



					//            x = max(     XSLOPE, xmin +      i +      max(((dy==-1) ? ymin-dy*i-(NY+YSLOPE-1) : 0),0));//B12的x不用位移
					//          y = min(NY+YSLOPE-1, ymin - dy * i + dy * max(XSLOPE-(xmin+i),0));
					//           ny = min(max(0, (dy==-1) ? min(ymin-ymax-max(XSLOPE-(xmin+i),0), NY+YSLOPE-1-(ymax-dy*i)) : min(ymax-dy*i, NY+YSLOPE)-(ymin-dy*i+max(XSLOPE-(xmin+i), 0))), NX+XSLOPE-x);

					for(; x<xr; x+=XSLOPE, y+=dy*YSLOPE)
					{
						if(t>tt+tb)
						{
#pragma simd
							for (z = ZSLOPE; z < NZ+ZSLOPE; z++)
							{
								kernel(A,t-1);
							}
						}
						if(t<min(T,tt+2*tb))
						{
#pragma simd
							for (z = ZSLOPE; z < NZ+ZSLOPE; z++)
							{
								kernel(A,t);
							}
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
					kernel(B,t);
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
