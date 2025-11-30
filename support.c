//#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <stdbool.h>
#include <ctype.h>
#define NR_END 1
#define FREE_ARG char*
#define Pi 3.14159265358979323


void nrerror(char error_text[])/* Numerical Recipes standard error handler */{
	printf("...now exiting to system...\n");
	exit(1);
}


double* vector(long nl, long nh)/* allocate a float vector with subscript range v[nl..nh] */ {
	double* v, * t;
	int i;
	v = (double*)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(double)));
	
	if (!v) nrerror("allocation failure in vector()");
	t = v - nl + NR_END;
	for (i = nl; i <= nh; i++) t[i] = 0.0;
	return v - nl + NR_END;
}

long double* ldvector(long nl, long nh)/* allocate a float vector with subscript range v[nl..nh] */ {
	long double* v, * t;
	int i;
	v = (long double*)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(long double)));
	if (!v) nrerror("allocation failure in vector()");
	t = v - nl + NR_END;
	for (i = nl; i <= nh; i++) t[i] = 0.0;
	return v - nl + NR_END;
}


char* cvector(long nl, long nh)/* allocate a float vector with subscript range v[nl..nh] */ {
	char* v, * t;
	int i;
	v = (char*)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(char)));
	if (!v) nrerror("allocation failure in vector()");
	t = v - nl + NR_END;
	for (i = nl; i <= nh; i++) t[i] = 0.0;
	return v - nl + NR_END;
}


float* fvector(long nl, long nh)/* allocate a float vector with subscript range v[nl..nh] */ {
	float* v, * t;
	int i;
	v = (float*)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(float)));
	if (!v) nrerror("allocation failure in vector()");
	t = v - nl + NR_END;
	for (i = nl; i <= nh; i++) t[i] = 0.0;
	return v - nl + NR_END;
}

int* ivector(int nl, int nh)/* allocate a int vector with subscript range v[nl..nh] */ {
	int* v, * t;
	int i;
	v = (int*)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(long)));
	if (!v) nrerror("allocation failure in ivector()");
	t = v - nl + NR_END;
	for (i = nl; i <= nh; i++) t[i] = 0;
	return v - nl + NR_END;
}

float** fmatrix(long nrl, long nrh, long ncl, long nch)/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */ {
	long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	float** m;
	int j, z;
	/* allocate pointers to rows */
	m = (float**)malloc((size_t)((nrow + NR_END) * sizeof(float*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;
	/* allocate rows and set pointers to them */
	m[nrl] = (float*)malloc((size_t)((nrow * ncol + NR_END) * sizeof(float)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	for (i = nrl + 1;i <= nrh;i++) m[i] = m[i - 1] + ncol;
	/* return pointer to array of pointers to rows */
	for (j = nrl; j <= nrh; j++) for (z = ncl; z <= nch; z++) m[j][z] = 0.0;
	return m;
}
double** matrix(long nrl, long nrh, long ncl, long nch)/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */ {
	long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	double** m;
	int j, z;
	/* allocate pointers to rows */
	m = (double**)malloc((size_t)((nrow + NR_END) * sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;
	/* allocate rows and set pointers to them */
	m[nrl] = (double*)malloc((size_t)((nrow * ncol + NR_END) * sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	for (i = nrl + 1;i <= nrh;i++) m[i] = m[i - 1] + ncol;
	/* return pointer to array of pointers to rows */
	for (j = nrl; j <= nrh; j++) for (z = ncl; z <= nch; z++) m[j][z] = 0.0;
	return m;
}

double** matrix2(long nrl, long nrh, long ncl, long nch)/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */ {
	long long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	double** m;
	int j, z;
	/* allocate pointers to rows */
	m = (double**)malloc((size_t)((nrow + NR_END) * sizeof(double*)));
	if (!m) {
		printf("\t1\t");
		nrerror("allocation failure 1 in matrix()");
	}
	m += NR_END;
	m -= nrl;
	/* allocate rows and set pointers to them */
	printf("\nsize=%llu\n", ((nrow * ncol + NR_END) * sizeof(double)));
	m[nrl] = (double*)malloc((size_t)((((long long)(nrow)) * ncol + NR_END) * sizeof(double)));
	if (!m[nrl]) {
		printf("\t2\t");
		nrerror("allocation failure 1 in matrix()");
	}
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	for (i = nrl + 1;i <= nrh;i++) m[i] = m[i - 1] + ncol;
	/* return pointer to array of pointers to rows */
	for (j = nrl; j <= nrh; j++) for (z = ncl; z <= nch; z++) m[j][z] = 0.0;
	return m;
}

long double** ldmatrix(long nrl, long nrh, long ncl, long nch)/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */ {
	long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	long double** m;
	int j, z;
	/* allocate pointers to rows */
	m = (long double**)malloc((size_t)((nrow + NR_END) * sizeof(long double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;
	/* allocate rows and set pointers to them */
	m[nrl] = (long double*)malloc((size_t)((nrow * ncol + NR_END) * sizeof(long double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	for (i = nrl + 1;i <= nrh;i++) m[i] = m[i - 1] + ncol;
	/* return pointer to array of pointers to rows */
	for (j = nrl; j <= nrh; j++) for (z = ncl; z <= nch; z++) m[j][z] = 0.0;
	return m;
}

char** cmatrix(long nrl, long nrh, long ncl, long nch)/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */ {
	long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	char** m;
	int j, z;
	/* allocate pointers to rows */
	m = (char**)malloc((size_t)((nrow + NR_END) * sizeof(char*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;
	/* allocate rows and set pointers to them */
	m[nrl] = (char*)malloc((size_t)((nrow * ncol + NR_END) * sizeof(char)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	for (i = nrl + 1;i <= nrh;i++) m[i] = m[i - 1] + ncol;
	/* return pointer to array of pointers to rows */
	for (j = nrl; j <= nrh; j++) for (z = ncl; z <= nch; z++) m[j][z] = 0.0;
	return m;
}



int** nmatrix(long nrl, long nrh, long ncl, long nch)/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */ {
	long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	int** m;
	int j, z;
	/* allocate pointers to rows */
	m = (int**)malloc((size_t)((nrow + NR_END) * sizeof(int*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;
	/* allocate rows and set pointers to them */
	m[nrl] = (int*)malloc((size_t)((nrow * ncol + NR_END) * sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	for (i = nrl + 1;i <= nrh;i++) m[i] = m[i - 1] + ncol;
	/* return pointer to array of pointers to rows */
	for (j = nrl; j <= nrh; j++) for (z = ncl; z <= nch; z++) m[j][z] = 0.0;
	return m;
}

void free_vector(double* v, long nl, long nh)/* free a float vector allocated with vector() */ {
	free((FREE_ARG)(v + nl - NR_END));
}

void free_ldvector(long double* v, long nl, long nh)/* free a float vector allocated with vector() */ {
	free((FREE_ARG)(v + nl - NR_END));
}

void free_cvector(char* v, long nl, long nh)/* free a float vector allocated with vector() */ {
	free((FREE_ARG)(v + nl - NR_END));
}


void free_fvector(float* v, long nl, long nh)/* free a float vector allocated with vector() */ {
	free((FREE_ARG)(v + nl - NR_END));
}

void free_ivector(int* v, long nl, long nh)/* free a int vector allocated with vector() */ {
	free((FREE_ARG)(v + nl - NR_END));
}

void free_matrix(double** m, long nrl, long nrh, long ncl, long nch)/* free a float matrix allocated by matrix() */ {
	free((FREE_ARG)(m[nrl] + ncl - NR_END));
	free((FREE_ARG)(m + nrl - NR_END));
}

void free_fmatrix(float** m, long nrl, long nrh, long ncl, long nch)/* free a float matrix allocated by matrix() */ {
	free((FREE_ARG)(m[nrl] + ncl - NR_END));
	free((FREE_ARG)(m + nrl - NR_END));
}

void free_ldmatrix(long double** m, long nrl, long nrh, long ncl, long nch)/* free a float matrix allocated by matrix() */ {
	free((FREE_ARG)(m[nrl] + ncl - NR_END));
	free((FREE_ARG)(m + nrl - NR_END));
}


void free_cmatrix(char** m, long nrl, long nrh, long ncl, long nch)/* free a float matrix allocated by matrix() */ {
	free((FREE_ARG)(m[nrl] + ncl - NR_END));
	free((FREE_ARG)(m + nrl - NR_END));
}
void free_nmatrix(int** m, long nrl, long nrh, long ncl, long nch)/* free a float matrix allocated by matrix() */ {
	free((FREE_ARG)(m[nrl] + ncl - NR_END));
	free((FREE_ARG)(m + nrl - NR_END));
}

//=======================================================================================================================
//================================================== Amoeba =============================================================
//=======================================================================================================================

#define NMAX 50000000
#define aSWAP(a,b) { tempr=(a);(a)=(b);(b)=tempr; }


double amotry(double** p, double y[], double psum[], int ndim, double (*funk)(double[]), int ihi, double fac)
/*Extrapolates by a factor fac through the face of the simplex across from the high point, tries
it, and replaces the high point if the new point is better.*/
{
	int j;
	double fac1, fac2, ytry, * ptry;
	ptry = vector(1, ndim);
	fac1 = (1.0 - fac) / ndim;
	fac2 = fac1 - fac;
	for (j = 1;j <= ndim;j++) ptry[j] = psum[j] * fac1 - p[ihi][j] * fac2;
	ytry = (*funk)(ptry); //Evaluate the function at the trial point.
	if (ytry < y[ihi]) { //If it's better than the highest, then replace the highest.
		y[ihi] = ytry;
		for (j = 1;j <= ndim;j++) {
			psum[j] += ptry[j] - p[ihi][j];
			p[ihi][j] = ptry[j];
		}
	}
	free_vector(ptry, 1, ndim);
	return ytry;
}


void amoeba(double** p, double y[], int ndim, double ftol, double (*funk)(double[]), int* nfunk)
/*Multidimensional minimization of the function funk(x) where x[1..ndim] is a vector in ndim
dimensions, by the downhill simplex method of Nelder and Mead. The matrix p[1..ndim+1]
[1..ndim] is input. Its ndim+1 rows are ndim-dimensional vectors which are the vertices of
the starting simplex. Also input is the vector y[1..ndim+1], whose components must be preinitialized
to the values of funk evaluated at the ndim+1 vertices (rows) of p; and ftol the
fractional convergence tolerance to be achieved in the function value (n.b.!). On output, p and
y will have been reset to ndim+1 new points all within ftol of a minimum function value, and
nfunk gives the number of function evaluations taken.*/
{
	double amotry(double** p, double y[], double psum[], int ndim, double (*funk)(double[]), int ihi, double fac);
	int i, ii, ihi, ilo, inhi, j, mpts = ndim + 1;
	double rtol, sum, swap, ysave, ytry, * psum, tempr;
	psum = vector(1, ndim);
	*nfunk = 0;
	for (j = 1;j <= ndim;j++) {
		for (sum = 0.0, i = 1;i <= mpts;i++) sum += p[i][j];
		psum[j] = sum;
	}

	for (;;) {
		ilo = 1;
		ihi = y[1] > y[2] ? (inhi = 2, 1) : (inhi = 1, 2);
		for (i = 1;i <= mpts;i++) {
			if (y[i] <= y[ilo]) ilo = i;
			if (y[i] > y[ihi]) {
				inhi = ihi;
				ihi = i;
			}
			else if (y[i] > y[inhi] && i != ihi) inhi = i;
		}
		rtol = 2.0 * fabs(y[ihi] - y[ilo]) / (fabs(y[ihi]) + fabs(y[ilo]));
		printf("\nrtol=%e", rtol);

		if (rtol < ftol) {
			aSWAP(y[1], y[ilo]);
			for (i = 1;i <= ndim;i++) aSWAP(p[1][i], p[ilo][i]);
			break;
		}
		if (*nfunk >= NMAX) {
			printf("\tNMAX exceeded");
			break;
		}
		*nfunk += 2;
		ytry = amotry(p, y, psum, ndim, funk, ihi, -1.0);
		if (ytry <= y[ilo])
			ytry = amotry(p, y, psum, ndim, funk, ihi, 2.0);
		else if (ytry >= y[inhi]) {
			ysave = y[ihi];
			ytry = amotry(p, y, psum, ndim, funk, ihi, 0.5);
			if (ytry >= ysave) {
				for (i = 1;i <= mpts;i++) {
					if (i != ilo) {
						for (j = 1;j <= ndim;j++)
							p[i][j] = psum[j] = 0.5 * (p[i][j] + p[ilo][j]);
						y[i] = (*funk)(psum);
					}
				}
				*nfunk += ndim;
				for (j = 1;j <= ndim;j++) {
					for (sum = 0.0, i = 1;i <= mpts;i++) sum += p[i][j];
					psum[j] = sum;
				}
			}
		}
		else --(*nfunk);
	}
	free_vector(psum, 1, ndim);
}



//==================================================================================================================
//============================================ Poly Clip ===========================================================
//==================================================================================================================


#define MAXITE 30 //Maximum allowed number of iterations.

float rtsec(double (*func)(double), double x1, double x2, double xacc, int report)
/*Using the secant method, nd the root of a function func thought to lie between x1 and x2.
The root, returned as rtsec, is rened until its accuracy is xacc.*/
{
	void nrerror(char error_text[]);
	int j;
	double fl, f, dx, swap, xl, rts;
	fl = (*func)(x1);
	f = (*func)(x2);
	if (fabs(fl) < fabs(f)) {
		//Pick the bound with the smaller function value as
		rts = x1; //the most recent guess.
		xl = x2;
		swap = fl;
		fl = f;
		f = swap;
	}
	else {
		xl = x1;
		rts = x2;
	}
	for (j = 1;j <= MAXITE;j++) {
		//Secant loop.
		dx = (xl - rts) * f / (f - fl); //Increment with respect to latest value.
		xl = rts;
		fl = f;
		rts += dx;
		f = (*func)(rts);
		if (fabs(dx) < xacc || f == 0.0) return rts; //Convergence.
		if(report) printf("\n%e",fabs(dx));
	}
	nrerror("Maximum number of iterations exceeded in rtsec");
	return 0.0; //Never get here.
}


#define IM1 2147483563
#define IM2 2147483399
#define AM1 (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV1 (1+IMM1/NTAB)
#define EPSS 1.2e-7
#define RNMX (1.0-EPSS)

float ran2(long* idum)/*Long period (> 2  1018) random number generator of L'Ecuyer with Bays-Durham shue
and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
the endpoint values). Call with idum a negative integer to initialize; thereafter, do not alter
idum between successive deviates in a sequence. RNMX should approximate the largest floating
value that is less than 1.*/ {
	int j;
	long k;
	static long idum2 = 123456789;
	static long iy = 0;
	static long iv[NTAB];
	float temp;
	if (*idum <= 0) { //Initialize.
		if (-(*idum) < 1) *idum = 1; //Be sure to prevent idum = 0.
		else *idum = -(*idum);
		idum2 = (*idum);
		for (j = NTAB + 7;j >= 0;j--) { //Load the shue table (after 8 warm-ups).
			k = (*idum) / IQ1;
			*idum = IA1 * (*idum - k * IQ1) - k * IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy = iv[0];
	}
	k = (*idum) / IQ1; //Start here when not initializing.
	*idum = IA1 * (*idum - k * IQ1) - k * IR1; //Compute idum=(IA1*idum) % IM1 without
	if (*idum < 0) *idum += IM1; //overflows by Schrage's method.
	k = idum2 / IQ2;
	idum2 = IA2 * (idum2 - k * IQ2) - k * IR2; //Compute idum2=(IA2*idum) % IM2 likewise.
	if (idum2 < 0) idum2 += IM2;
	j = iy / NDIV1; //Will be in the range 0..NTAB-1.
	iy = iv[j] - idum2; //Here idum is shued, idum and idum2 are
	iv[j] = *idum; //combined to generate output.
	if (iy < 1) iy += IMM1;
	if ((temp = AM1 * iy) > RNMX) return RNMX; //Because users don't expect endpoint values.
	else return temp;
}

double gasdev(long* idum)/*Returns a normally distributed deviate with zero mean and unit variance, using ran1(idum)
as the source of uniform deviates.*/ {
	static int iset = 0;
	static float gset;
	double fac, rsq, v1, v2;
	if (*idum < 0) iset = 0; //Reinitialize.
	if (iset == 0) { //We don't have an extra deviate handy, so
		do {
			v1 = 2.0 * ran2(idum) - 1.0; //pick two uniform numbers in the square ex
			v2 = 2.0 * ran2(idum) - 1.0; //tending from -1 to +1 in each direction,
			rsq = v1 * v1 + v2 * v2; //see if they are in the unit circle,
		} while (rsq >= 1.0 || rsq == 0.0); //and if they are not, try again.
		fac = sqrt(-2.0 * log(rsq) / rsq);
		/*Now make the Box-Muller transformation to get two normal deviates. Return one and
		save the other for next time.*/
		gset = v1 * fac;
		iset = 1;
		//Set flag.
		return v2 * fac;
	}
	else { //We have an extra deviate handy,
		iset = 0; //so unset the flag,
		return gset; //and return it.
	}
}


#define EPS_S 1.0e-20
//Here EPS is the absolute precision, which should be adjusted to the scale of your variables.
#define FREEALL free_ivector(l3,1,m);free_ivector(l1,1,n+1);

void simp1(double **a, int mm, int ll[], int nll, int iabf, int *kp, double *bmax){
	int k;
	double test;
	if (nll <= 0) *bmax=0.0;
	else {
		*kp=ll[1];
		*bmax=a[mm+1][*kp+1];
		for (k=2;k<=nll;k++) {
			if (iabf == 0)test=a[mm+1][ll[k]+1]-(*bmax);
			else
				test=fabs(a[mm+1][ll[k]+1])-fabs(*bmax);
			if (test > 0.0) {
				*bmax=a[mm+1][ll[k]+1];
				*kp=ll[k];
			}
		}
	}
}

void simp2(double **a, int n, int l2[], int nl2, int *ip, int kp, double *q1) {
	int k,ii, i;
	double qp,q0,q;
	*ip=0;
	for (i=1;i<=nl2;i++)
		if (a[i+1][kp+1] < -EPS_S) break;
	if (i>nl2) return;
	*q1 = -a[l2[i]+1][1]/a[l2[i]+1][kp+1];
	*ip=l2[i];
	for (i=*ip+1;i<=nl2;i++) {
		ii=l2[i];
		if (a[ii+1][kp+1] < -EPS_S) {
			q = -a[ii+1][1]/a[ii+1][kp+1];
			if (q < *q1) {
				*ip=ii;
				*q1=q;
			}
			else if (q == *q1) { //We have a degeneracy.
				for (k=1;k<=n;k++) {
					qp = -a[*ip+1][k+1]/a[*ip+1][kp+1];
					q0 = -a[ii+1][k+1]/a[ii+1][kp+1];
					if (q0 != qp) break;
				}
				if (q0 < qp) *ip=ii;
			}
		}
	}
}

void simp3(double **a, int i1, int k1, int ip, int kp)//Matrix operations to exchange a left-hand and right-hand variable (see text).
{
	int kk,ii;
	double piv;
	piv=1.0/a[ip+1][kp+1];
	for (ii=1;ii<=i1+1;ii++)
		if (ii-1 != ip) {
			a[ii][kp+1] *= piv;
			for (kk=1;kk<=k1+1;kk++)
				if (kk-1 != kp) a[ii][kk] -= a[ip+1][kk]*a[ii][kp+1];
		}
	for (kk=1;kk<=k1+1;kk++) if (kk-1 != kp) a[ip+1][kk] *= -piv;
	a[ip+1][kp+1]=piv;
}

void simplx(double **a, int m, int n, int m1, int m2, int m3, int *icase, int izrov[], int iposv[]) {
	int i,ip,ir,is,k,kh,kp,m12,nl1, nl2;
	int *l1,*l2,*l3;
	double q1,bmax;
	if (m != (m1+m2+m3)) {
		printf("Bad input constraint counts in simplx");
		return;
	}
	l1=ivector(1,n+1);
	l2=ivector(1,m);
	l3=ivector(1,m);
	nl1=n;
	for (k=1;k<=n;k++) l1[k]=izrov[k]=k;
	nl2=m;
	for (i=1;i<=m;i++) {
		if (a[i+1][1] < 0.0) {
			printf("Bad input tableau in simplx");
			return;
		}
		l2[i]=i;
		iposv[i]=n+i;
	}
	for (i=1;i<=m2;i++) l3[i]=1;
	ir=0;
	if (m2+m3) {
		ir=1;
		for (k=1;k<=(n+1);k++) {
			q1=0.0;
			for (i=m1+1;i<=m;i++) q1 += a[i+1][k];
			a[m+2][k] = -q1;
		}
		do {
			simp1(a,m+1,l1,nl1,0,&kp,&bmax);
			if(bmax <= EPS_S && a[m+2][1] < -EPS_S) {
				*icase = -1;
				FREEALL return;
			}
			else if (bmax <= EPS_S && a[m+2][1] <= EPS_S) {
				m12=m1+m2+1;
				for (ip=m12;ip<=m;ip++) {
					if (iposv[ip] == (ip+n)) {
						simp1(a,ip,l1,nl1,1,&kp,&bmax);
						if (bmax > 0.0) goto one;
					}
				}
				ir=0;
				--m12;
				for (i=m1+1;i<=m1+m2;i++) if (l3[i-m1] == 1) for(k=1;k<=n+1;k++) a[i+1][k] = -a[i+1][k];
				break;
			}
			simp2(a,n,l2, nl2, &ip,kp, &q1);
			if (ip == 0) { //Maximum of auxiliary objective function is unbounded, so no feasible solution exists.
				*icase = -1;
				FREEALL return;
			}
one:
			simp3(a,m+1,n,ip,kp);
			if (iposv[ip] >= (n+m1+m2+1)) {
				for (k=1;k<=nl1;k++)
					if (l1[k] == kp) break;
				--nl1;
				for (is=k;is<=nl1;is++) l1[is]=l1[is+1];
				++a[m+2][kp+1];
				for(i=1;i<=m+2;i++) a[i][kp+1]=-a[i][kp+1];
			}
			else {
				if (iposv[ip]>=(n+m1+1)) {
					kh=iposv[ip]-m1-n;
					if(l3[kh]) {
						l3[kh]=0;
						++a[m+2][kp+1];
						for (i=1;i<=m+2;i++)
							a[i][kp+1] = -a[i][kp+1];

					}
				}
			}
			is=izrov[kp];
			izrov[kp]=iposv[ip];
			iposv[ip]=is;
		}
		while(ir);
	}
	for (;;) {
		simp1(a,0,l1,nl1,0,&kp,&bmax);
		if (bmax <= 0.0) {
			*icase=0;
			FREEALL return;
		}
		simp2(a, n,l2, nl2, &ip,kp, &q1);
		if (ip == 0) {
			*icase=1;
			FREEALL return;
		}
		simp3(a,m,n,ip,kp);
		is=izrov[kp];
		izrov[kp]=iposv[ip];
		iposv[ip]=is;
	}
}

