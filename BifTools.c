#include "BifSysDepend.h"

#define min(x,y)    (x)>(y)?(y):(x)
#define max(x,y)    (x)<(y)?(y):(x)

static char rcsid[] = "$Id: BifTools.c,v 1.7 2019/11/30 03:44:14 tetsushi Exp $";

/*************** tools *********************/

int gauss(int n, double a[NNM][NNM], double h[NNM])
{
	int i, j;
	double sum;

	for (i = 0; i < n - 1; i++){
		int dummy;
		dummy = pivoting(a, i, n);
		dummy = elimination(a, i, n);
	}

	for (i = n - 1; i >= 0; i--){
		sum = 0.0;
		for ( j = n -1; j > i; j--) sum += a[i][j] * h[j];
		h[i] = (a[i][n] - sum) / a[i][i];
	}
	return 0;
}

int pivoting(double a[NNM][NNM], int k, int n)
{
	int i,j, posess = 0;
	int swap = 0;
	double max = 0.0;
	double v0[NNM],v1[NNM];
	double xx; 

	for (i = k; i < n; i++){
		xx = sqrt(a[i][k] * a[i][k]); 	
		if (max < xx) { 
			max = xx; 
			posess = i; 
		}
	}
	if (posess != k) swap = 1;
	for (i = 0; i < n + 1; i++){
		v0[i] = a[posess][i];
		v1[i] = a[k][i];
	}
	for (i = 0; i < n + 1; i++){
		a[posess][i] = v1[i];
		a[k][i] = v0[i];
	}
	return swap;
}


int elimination(double a[NNM][NNM], int k, int n)
{
	double mp[NNM][NNM];
	int i, j;

	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			if (i == j) mp[i][j] = 1.0;
			else mp[i][j] = 0.0;
		}
	}
	for ( i = k + 1; i < n; i++){
		if (fabs(a[k][k]) < DBL_EPSILON){
			return -1;
		}
		mp[i][k] = -a[i][k] / a[k][k];
	}
	product(mp,a,n);
	return 0;
}

/* mp * a -> a */
void product( double mp[NNM][NNM], double a[NNM][NNM], int n)
{
	int i, j, k;
	double c[NNM][NNM];

	for (i = 0; i < NNM; i++)
		for (j = 0; j < NNM; j++)
			c[i][j] = 0.0;

	for (i = 0; i < n; i++){
		for (j = 0; j < n + 1; j++){
			c[i][j] = 0.0;
			for (k = 0; k < n; k++){
				c[i][j] += mp[i][k] * a[k][j];
			}
		}
	}
	for (i = 0; i < n; i++){
		for (j = 0; j < n+1; j++){
			a[i][j] = c[i][j];
		}
	}
}

void product1( double a[NNM][NNM], double b[NNM][NNM], 
	double c[NNM][NNM], int n, int m, int l)
{
	int i, j, k;
	for (i = 0; i < n; i++){
		for (j = 0; j < l; j++){
			c[i][j] = 0.0;
			for (k = 0; k < m; k++){
				c[i][j] += a[i][k] * b[k][j];
			}
		}
	}
}


void matpromatmn( double a[NDE][NDE], double b[NDE][NDE], 
	double c[NDE][NDE], int n, int m, int l)
{
    int i, j, k;

    for (i = 0; i < n; i++){
        for (j = 0; j < l; j++){
            c[i][j] = 0.0;
            for (k = 0; k < m; k++){
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
}

/* a * b -> c */
void matpromat(double a[NDE][NDE], double b[NDE][NDE], 
	double c[NDE][NDE], int n)
{
    int i, j, k;

    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            c[i][j] = 0.0;
            for (k = 0; k < n; k++){
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
}

/* a -> c copy */
void matcopy(double a[NDE][NDE], double b[NDE][NDE], int n)
{
    int i, j;

    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
        b[i][j] = a[i][j];
        }
    }
}

/* a * b -> c */
void matprovec(double a[NDE][NDE], double b[NDE], 
	double c[NDE], int n)
{
    int i, j, k;

    for (i = 0; i < n; i++){
        c[i]= 0.0;
        for (k = 0; k < n; k++){ c[i] += a[i][k] * b[k]; }
    }
}

/* a * b -> c */
void vecpromat(double a[NDE], double b[NDE][NDE], 
	double c[NDE], int n, int m)
{
    int i, j;

    for (i = 0; i < m; i++){
        c[i] = 0.0;
        for (j = 0; j < n; j++){
            c[i] += a[j] * b[j][i];
        }
    }
}


/* double QR algorhithm by Dr. Ikeda */

void DQR( double a[NDE][NDE], double rel[NDE], double iml[NDE], int n)
{
	int	i, j, ns, nt;

	hessenberg(a);

	ns = 0;
	nt = NDE;
	for (j = 0; j < NDE; j++)
		rel[j] = iml[j] = 0.0;
	for (i = 0; i < 1000; i++) {
		for (j=nt-2; j>=ns; j--)
			if (fabs(a[j+1][j]) < DBL_EPSILON 
				* (fabs(a[j][j])+fabs(a[j+1][j+1]))) {
				if (j == nt-2) {
					singlel(a, rel, iml, j + 1);
				} else if (j == nt-3) {
					doublel(a, rel, iml, j+1);
				}
				nt = j+1;
				break;
			}
		if (nt <= ns+2) break;
		qrhessenberg(a,ns,nt);
	}
	if (nt==1)
		singlel(a, rel, iml, 0);
	else if (nt==2)
		doublel(a, rel, iml, 0);

}

void singlel(double a[NDE][NDE], 
	double re[NDE], double im[NDE], int j)
{
	re[j] = a[j][j];
	im[j] = 0.0;
}

void doublel(double a[NDE][NDE], 
	double re[NDE], double im[NDE], int j)
{
	double	b, c, x, y;

	b = -a[j][j]-a[j+1][j+1];
	c = a[j][j]*a[j+1][j+1]-a[j+1][j]*a[j][j+1];
	x = -b/2.0;
	y = sqrt(fabs(b * b - 4.0 * c))/2.0;
	re[j] = re[j+1] = x;
	if ( b * b - 4.0 * c < 0.0) {
		im[j] = y;
		im[j+1] = -y;
	} else {
		im[j] = 0.0;
		im[j+1] = 0.0;
		if (x < 0.0)
			re[j] -= y;
		else
			re[j] += y;
		re[j+1] = c/re[j];
	}
}

void hessenberg(double a[NDE][NDE])
{
	double	s,t,c,v[NDE];
	int	i,j,k;

	for (k = 0; k < NDE - 1; k++) {
		t = 0.0;
		for (i = k + 1; i < NDE; i++) {
			v[i] = a[i][k];
			t += v[i] * v[i];
			a[i][k] = 0.0;
		}
		s = sqrt(t);
		if (v[k+1] < 0.0) s = -s;
		v[k+1] += s;
		c = 1.0 / (v[k+1]*s);
		a[k+1][k] = -s;
		for (j = k + 1; j < NDE; j++) {
			t = 0.0;
			for (i = k + 1; i < NDE; i++)
				t += v[i] * a[i][j];
			t *= c;
			for (i = k + 1; i < NDE; i++)
				a[i][j] -= t * v[i];
		}
		for (i = 0; i < NDE; i++) {
			t = 0.0;
			for (j = k + 1; j < NDE; j++)
				t += a[i][j] * v[j];
			t *= c;
			for (j = k + 1; j < NDE; j++)
				a[i][j] -= t * v[j];
		}
	}
}

void qrhessenberg(double a[NDE][NDE], int ns, int nt)
{
	double	x,y,z,s,t,c;
	int	i,j,k,m,n;

	n=nt-1; m=n-1;
	s = a[m][m]+a[n][n];
	t = a[m][m]*a[n][n]-a[m][n]*a[n][m];

	n = ns; m = ns+1;
	x = a[n][n]*a[n][n]+a[n][m]*a[m][n]-s*a[n][n]+t;
	y = a[m][n]*(a[n][n]+a[m][m]-s);
	z = a[n+2][n+1]*a[n+1][n];

	n = nt-2;
	for (k=ns; k<nt-2; k++) {
		s = x*x+y*y+z*z;
		s = sqrt(s);
		if (x<0.0) s = -s;
		x += s;
		c = 1.0/(x*s);
		for (j=max(k-1,0); j<nt; j++) {
			s = c*(x*a[k][j]+y*a[k+1][j]+z*a[k+2][j]);
			a[k][j] -= s*x;
			a[k+1][j] -= s*y;
			a[k+2][j] -= s*z;
		}
		m = min(nt,k+4);
		for (i=ns; i<m; i++) {
			s = c*(a[i][k]*x+a[i][k+1]*y+a[i][k+2]*z);
			a[i][k] -= x*s;
			a[i][k+1] -= y*s;
			a[i][k+2] -= z*s;
		}
		x = a[k+1][k];
		y = a[k+2][k];
		if (k<nt-3) {
			z = a[k+3][k];
		}
	}
	s = x*x+y*y;
	s = sqrt(s);
	if (x<0.0) s = -s;
	x += s;
	c = 1.0/(x*s);
	for (j=max(k-1,0); j<nt; j++) {
		s = c*(x*a[k][j]+y*a[k+1][j]);
		a[k][j] -= s*x;
		a[k+1][j] -= s*y;
	}
	for (i=ns; i<nt; i++) {
		s = c*(a[i][k]*x+a[i][k+1]*y);
		a[i][k] -= x*s;
		a[i][k+1] -= y*s;
	}
}

/* end of Ikeda's contribution */


double inner(double a[NDE], double b[NDE], int n)
{
    int k;
    double c = 0.0;

    for (k = 0; k < n; k++){ c += a[k] * b[k]; }
    return c;
}

void printvec(double a[NDE], int n)
{
    int i;
    for (i = 0 ; i < n; i++){ printf("%lf\n", a[i]); }
    printf("\n");
}


/**********************************************************
  Calculating a determinant using LU decomposion. 
    	Reference: Introduction to Numerical Computations
    	J. S. Vandergraft, Academic Press, 1978
***********************************************************/ 
double det( double x[NDE][NDE], int n)
{
    double d =1.0;
    double q;
    double a[NNM][NNM];
    int count = 0;
    int result;
    int i, j, k;

	/* copy x -> a */
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            a[i][j] = x[i][j];
        }
    }
    /* a prosess for an extensional matrix */
    for (i = 0; i < n ; i++){  
    	a[i][n] = 0.0;
	}

    /* Using LU decomposition */
    for (i = 0; i < n - 1; i++){
    	count += pivoting(a, i, n);
    	result = elimination(a, i, n);
    	if (result < 0){ 
    		/* not a full rank */
			return 0.0;
		}
	}
    for (k = 0; k < n - 1; k++){
        d *= a[k][k];
        /* below routine is unneeded */
        /****
        for (i = k + 1; i < n ; i++){
            q = a[i][k] / a[k][k];
            for (j = k + 1; j < n; j++){
                a[i][j] -= (q * a[k][j]);
            }
        }
        *****/
    }
    d *= a[n-1][n-1];
    for (i = 0; i < count; i++){ d *= -1.0; }
    return d;
}


/* Expansion of a determinant with a complex eigenvalue. */


DecompMat *ExpandDet(double dx[][NDE], double RE, double IM)
{
	int i, j;
	int countup, fig, loop, found = 0;
	int loc[NDE];
	DecompMat *now_decompmat, *prev_decompmat, *root;

	root = NULL;
	for (fig = 0; fig < NDE + 1; fig++){
		for (countup = 0; countup < (1 << NDE); countup++){
			if (choosecomb(countup, fig, NDE, loc)) {
				for (loop = 0; loop < (1 << fig); loop++){
					now_decompmat = (DecompMat *)calloc(sizeof (DecompMat), 1);
					now_decompmat->next = NULL;
					now_decompmat->mu_real = RE;
					now_decompmat->mu_image = IM;
					if (countup == 0){
						root = now_decompmat;
					} else {
						prev_decompmat->next = now_decompmat;
					}
					now_decompmat->number = found;
					now_decompmat->n = NDE;
					now_decompmat->r = fig;
					switch (treatcomplex(loop, fig)){
						case ZERO:
							now_decompmat->sign = 1; 
							now_decompmat->has_image = 0;
						break;
						case POS_REAL:
							now_decompmat->sign = 1; 
							now_decompmat->has_image = 0;
						break;
						case POS_IMAGE:
							now_decompmat->sign = 1; 
							now_decompmat->has_image = 1;
						break;
						case NEG_REAL:
							now_decompmat->sign = -1; 
							now_decompmat->has_image = 0;
						break;
						case NEG_IMAGE:
							now_decompmat->sign = -1; 
							now_decompmat->has_image = 1;
						break;
					}
					for (j = 0; j < fig; j++){
						now_decompmat->locator[j].pos = loc[j];
						now_decompmat->locator[j].property = 
								(loop >> j) & IMAGEVAL; 
					}
					found++;
					prev_decompmat = now_decompmat;
				}
			}
		}
	}
	return root;
}


void SumDecompMat(DecompMat *decompmat, double *re, double *im)
{
	*re += ((decompmat->has_image) ? 0.0 : decompmat->ret_real ); 
	*im += ((decompmat->has_image) ? decompmat->ret_image : 0.0); 
	if (decompmat->next != NULL){
		SumDecompMat(decompmat->next, re, im);
	} 
}


void putbit(char *s, int n, int w)
{
	int i;

	for (i = 0; i < w; i++) s[w - 1 - i] = ((n >> i) & 1) + '0'; 
	s[w] = '\0';
}

int choosecomb(int k, int r, int n, int location[])
{
	unsigned int x;
	int i, count = 0;

	for (i = 0; i < n; i++){
		x = (k >> i) & 1; 
		if (x){
			 location[count] = i;
			 count++;
		}
	}
	if (count == r) return 1;
	else return 0;
}


int countcomb(int n, int r)
{ /* calculate nCn */
	int ret, fact(int);

	if (r == 0) ret = 1;
	else if (n == r) ret = 1;
	else ret = fact(n)/(fact(n-r) * fact(r));
	return ret;
}

int fact(int n)
{
	/* calculate n! */
	int i = n;
	
	while (--n){ i *= n; }
	return i;
}

int treatcomplex(int n, int w)
{	
	int countbit(int, int);
	/* returns 0 to 4 */
	/* 0 --- zero */
	/* 1 --- i^0 */
	/* 2 --- i^1 */
	/* 3 --- i^2 */
	/* 4 --- i^3 */

	if (n == 0) return 0;
	else return countbit(n, w) % 4 + 1;
}


int countbit(int s, int figures)
{
	int i, x, count = 0;

	for (i = 0; i < figures; i++){
		x = (s >> i) & 1;
		if (x) count++;
	}
	return count;
}

int print_decompmat_data(DecompMat *decompmat)
{
	int i;

	printf("--------\n");
	printf("No.%d \n", decompmat->number + 1);
	printf("combination = %dC%d\n", decompmat->n, decompmat->r);
	printf("sign = %s\n", (decompmat->sign < 0) ? "-" : "+" );
	printf("%s\n", (decompmat->has_image == 1) ? "IMAGE" : "REAL");
	printf("location : ");
	for (i = 0; i < decompmat->r; i++){
		printf("[%d] ", decompmat->locator[i].pos);
	}
	printf("\n");
	if (decompmat->next == NULL) {	
		return -1; 
	} else {
		print_decompmat_data(decompmat->next);
	}
	return 0;
}

void ComplexDet(double dx[][NDE], 
	DecompMat *decompmat, double *re, double *im) 
{
	int i, j;
	static int count = 0;
	double result;
	double a[NDE][NDE];

	count++;
	if (count == 1){ 
		*re = *im = 0.0;
	}

	matcopy(dx, a, NDE);
	for (i = 0; i < decompmat->r; i++){
		for (j = 0; j < NDE; j++){
			if (j == decompmat->locator[i].pos){
				a[j][decompmat->locator[i].pos] = 
					(decompmat->locator[i].property) ? 
						-decompmat->mu_image : -decompmat->mu_real;
			} else {
				a[j][decompmat->locator[i].pos] = 0.0;
			}
		}
	}
	result = det(a, NDE);
	result *= decompmat->sign;

	if (decompmat->has_image){
		*im += result;
	} else {
		*re += result;
	}

	if (decompmat->next != NULL){
		ComplexDet(dx, decompmat->next, re, im);
	} 
	count = 0;
	return;
}

void DerivDetByX(DecompMat *decompmat, 
	int varno, double *re, double *im, SysData *sys)
{
	int column, i, j;
	double *x;
	double a[NDE][NDE];
	double sum = 0.0;
	int row;
	double result;
	static int count  = 0;

	x = sys->x;
	sum = 0.0;

	count++;
	if (count == 1){  /* initialization */
		*re = 0.0; 
		*im = 0.0;
	}

	for (column = 0; column < NDE; column++){
		matcopy(sys->Dx, a, NDE);
		for (row = 0; row < NDE; row++){
			a[row][column] = 
				x[NE1 + NDE * NDE * varno + NDE * column + row];
		}
		for (i = 0; i < decompmat->r; i++){
			for (j = 0; j < NDE; j++){
				if (j == decompmat->locator[i].pos){
					a[j][decompmat->locator[i].pos] = 
						(decompmat->locator[i].property) ? 
							-decompmat->mu_image : -decompmat->mu_real;
				} else {
					a[j][decompmat->locator[i].pos] = 0.0;
				}
			}
		}
		for (i = 0; i < decompmat->r; i++){
			if (decompmat->locator[i].pos == column){ 
				for (j = 0; j < NDE; j++){
					a[j][decompmat->locator[i].pos] = 0.0;
				}
			}
		}
#ifdef DEBUG_X
	printf("varno = %d\n", varno);
	printf("NO = %d\n", decompmat->number);
	printf("(n, r) = (%d, %d)\n", decompmat->n, decompmat->r);
	printf("position =");
	for (i = 0 ; i < decompmat->r; i++){
		printf("[%d:%c] ", decompmat->locator[i].pos, 
			(decompmat->locator[i].property == 1) ? 'I' : 'R'); 
	}
	printf("\n");
		for (i = 0; i < NDE; i++){
			for (j = 0; j < NDE; j++){
				printf("%lf ", a[i][j] );
			}
			printf("\n");
		}
#endif
		result = det(a, NDE);
		result *= decompmat->sign;
		if (decompmat->has_image){
			*im += result;
		} else {
			*re += result;
		}
	}

	if (decompmat->next != NULL){
		DerivDetByX(decompmat->next, varno, re, im, sys);
	}
	count = 0;
	return;
}

void DerivDetByL(DecompMat *decompmat, 
	double *re, double *im, SysData *sys)
{
	int column, i, j;
	double *x;
	double a[NDE][NDE];
	double sum = 0.0;
	int row;
	double result;
	static int count  = 0;

	x = sys->x;
	sum = 0.0;

	count++;
	if (count == 1){  /* initialization */
		*re = 0.0; 
		*im = 0.0;
	}

	for (column = 0; column < NDE; column++){
		matcopy(sys->Dx, a, NDE);
		for (row = 0; row < NDE; row++){
			a[row][column] = 
				x[NE1 + NDE * NDE * NDE + NDE * column + row];
		}
		for (i = 0; i < decompmat->r; i++){
			for (j = 0; j < NDE; j++){
				if (j == decompmat->locator[i].pos){
					a[j][decompmat->locator[i].pos] = 
						(decompmat->locator[i].property) ? 
							-decompmat->mu_image : -decompmat->mu_real;
				} else {
					a[j][decompmat->locator[i].pos] = 0.0;
				}
			}
		}
		for (i = 0; i < decompmat->r; i++){
			if (decompmat->locator[i].pos == column){ 
				for (j = 0; j < NDE; j++){
					a[j][decompmat->locator[i].pos] = 0.0;
				}
			}
		}
#ifdef DEBUG_X
	printf("NO = %d\n", decompmat->number);
	printf("(n, r) = (%d, %d)\n", decompmat->n, decompmat->r);
	printf("position =");
	for (i = 0 ; i < decompmat->r; i++){
		printf("[%d:%c] ", decompmat->locator[i].pos, 
			(decompmat->locator[i].property == 1) ? 'I' : 'R'); 
	}
	printf("\n");
		for (i = 0; i < NDE; i++){
			for (j = 0; j < NDE; j++){
				printf("%lf ", a[i][j] );
			}
			printf("\n");
		}
#endif
		result = det(a, NDE);
		result *= decompmat->sign;
		if (decompmat->has_image){
			*im += result;
		} else {
			*re += result;
		}
	}

	if (decompmat->next != NULL){
		DerivDetByL(decompmat->next,  re, im, sys);
	}
	count = 0;
	return;
}


void DerivDetByTheta(DecompMat *decompmat, 
	double *re, double *im, SysData *sys)
{
	int loop, i, j;
	double *x;
	double a[NDE][NDE];
	double res = 0.0;
	double st, ct, dst, dct;
	static int count = 0;

	x = sys->x;
	res = 0.0;

	count++;
	if (count == 1){ 
		*re = 0.0; 
		*im = 0.0;
	}

	ct = decompmat->mu_real; 	/* cos \theta */
	st = decompmat->mu_image; 	/* sin \theta */

	dst = ct;				/* derivatives */
	dct = -st;

	for (loop = 0; loop < decompmat->r; loop++){
		matcopy(sys->Dx, a, NDE);
		for (i = 0; i < decompmat->r; i++){
			for (j = 0; j < NDE; j++){
				if (j == decompmat->locator[i].pos){
					a[j][decompmat->locator[i].pos] = 
						(decompmat->locator[i].property) ? 
							-st : -ct;
				} else {
					a[j][decompmat->locator[i].pos] = 0.0;
				}
			}
		}
		if (decompmat->locator[loop].property == IMAGEVAL){
			for (i = 0; i < NDE; i++){
				if (i == decompmat->locator[loop].pos){
					a[i][decompmat->locator[loop].pos] = -dst;
				} else {
					a[i][decompmat->locator[loop].pos] = 0.0;
				}
			}
		} else {
			for (i = 0; i < NDE; i++){
				if (i == decompmat->locator[loop].pos){
					a[i][decompmat->locator[loop].pos] = -dct;
				} else {
					a[i][decompmat->locator[loop].pos] = 0.0;
				}
			}
		}
#ifdef DEBUG_THETA
	printf("NO = %d\n", decompmat->number);
	printf("(n, r) = (%d, %d)\n", decompmat->n, decompmat->r);
	printf("position =");
	for (i = 0 ; i < decompmat->r; i++){
		printf("[%d:%c] ", decompmat->locator[i].pos, 
			(decompmat->locator[i].property == 1) ? 'I' : 'R'); 
	}
	printf("\n");
	for (i = 0; i < NDE; i++){
		for (j = 0; j < NDE; j++){
			printf("%lf ", a[i][j] );
		}
		printf("\n");
	}
	printf("RES = (%s)%.15f\n", 
		decompmat->has_image ? "image" : "real", det(a,NDE));
#endif
		res = det(a, NDE);
		res *= decompmat->sign;
		if (decompmat->has_image){
			*im += res;
		} else {
			*re += res;
		}
	}

	if (decompmat->next != NULL){
		DerivDetByTheta(decompmat->next, re, im, sys);
	}
	count = 0;
	return;
}



void DerivDetByTau(DecompMat *decompmat, double *f, 
	double *re, double *im, SysData *sys)
{
	int column, i, j;
	double a[NDE][NDE];
	double sum = 0.0;
	double result;
	int row;
	static int count = 0;

	count++;
	if (count == 1){ 
		*re = 0.0; 
		*im = 0.0;
	}
	sum = 0.0;

	for (column = 0; column < NDE; column++){
		matcopy(sys->Dx, a, NDE);
		for (row = 0; row < NDE; row++){
			a[row][column] = f[NDE + NDE * column + row];
		}
		for (i = 0; i < decompmat->r; i++){
			for (j = 0; j < NDE; j++){
				if (j == decompmat->locator[i].pos){
					a[j][decompmat->locator[i].pos] = 
						(decompmat->locator[i].property) ? 
							-decompmat->mu_image : -decompmat->mu_real;
				} else {
					a[j][decompmat->locator[i].pos] = 0.0;
				}
		}
		}
		for (i = 0; i < decompmat->r; i++){
				if (decompmat->locator[i].pos == column){ 
					for (j = 0; j < NDE; j++){
						a[j][decompmat->locator[i].pos] = 0.0;
					}
				}
		}
#ifdef DEBUG_TAU
		printf("differentiated by tau\n");
		printf("NO = %d\n", decompmat->number);
		printf("(n, r) = (%d, %d)\n", decompmat->n, decompmat->r);
		printf("position =");
		for (i = 0 ; i < decompmat->r; i++){
			printf("[%d:%c] ", decompmat->locator[i].pos, 
				(decompmat->locator[i].property == 1) ? 'I' : 'R'); 
		}
		printf("\n");
			for (i = 0; i < NDE; i++){
				for (j = 0; j < NDE; j++){
					printf("%lf ", a[i][j] );
				}
				printf("\n");
			}
#endif
		result = det(a, NDE);
		result *= decompmat->sign;
		if (decompmat->has_image){
			*im += result;
		} else {
			*re += result;
		}
	}

	if (decompmat->next != NULL){
		DerivDetByTau(decompmat->next, f, re, im, sys);
	}
	count = 0;
	return;
}


void FreeDecompMat(DecompMat *decompmat)
{
	if (decompmat->next != NULL){
		FreeDecompMat(decompmat->next);
	} 
	free(decompmat);
}

/* _____________________________________________________________________*/

/* Bifurcation of an Equilibrium */
/* derivative of the characteristic equaiton for the equilibria */
/* 2001 Feb. 24 */

void DerivDetEqByX(DecompMat *decompmat, int varno, 
	double *re, double *im, SysData *sys)
{
	int column, i, j;
	double *x;
	double a[NDE][NDE];
	double sum = 0.0;
	int row;
	double result;
	static int count  = 0;

	x = sys->x;
	sum = 0.0;

	count++;
	if (count == 1){  /* initialization */
		*re = 0.0; 
		*im = 0.0;
	}

	for (column = 0; column < NDE; column++){
		matcopy(sys->Df, a, NDE);
		for (row = 0; row < NDE; row++){
			a[row][column] = sys->DF2[varno][row][column];
		}
		for (i = 0; i < decompmat->r; i++){
			for (j = 0; j < NDE; j++){
				if (j == decompmat->locator[i].pos){
					a[j][decompmat->locator[i].pos] = 
						(decompmat->locator[i].property) ? 
							-decompmat->mu_image : -decompmat->mu_real;
				} else {
					a[j][decompmat->locator[i].pos] = 0.0;
				}
			}
		}
		for (i = 0; i < decompmat->r; i++){
			if (decompmat->locator[i].pos == column){ 
				for (j = 0; j < NDE; j++){
					a[j][decompmat->locator[i].pos] = 0.0;
				}
			}
		}
#ifdef DEBUG_EQ_X
	printf("varno = %d\n", varno);
	printf("NO = %d\n", decompmat->number);
	printf("(n, r) = (%d, %d)\n", decompmat->n, decompmat->r);
	printf("position =");
	for (i = 0 ; i < decompmat->r; i++){
		printf("[%d:%c] ", decompmat->locator[i].pos, 
			(decompmat->locator[i].property == 1) ? 'I' : 'R'); 
	}
	printf("\n");
		for (i = 0; i < NDE; i++){
			for (j = 0; j < NDE; j++){
				printf("%lf ", a[i][j] );
			}
			printf("\n");
		}
#endif
		result = det(a, NDE);
		result *= decompmat->sign;
		if (decompmat->has_image){
			*im += result;
		} else {
			*re += result;
		}
	}

	if (decompmat->next != NULL){
		DerivDetEqByX(decompmat->next, varno, re, im, sys);
	}
	count = 0;
	return;
}

void DerivDetEqByL(DecompMat *decompmat, double *re, double *im, SysData *sys)
{
	int column, i, j;
	double *x;
	double a[NDE][NDE];
	double sum = 0.0;
	int row;
	double result;
	static int count  = 0;

	x = sys->x;
	sum = 0.0;

	count++;
	if (count == 1){  /* initialization */
		*re = 0.0; 
		*im = 0.0;
	}

	for (column = 0; column < NDE; column++){
		matcopy(sys->Df, a, NDE);
		for (row = 0; row < NDE; row++){
			a[row][column] = sys->DFL[row][column];
		}
		for (i = 0; i < decompmat->r; i++){
			for (j = 0; j < NDE; j++){
				if (j == decompmat->locator[i].pos){
					a[j][decompmat->locator[i].pos] = 
						(decompmat->locator[i].property) ? 
							-decompmat->mu_image : -decompmat->mu_real;
				} else {
					a[j][decompmat->locator[i].pos] = 0.0;
				}
			}
		}
		for (i = 0; i < decompmat->r; i++){
			if (decompmat->locator[i].pos == column){ 
				for (j = 0; j < NDE; j++){
					a[j][decompmat->locator[i].pos] = 0.0;
				}
			}
		}
#ifdef DEBUG_EQ_L
	printf("NO = %d\n", decompmat->number);
	printf("(n, r) = (%d, %d)\n", decompmat->n, decompmat->r);
	printf("position =");
	for (i = 0 ; i < decompmat->r; i++){
		printf("[%d:%c] ", decompmat->locator[i].pos, 
			(decompmat->locator[i].property == 1) ? 'I' : 'R'); 
	}
	printf("\n");
		for (i = 0; i < NDE; i++){
			for (j = 0; j < NDE; j++){
				printf("%lf ", a[i][j] );
			}
			printf("\n");
		}
#endif
		result = det(a, NDE);
		result *= decompmat->sign;
		if (decompmat->has_image){
			*im += result;
		} else {
			*re += result;
		}
	}

	if (decompmat->next != NULL){
		DerivDetEqByL(decompmat->next, re, im, sys);
	}
	count = 0;
	return;
}


void inverse(double inmat[][NDE], double *x, int n)
{
	int  nn = NDE, mm = NDE;     
	int m = n;
	int     i,j,k,kk,i1,k1,imax,jmax ;
	double  s,p,q,r ;
	double  *ai,*aj,*ak,*aik,*aij ;
	double  *bi,*bj,*bk,*bik,*bij ;
	double  *xi,*xj,*xk,*xik,*xij ;
	double aa[NDE][NDE];
	double b[NDE][NDE];
	int	*jun;
	double	*sf;
	double *a;

	a = &aa[0][0];

	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			aa[i][j] = inmat[i][j];
			b[i][j] = (i == j) ? 1.0 : 0.0; 				
		}
	}

	jun=(int *)malloc(4*n);
	
	sf=(double *)malloc(8*n);
	/* balancing */

	for ( i=0 , ai=a , bi=&b[0][0] , xi=x ; i<n ;
	      ++i , ai+=nn , bi+=mm , xi+=mm )
		for ( j=0 ; j<m ; ++j )
			xi[j]=bi[j];

	for ( i=0 , ai=a , bi=&b[0][0] , xi=x ; i<n ;
	      ++i , ai+=nn , bi+=mm , xi+=mm )
	{
		p=0.0;
		for ( j=0 ; j<n ; ++j )
		{
			if ( (q=fabs(ai[j]))>p )
				p=q;	
		}
		for ( j=0 ; j<n ; ++j )
			ai[j]/=p;
		for ( j=0 ; j<m ; ++j )
			xi[j]=bi[j]/p;
	}
	for ( j=0 , aj=a ; j<n ; ++j , ++aj )
	{
		p=0.0;
		for ( i=0 , aij=aj ; i<n ; ++i , aij+=nn )
		{
			if ( (q=fabs(*aij))>p )
				p=q;	
		}
		for ( i=0 , aij=aj ; i<n ; ++i , aij+=nn )
			(*aij)/=p ;
		sf[j]=p;
	}

	for ( k=0 , ak=a   , xk=x ; k<n-1 ;
	      ++k , ak+=nn , xk+=mm )
	{
		k1=k+1 ;
		/* complete pivotting */
		imax=k;
		jmax=k;
		p=0.0;
		for ( i=k , ai=ak ; i<n ; ++i , ai+=nn )
		{
			for ( j=k ; j<n ; ++j )
			{
				if ( (q=fabs(ai[j]))>p )
				{
					p=q;	
					imax=i;
					jmax=j;
				}
			}
		}
		jun[k]=jmax;
		if ( imax!=k ) 
		{
		    ai=a+imax*nn;
		    for ( j=k ; j<n ; ++j ) 
		    {
			s=ak[j];
			ak[j]=ai[j];
			ai[j]=s;
		    }
		    xi=x+imax*mm;
		    for ( j=0 ; j<m ; ++j ) 
		    {
			s=xk[j];
			xk[j]=xi[j];
			xi[j]=s;
		    }
		}
		if ( jmax!=k )
		{
		    for ( i=0 , aik=a+k , aij=a+jmax ; i<n ;
		          ++i , aik+=nn , aij+=nn )
		    {
			s=(*aik);
			(*aik)=(*aij);
			(*aij)=s;
		    }
		}
		/*  k  */
		p=ak[k] ;
		for ( j=k1 ; j<n ; ++j )
			ak[j]/=p ;
		for ( j=0 ; j<m ; ++j )
			xk[j]/=p ;
		/*  i  */
		for ( i=k1 , ai=ak+nn , xi=xk+mm ; i<n ;
		       ++i , ai+=nn   , xi+=mm )
		{
			q=ai[k] ;
			for ( j=k1 ; j<n ; ++j )
				ai[j]-=q*ak[j] ;
			for ( j=0 ; j<m ; ++j )
				xi[j]-=q*xk[j] ;
		} 
	} 
	p=ak[n-1] ;
	xi=x+(n-1)*mm ;
	for ( j=0 ; j<m ; ++j )
		xi[j]/=p ;
	for ( k=n-1 ; k>0 ; --k )
	{
		xk=x+k*mm ;
		for ( i=0 ; i<k ; ++i ) 
		{
			xi=x+i*mm ;
			q=(*(a+nn*i+k)) ;
			for ( j=0 ; j<m ; ++j ) 
				xi[j]-=q*xk[j] ;
		}
	}

	for ( k=n-2 ; k>=0 ; --k )
	{
		i=jun[k];
		xi=x+i*mm;
		xk=x+k*mm;
		for ( j=0 ; j<m ; ++j ) 
		{
			s=xi[j];
			xi[j]=xk[j];
			xk[j]=s;
		}
	}
	for ( i=0 , xi=x ; i<n ; ++i , xi+=mm )
	{
		p=sf[i];
		for ( j=0 ; j<m ; ++j )
		{
			xi[j]/=p;
		}
	}
	free(jun);
	free(sf);
}

void DerivDetEqByOmega(DecompMat *decompmat, 
	double *re, double *im, SysData *sys)
{
	int loop, i, j;
	double *x;
	double a[NDE][NDE];
	double res = 0.0;
	double st, ct, dst, dct;
	static int count = 0;

	x = sys->x;
	res = 0.0;

	count++;
	if (count == 1){ 
		*re = 0.0; 
		*im = 0.0;
	}

	for (loop = 0; loop < decompmat->r; loop++){
		matcopy(sys->Df, a, NDE);
		for (i = 0; i < decompmat->r; i++){
			for (j = 0; j < NDE; j++){
				if (j == decompmat->locator[i].pos){
					a[j][decompmat->locator[i].pos] = 
						(decompmat->locator[i].property) ? 
							-sys->omega : -sys->sigma;
				} else {
					a[j][decompmat->locator[i].pos] = 0.0;
				}
			}
		}
		if (decompmat->locator[loop].property == IMAGEVAL){
			for (i = 0; i < NDE; i++){
				if (i == decompmat->locator[loop].pos){
					a[i][decompmat->locator[loop].pos] = -1.0;
				} else {
					a[i][decompmat->locator[loop].pos] = 0.0;
				}
			}
		} else {
			for (i = 0; i < NDE; i++){
				if (i == decompmat->locator[loop].pos){
					a[i][decompmat->locator[loop].pos] = 0.0;
				} else {
					a[i][decompmat->locator[loop].pos] = 0.0;
				}
			}
		}

/*
	printf("NO = %d\n", decompmat->number);
	printf("(n, r) = (%d, %d)\n", decompmat->n, decompmat->r);
	printf("position =");
	for (i = 0 ; i < decompmat->r; i++){
		printf("[%d:%c] ", decompmat->locator[i].pos, 
			(decompmat->locator[i].property == 1) ? 'I' : 'R'); 
	}
	printf("\n");
	for (i = 0; i < NDE; i++){
		for (j = 0; j < NDE; j++){
			printf("%lf ", a[i][j] );
		}
		printf("\n");
	}
	printf("RES = (%s)%.15f\n", 
		decompmat->has_image ? "image" : "real", det(a,NDE));
*/

		res = det(a, NDE);
		res *= decompmat->sign;
		if (decompmat->has_image){
			*im += res;
		} else {
			*re += res;
		}
	}

	if (decompmat->next != NULL){
		DerivDetEqByOmega(decompmat->next, re, im, sys);
	}
	count = 0;
	return;
}


void runge(int n, double h, double x[], double t,  SysData *Dat)
{
    double k1[NVALL], k2[NVALL], k3[NVALL], k4[NVALL];
    double xtemp[NVALL];
    int i;

    function(x, k1, t, Dat);
    for (i = 0; i < n; i++) xtemp[i] = x[i] + k1[i] / 2.0 * h;
    t += h / 2.0;

    function(xtemp, k2, t, Dat);
    for (i = 0; i < n; i++) xtemp[i] = x[i] + k2[i] / 2.0 * h;

    function(xtemp, k3, t, Dat);
    for (i = 0; i < n; i++) xtemp[i] = x[i] + k3[i] * h;
    t += h / 2.0;

    function(xtemp, k4, t, Dat);
    for (i = 0; i < n; i++) x[i] += 
		(k1[i] + 2.0 * (k2[i] + k3[i]) + k4[i]) / 6.0 * h;

}
