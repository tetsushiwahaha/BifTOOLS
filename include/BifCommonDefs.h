/* $Id: BifCommonDefs.h,v 1.2 2019/11/22 09:24:19 tetsushi Exp tetsushi $ */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>

#define EPS 1.0e-8
#define NVE1 ((NDE) * (NDE))
#define NDVE1 NDE + NVE1 
#define NVEP1 NDE
#define NVE2 ((NDE) * (NDE) * (NDE))
#define NVEP2 ((NDE) * (NDE))
#define NE1 ((NDE) + (NVE1) + (NVEP1))
#define NVALL ((NE1) + (NVE2) + (NVEP2))
#define NNM 30
#define BUF 50
#define PMAX 20
#define sign(x) ((x) < 0.0 ? -1.0 : 1.0 )
#define LE_IDLING 500
#define LE_ITERATION 1000

#define ROW_Q ((NDE) - 1)
#define ROW_CHIR (NDE)
#define ROW_CHII ((NDE) + 1)
#define CLM_T ((NDE)-1)
#define CLM_L ((NDE))
#define CLM_TH ((NDE) + 1)

#define MODE_F 0
#define MODE_D 1
#define MODE_T 2
#define MODE_N 3
#define MODE_P 4
#define MODE_E 5
#define MODE_H 6
#define MODE_ET 7
#define MODE_L 8
#define MODE_EP 9

#define REALVAL 0
#define IMAGEVAL 1

#define ZERO      0
#define POS_REAL  1
#define POS_IMAGE 2
#define NEG_REAL  3
#define NEG_IMAGE 4


typedef struct { 
	double param[NP];
	double dparam[NP];

	double tau;

	int data;			/* number of data */

	double h;				/* Integral strip */
	int l;				/* period */
	int m;				/* number of strips */
	int section;

	double flag;	/* file write, draw graph flags */

	double x[NVALL];
	double xx[NVALL];
	double xe;
	double xn[PMAX][NDE];
	double Df[NDE][NDE];
	double Dx[NDE][NDE];
	double Dl[NDE];
	double DF2[NDE][NDE][NDE];
	double DFL[NDE][NDE];
	double P[NDE][NDE]; /* for pitchfork bif. */

	double LE[NDE]; /* Laypunov exponent */
	double Tau;
	double gamma;

	double emax;
	double eps0;
	double eps1;

	double mu;
	double sigma;
	double omega;
	double theta;
	double sign_of_sigma;

	int on_real_axis;

	int ite_max;
	int variable;
	int increment;
	int mode;

	int outflag;
	FILE *fpout, *fpout0;

} SysData;

void StoreDf(double [], double [], SysData *);
void function(double [], double [], double, SysData *);
void CalcEquil(SysData *);
void runge( int, double, double [], double, SysData *);

int fixed(SysData *);
void sysvar(SysData *);
void newton(SysData *);

double det(double [][NDE], int);
double Det(double [][NDE], int);
double inner(double [], double [], int);

/* a keyword 'Complex' cannot be used in MacOS */
typedef struct _ComplexV {
	double real;
	double image;
} ComplexV;

typedef struct _Locator {
	int pos;
	int property;
} Locator;

typedef struct _DecompMat {
	int number;
	int n;
	int r;
	int sign;
	int has_image;
	struct _DecompMat *next;
	double mu_real;
	double mu_image;
	double ret_real;
	double ret_image;
	struct _Locator locator[NDE];
} DecompMat;

DecompMat *ExpandDet(double [][NDE], double, double);
void ComplexDet(double [][NDE], DecompMat *, double *, double *); 
void SumDecompMat(DecompMat *, double *, double *);
void DerivDetByX(DecompMat *, int, double *, double *, SysData *);
void DerivDetByL(DecompMat *, double *, double *, SysData *);
void DerivDetByTheta(DecompMat *, double *, double *, SysData *);
void DerivDetByTau(DecompMat *, double *, double *, double *, SysData *);
void DerivDetEqByX(DecompMat *, int, double *, double *, SysData *);
void DerivDetEqByL(DecompMat *, double *, double *, SysData *);
void DerivDetEqByOmega(DecompMat *, double *, double *, SysData *);

void FreeDecompMat(DecompMat *);

void utox(double u[], double x[], SysData *);
void xtou(double x[], double u[], SysData *);
void conv_param_real_to_model(SysData *);
void conv_param_model_to_real(SysData *);
void show_param(SysData *);
double tautot(double tau, SysData *s);
double ttotau(double t, SysData *s);

SysData *BifInit(int, char **);
void inverse(double [][NDE], double *, int); 


void DQR( double [][NDE], double [NDE], double [NDE], int);
int pivoting( double [][NNM], int, int);
int elimination( double [][NNM], int, int);
void product(double [][NNM], double [][NNM], int);
void product1(double [][NNM], double [][NNM], double [][NNM], int, int, int);
void matpromatmn( double [][NDE], double [][NDE], 
	double [][NDE], int, int, int);
void matprovec(double [][NDE], double [], double [], int);
void vecpromat(double [], double [][NDE], double [], int, int);
void matcopy(double [][NDE], double [][NDE], int);
void hessenberg(double [][NDE]);
void qrhessenberg(double [][NDE], int, int);

void singlel(double [][NDE], double [], double [], int);
void doublel(double [][NDE], double [], double [], int);
void printvec(double [], int);
int choosecomb(int, int, int, int []);
int treatcomplex(int, int);
void putbit(char *, int, int);
int gauss(int, double [][NNM], double []);
