#include "BifSysDepend.h"

static char rcsid[] = "$Id: BifFunc.c,v 1.2 2006/08/25 03:39:08 tetsushi Exp $";

void function(double x[], double f[], double t, SysData *sys)
{
	int i, j, k;
	int pos;
	double Dx[NDE][NDE];
	double DF2[NDE][NDE][NDE], Dfx2[NDE][NDE][NDE];
	double Dl[NDE];
	double DFL[NDE][NDE], Dfl[NDE][NDE];
	double buf[NDE], res[NDE];

	/* LOCAL DEFINITIONS */
	/* ############################################################ */
	double sech(double);
	double k1, k2, delta, gamma;
	/* ############################################################ */

	/* initialize */

    for (i = 0; i < NDE ; i++){
        for (j = 0; j < NDE; j++){
            for (k = 0; k < NDE; k++){
                DF2[i][j][k] = 0.0;
                Dfx2[i][j][k] = 0.0;
            }
        }
    }
    for (i = 0; i < NDE ; i++){
        for (j = 0; j < NDE; j++){
            DFL[i][j] = 0.0;
            Dfl[i][j] = 0.0;
        }
    }
    for (i = 0; i < NDE; i++) Dl[i] = 0.0;

	/* ############################################################ */

	/* store local variables */

	gamma = sys->gamma = 1.6369909;
	k1 = sys->param[K1];
	k2 = sys->param[K2];
	delta = sys->param[DELTA];

	/* store system equation here */

    f[0] =  -x[1] + tanh(gamma * x[0]) - delta * (x[0] - x[2]);
	f[1] =  x[0] - k1 * x[1];
   	f[2] =  -x[3] + tanh(gamma * x[2]) - delta * (x[2] - x[0]);
   	f[3] =  x[2] - k2 * x[3];

	/* STORE JACOBI MATRIX HERE */

    sys->Df[0][0] = -delta + gamma * sech(gamma * x[0]) * sech(gamma * x[0]);
    sys->Df[0][1] = -1.0;
    sys->Df[0][2] =  delta;
    sys->Df[0][3] = 0.0;

    sys->Df[1][0] = 1.0;
    sys->Df[1][1] = -k1;
    sys->Df[1][2] = 0.0;
    sys->Df[1][3] = 0.0;

    sys->Df[2][0] = delta;
    sys->Df[2][1] =  0.0;
    sys->Df[2][2] = -delta + gamma * sech(gamma * x[2]) * sech(gamma * x[2]);
    sys->Df[2][3] = -1.0;

    sys->Df[3][0] = 0.0;
    sys->Df[3][1] = 0.0;
    sys->Df[3][2] = 1.0;
    sys->Df[3][3] = -k2;

	/* STORE Second Derivative of Jacobi Matrix */

	DF2[0][0][0] = -2.0 * gamma * gamma 
		* sech(gamma * x[0])
		* sech(gamma * x[0])
		* tanh(gamma * x[0]);
	DF2[2][2][2] = -2.0 * gamma * gamma 
		* sech(gamma * x[2])
		* sech(gamma * x[2])
		* tanh(gamma * x[2]);

	/* STORE variations related to parameters */

	switch(sys->variable){
		case DELTA:
			Dl[0] = -(x[0] - x[2]);
			Dl[1] = 0.0;
			Dl[2] = -(x[2] - x[0]);
			Dl[3] = 0.0;

			DFL[0][0] =  -1.0;
			DFL[0][2] =  1.0;
			DFL[2][0] =  1.0;
			DFL[2][2] =  -1.0;
		break;
		case K1:
			Dl[1] = -x[1];
			DFL[1][1] = -1.0;
		break;
		case K2:
			Dl[3] = -x[3];
			DFL[3][3] = -1.0;
		break;

	}

	/* ############################################################ */

	/* DO NOT EDIT BELOW */

	/* check routine the 2nd. variations (symmetry) */ 
	for (i = 0; i < NDE; i++){
		for (j = 0; j < NDE; j++){
			for (k = j + 1; k < NDE; k++){
				double diff;
				diff = fabs(DF2[i][j][k] - DF2[i][k][j]);
				if (diff > 1e-8){
					printf("something wrong with 2nd variations\n");
					printf("i = %d, j = %d, k = %d : val = %.15f\n", 
						i, j, k, diff);
				}
			}
		}
	}

	/* complete 1st-order variational eqs for states */
	/* i = NDE means variations for the parameter    */

	for (pos = NDE, i = 0; i < NDE + 1; i++, pos += NDE){
		for (j = 0; j < NDE; j++){ buf[j] = x[pos + j]; }
		matprovec(sys->Df, buf, res, NDE);
		for (j = 0; j < NDE; j++){ f[pos + j] = res[j]; }
	}

	/* complete 1st variational eqs for the parameter */

    for (j = 0; j < NDE; j++){ f[NDE + NDE * NDE + j] += Dl[j]; }

	/* in the following, prepare information for 2nd variational eqs. */

    for (i = 0; i < (NDE + 1); i++){
        for (j = 0; j < NDE; j++){
            for (k = 0; k < NDE; k++){ 
                buf[k] = x[NE1 + NDE * NDE * i + NDE * j + k]; 
            }
            matprovec(sys->Df, buf, res, NDE);
            for (k = 0; k < NDE; k++){ 
                f[NE1 + NDE * NDE * i + NDE * j + k] = res[k]; 
            }
        }
    }

    for (i = 0; i < NDE; i++){
        for (j = 0; j < NDE; j++){ buf[j] = x[NDE + NDE * i + j]; }
        for (j = 0; j < NDE; j++){ 
            for (k = 0; k < NDE; k++){
                Dfx2[i][j][k] = inner(DF2[j][k], buf, NDE);
            }
        }
    }

    for (i = 0; i < NDE; i++){
        for (j = 0; j < NDE; j++){
            buf[j] = x[NDE + NDE * i + j]; 
        }
        for (j = 0; j < NDE; j++){
            Dfl[i][j] = inner(DFL[j], buf, NDE);
        }
    }

    /* store complete 2nd variations. */
    /* attention the subscription! */
    /* subscript i and j is commutative (Jul.31, 1998)*/

    for (i = 0; i < NDE; i++){
        for (j = 0; j < NDE; j++){
            /* remark the subscription */
            for (k = 0; k < NDE; k++){ buf[k] = x[NDE + NDE * j + k]; }
            matprovec(Dfx2[i], buf, res, NDE);
            for (k = 0; k < NDE; k++){
                f[NE1 + NDE * NDE * i + NDE * j + k ] += res[k];
                /*
                f[NE1 + NDE * NDE * j + NDE * i + k ] += res[k];
                */
            }
        }
    }

    /* complete 2nd variations for parameters (1)*/

    for (i = 0; i < NDE; i++){
        for (j = 0; j < NDE; j++){ buf[j] = x[NDE + NDE * NDE + j]; }
        matprovec(Dfx2[i], buf, res, NDE);
        for (j = 0; j < NDE; j++){
            f[NE1 + NDE * NDE * NDE + NDE * i + j] += res[j];
        }
    }

    /* store also 2nd variations for the parameter (2)*/

    for (i = 0; i < NDE; i++){
        for (j = 0; j < NDE; j++){
            f[NE1 + NDE * NDE * NDE + NDE * i + j] += Dfl[i][j];
		}
	}
    for (i = 0; i < NDE; i++){
        sys->Dl[i] = Dl[i];
    }

    for (i = 0; i < NDE; i++){
        for (j = 0; j < NDE; j++){
            sys->DFL[i][j] = DFL[i][j];
        }
    }
}

double sech(double x)
{
    return 1.0 / cosh(x);
}

