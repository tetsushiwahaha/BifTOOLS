#include <BifSysDepend.h>

void function(double x[], double f[], double t, SysData *sys)
{
	int i, j, k;
	int pos;
	double Dx[NDE][NDE];
    double Dfx2[NDE][NDE][NDE];
	double Dl[NDE];
	double Dfl[NDE][NDE];
	double buf[NDE], res[NDE];



	double a, b, c, d, rx, ry, delta;
	double q[NDE], qd[NDE], qdd[NDE];
	void StoreDf(double *, double *, SysData *);



	/* initialize */

    for (i = 0; i < NDE ; i++){
        for (j = 0; j < NDE; j++){
            for (k = 0; k < NDE; k++){
                sys->DF2[i][j][k] = 0.0;
				Dfx2[i][j][k] = 0.0;
            }
        }
    }
    for (i = 0; i < NDE ; i++){
        for (j = 0; j < NDE; j++){
            sys->DFL[i][j] = 0.0;
			Dfl[i][j] = 0.0;
        }
    }


    for (i = 0; i < NDE; i++) { sys->Dl[i] = 0.0; }

	/* ############################################################ */

	/* store local variables */

	a = sys->param[A];
	b = sys->param[B];
	c = sys->param[C];
	d = sys->param[D];
	rx = sys->param[RX];
	ry = sys->param[RY];
	delta = sys->param[DELTA];





	StoreDf(x, q, sys);

	for (i = 0; i < NDE; i++){ qd[i] = (1.0 - q[i]) * q[i]; }
	for (i = 0; i < NDE; i++){ 
		qdd[i] = q[i] * (1.0 - q[i]) * (1.0 - 2.0 * q[i]);
	}



	sys->DF2[0][0][0] =  qdd[0] * a * a;
	sys->DF2[0][0][1] = -qdd[0] * a * b;
	sys->DF2[0][0][2] =  qdd[0] * a * delta;
	sys->DF2[0][0][3] =  0.0;
	sys->DF2[0][0][4] =  0.0;
	sys->DF2[0][0][5] =  0.0;
	sys->DF2[0][0][6] =  qdd[0] * a * delta;
	sys->DF2[0][0][7] =  0.0;

	sys->DF2[0][1][0] = -qdd[0] * a * b;
	sys->DF2[0][1][1] =  qdd[0] * b * b;
	sys->DF2[0][1][2] = -qdd[0] * b * delta;
	sys->DF2[0][1][3] =  0.0;
	sys->DF2[0][1][4] =  0.0;
	sys->DF2[0][1][5] =  0.0;
	sys->DF2[0][1][6] = -qdd[0] * b * delta;
	sys->DF2[0][1][7] =  0.0;

	sys->DF2[0][2][0] =  qdd[0] * a * delta;
	sys->DF2[0][2][1] = -qdd[0] * b * delta;
	sys->DF2[0][2][2] =  qdd[0] * delta * delta;
	sys->DF2[0][2][3] =  0.0;
	sys->DF2[0][2][4] =  0.0;
	sys->DF2[0][2][5] =  0.0;
	sys->DF2[0][2][6] =  qdd[0] * delta * delta;
	sys->DF2[0][2][7] =  0.0;

	sys->DF2[0][3][0] =  0.0;
	sys->DF2[0][3][1] =  0.0;
	sys->DF2[0][3][2] =  0.0;
	sys->DF2[0][3][3] =  0.0;
	sys->DF2[0][3][4] =  0.0;
	sys->DF2[0][3][5] =  0.0;
	sys->DF2[0][3][6] =  0.0;
	sys->DF2[0][3][7] =  0.0;

	sys->DF2[0][4][0] =  0.0;
	sys->DF2[0][4][1] =  0.0;
	sys->DF2[0][4][2] =  0.0;
	sys->DF2[0][4][3] =  0.0;
	sys->DF2[0][4][4] =  0.0;
	sys->DF2[0][4][5] =  0.0;
	sys->DF2[0][4][6] =  0.0;
	sys->DF2[0][4][7] =  0.0;

	sys->DF2[0][5][0] =  0.0;
	sys->DF2[0][5][1] =  0.0;
	sys->DF2[0][5][2] =  0.0;
	sys->DF2[0][5][3] =  0.0;
	sys->DF2[0][5][4] =  0.0;
	sys->DF2[0][5][5] =  0.0;
	sys->DF2[0][5][6] =  0.0;
	sys->DF2[0][5][7] =  0.0;

	sys->DF2[0][6][0] =  qdd[0] * a * delta;
	sys->DF2[0][6][1] = -qdd[0] * b * delta;
	sys->DF2[0][6][2] =  qdd[0] * delta * delta;
	sys->DF2[0][6][3] =  0.0;
	sys->DF2[0][6][4] =  0.0;
	sys->DF2[0][6][5] =  0.0;
	sys->DF2[0][6][6] =  qdd[0] * delta * delta;
	sys->DF2[0][6][7] =  0.0;

	sys->DF2[0][7][0] =  0.0;
	sys->DF2[0][7][1] =  0.0;
	sys->DF2[0][7][2] =  0.0;
	sys->DF2[0][7][3] =  0.0;
	sys->DF2[0][7][4] =  0.0;
	sys->DF2[0][7][5] =  0.0;
	sys->DF2[0][7][6] =  0.0;
	sys->DF2[0][7][7] =  0.0;

	sys->DF2[1][0][0] =  qdd[1] * c * c;
	sys->DF2[1][0][1] = -qdd[1] * c * d;
	sys->DF2[1][0][2] = 0.0;
	sys->DF2[1][0][3] = -qdd[1] * c * delta;
	sys->DF2[1][0][4] = 0.0;
	sys->DF2[1][0][5] = 0.0;
	sys->DF2[1][0][6] = 0.0;
	sys->DF2[1][0][7] = -qdd[1] * c * delta;

	sys->DF2[1][1][0] = -qdd[1] * c * d;
	sys->DF2[1][1][1] =  qdd[1] * d * d;
	sys->DF2[1][1][2] = 0.0;
	sys->DF2[1][1][3] = qdd[1] * d * delta;
	sys->DF2[1][1][4] = 0.0;
	sys->DF2[1][1][5] = 0.0;
	sys->DF2[1][1][6] = 0.0;
	sys->DF2[1][1][7] = qdd[1] * d * delta;

	sys->DF2[1][2][0] = 0.0;
	sys->DF2[1][2][1] = 0.0;
	sys->DF2[1][2][2] = 0.0;
	sys->DF2[1][2][3] = 0.0;
	sys->DF2[1][2][4] = 0.0;
	sys->DF2[1][2][5] = 0.0;
	sys->DF2[1][2][6] = 0.0;
	sys->DF2[1][2][7] = 0.0;

	sys->DF2[1][3][0] = -qdd[1] * c * delta;
	sys->DF2[1][3][1] =  qdd[1] * d * delta;
	sys->DF2[1][3][2] = 0.0;
	sys->DF2[1][3][3] = qdd[1] * delta * delta;
	sys->DF2[1][3][4] = 0.0;
	sys->DF2[1][3][5] = 0.0;
	sys->DF2[1][3][6] = 0.0;
	sys->DF2[1][3][7] = qdd[1] * delta * delta;

	sys->DF2[1][4][0] = 0.0;
	sys->DF2[1][4][1] = 0.0;
	sys->DF2[1][4][2] = 0.0;
	sys->DF2[1][4][3] = 0.0;
	sys->DF2[1][4][4] = 0.0;
	sys->DF2[1][4][5] = 0.0;
	sys->DF2[1][4][6] = 0.0;
	sys->DF2[1][4][7] = 0.0;

	sys->DF2[1][5][0] = 0.0;
	sys->DF2[1][5][1] = 0.0;
	sys->DF2[1][5][2] = 0.0;
	sys->DF2[1][5][3] = 0.0;
	sys->DF2[1][5][4] = 0.0;
	sys->DF2[1][5][5] = 0.0;
	sys->DF2[1][5][6] = 0.0;
	sys->DF2[1][5][7] = 0.0;

	sys->DF2[1][6][0] = 0.0;
	sys->DF2[1][6][1] = 0.0;
	sys->DF2[1][6][2] = 0.0;
	sys->DF2[1][6][3] = 0.0;
	sys->DF2[1][6][4] = 0.0;
	sys->DF2[1][6][5] = 0.0;
	sys->DF2[1][6][6] = 0.0;
	sys->DF2[1][6][7] = 0.0;

	sys->DF2[1][7][0] = -qdd[1] * c * delta;
	sys->DF2[1][7][1] =  qdd[1] * d * delta;
	sys->DF2[1][7][2] = 0.0;
	sys->DF2[1][7][3] = qdd[1] * delta * delta;
	sys->DF2[1][7][4] = 0.0;
	sys->DF2[1][7][5] = 0.0;
	sys->DF2[1][7][6] = 0.0;
	sys->DF2[1][7][7] = qdd[1] * delta * delta;

	sys->DF2[2][0][0] =  qdd[2] * delta * delta;
	sys->DF2[2][0][1] =  0.0;
	sys->DF2[2][0][2] =  qdd[2] * a * delta; 
	sys->DF2[2][0][3] = -qdd[2] * b * delta;
	sys->DF2[2][0][4] =  qdd[2] * delta * delta;
	sys->DF2[2][0][5] = 0.0;
	sys->DF2[2][0][6] = 0.0;
	sys->DF2[2][0][7] = 0.0;

	sys->DF2[2][1][0] = 0.0;
	sys->DF2[2][1][1] = 0.0;
	sys->DF2[2][1][2] = 0.0;
	sys->DF2[2][1][3] = 0.0;
	sys->DF2[2][1][4] = 0.0;
	sys->DF2[2][1][5] = 0.0;
	sys->DF2[2][1][6] = 0.0;
	sys->DF2[2][1][7] = 0.0;

	sys->DF2[2][2][0] =  qdd[2] * a * delta;
	sys->DF2[2][2][1] =  0.0;
	sys->DF2[2][2][2] =  qdd[2] * a * a; 
	sys->DF2[2][2][3] = -qdd[2] * a * b;
	sys->DF2[2][2][4] =  qdd[2] * a * delta;
	sys->DF2[2][2][5] = 0.0;
	sys->DF2[2][2][6] = 0.0;
	sys->DF2[2][2][7] = 0.0;

	sys->DF2[2][3][0] = -qdd[2] * b * delta;
	sys->DF2[2][3][1] =  0.0;
	sys->DF2[2][3][2] = -qdd[2] * a * b; 
	sys->DF2[2][3][3] =  qdd[2] * b * b;
	sys->DF2[2][3][4] = -qdd[2] * b * delta;
	sys->DF2[2][3][5] = 0.0;
	sys->DF2[2][3][6] = 0.0;
	sys->DF2[2][3][7] = 0.0;

	sys->DF2[2][4][0] =  qdd[2] * delta * delta;
	sys->DF2[2][4][1] =  0.0;
	sys->DF2[2][4][2] =  qdd[2] * a * delta; 
	sys->DF2[2][4][3] = -qdd[2] * b * delta;
	sys->DF2[2][4][4] =  qdd[2] * delta * delta;
	sys->DF2[2][4][5] = 0.0;
	sys->DF2[2][4][6] = 0.0;
	sys->DF2[2][4][7] = 0.0;

	sys->DF2[2][5][0] = 0.0;
	sys->DF2[2][5][1] = 0.0;
	sys->DF2[2][5][2] = 0.0;
	sys->DF2[2][5][3] = 0.0;
	sys->DF2[2][5][4] = 0.0;
	sys->DF2[2][5][5] = 0.0;
	sys->DF2[2][5][6] = 0.0;
	sys->DF2[2][5][7] = 0.0;

	sys->DF2[2][6][0] = 0.0;
	sys->DF2[2][6][1] = 0.0;
	sys->DF2[2][6][2] = 0.0;
	sys->DF2[2][6][3] = 0.0;
	sys->DF2[2][6][4] = 0.0;
	sys->DF2[2][6][5] = 0.0;
	sys->DF2[2][6][6] = 0.0;
	sys->DF2[2][6][7] = 0.0;

	sys->DF2[2][7][0] = 0.0;
	sys->DF2[2][7][1] = 0.0;
	sys->DF2[2][7][2] = 0.0;
	sys->DF2[2][7][3] = 0.0;
	sys->DF2[2][7][4] = 0.0;
	sys->DF2[2][7][5] = 0.0;
	sys->DF2[2][7][6] = 0.0;
	sys->DF2[2][7][7] = 0.0;

	sys->DF2[3][0][0] = 0.0;
	sys->DF2[3][0][1] = 0.0;
	sys->DF2[3][0][2] = 0.0;
	sys->DF2[3][0][3] = 0.0;
	sys->DF2[3][0][4] = 0.0;
	sys->DF2[3][0][5] = 0.0;
	sys->DF2[3][0][6] = 0.0;
	sys->DF2[3][0][7] = 0.0;

	sys->DF2[3][1][0] = 0.0;
	sys->DF2[3][1][1] =  qdd[3] * delta * delta;
	sys->DF2[3][1][2] = -qdd[3] * c * delta;
	sys->DF2[3][1][3] =  qdd[3] * d * delta;
	sys->DF2[3][1][4] = 0.0;
	sys->DF2[3][1][5] = qdd[3] * delta * delta;
	sys->DF2[3][1][6] = 0.0;
	sys->DF2[3][1][7] = 0.0;

	sys->DF2[3][2][0] = 0.0;
	sys->DF2[3][2][1] = -qdd[3] * c * delta;
	sys->DF2[3][2][2] =  qdd[3] * c * c;
	sys->DF2[3][2][3] = -qdd[3] * c * d;
	sys->DF2[3][2][4] = 0.0;
	sys->DF2[3][2][5] = -qdd[3] * c * delta;
	sys->DF2[3][2][6] = 0.0;
	sys->DF2[3][2][7] = 0.0;

	sys->DF2[3][3][0] = 0.0;
	sys->DF2[3][3][1] =  qdd[3] * d * delta;
	sys->DF2[3][3][2] = -qdd[3] * c * d;
	sys->DF2[3][3][3] =  qdd[3] * d * d;
	sys->DF2[3][3][4] = 0.0;
	sys->DF2[3][3][5] =  qdd[3] * d * delta;
	sys->DF2[3][3][6] = 0.0;
	sys->DF2[3][3][7] = 0.0;

	sys->DF2[3][4][0] = 0.0;
	sys->DF2[3][4][1] = 0.0;
	sys->DF2[3][4][2] = 0.0;
	sys->DF2[3][4][3] = 0.0;
	sys->DF2[3][4][4] = 0.0;
	sys->DF2[3][4][5] = 0.0;
	sys->DF2[3][4][6] = 0.0;
	sys->DF2[3][4][7] = 0.0;

	sys->DF2[3][5][0] = 0.0;
	sys->DF2[3][5][1] =  qdd[3] * delta * delta;
	sys->DF2[3][5][2] = -qdd[3] * c * delta;
	sys->DF2[3][5][3] =  qdd[3] * d * delta;
	sys->DF2[3][5][4] = 0.0;
	sys->DF2[3][5][5] = qdd[3] * delta * delta;
	sys->DF2[3][5][6] = 0.0;
	sys->DF2[3][5][7] = 0.0;

	sys->DF2[3][6][0] = 0.0;
	sys->DF2[3][6][1] = 0.0;
	sys->DF2[3][6][2] = 0.0;
	sys->DF2[3][6][3] = 0.0;
	sys->DF2[3][6][4] = 0.0;
	sys->DF2[3][6][5] = 0.0;
	sys->DF2[3][6][6] = 0.0;
	sys->DF2[3][6][7] = 0.0;

	sys->DF2[3][7][0] = 0.0;
	sys->DF2[3][7][1] = 0.0;
	sys->DF2[3][7][2] = 0.0;
	sys->DF2[3][7][3] = 0.0;
	sys->DF2[3][7][4] = 0.0;
	sys->DF2[3][7][5] = 0.0;
	sys->DF2[3][7][6] = 0.0;
	sys->DF2[3][7][7] = 0.0;

	sys->DF2[4][0][0] = 0.0;
	sys->DF2[4][0][1] = 0.0;
	sys->DF2[4][0][2] = 0.0;
	sys->DF2[4][0][3] = 0.0;
	sys->DF2[4][0][4] = 0.0;
	sys->DF2[4][0][5] = 0.0;
	sys->DF2[4][0][6] = 0.0;
	sys->DF2[4][0][7] = 0.0;

	sys->DF2[4][1][0] = 0.0;
	sys->DF2[4][1][1] = 0.0;
	sys->DF2[4][1][2] = 0.0;
	sys->DF2[4][1][3] = 0.0;
	sys->DF2[4][1][4] = 0.0;
	sys->DF2[4][1][5] = 0.0;
	sys->DF2[4][1][6] = 0.0;
	sys->DF2[4][1][7] = 0.0;

	sys->DF2[4][2][0] = 0.0;
	sys->DF2[4][2][1] = 0.0;
	sys->DF2[4][2][2] =  qdd[4] * delta * delta;
	sys->DF2[4][2][3] = 0.0;
	sys->DF2[4][2][4] =  qdd[4] * a * delta;
	sys->DF2[4][2][5] = -qdd[4] * b * delta;
	sys->DF2[4][2][6] =  qdd[4] * delta * delta;
	sys->DF2[4][2][7] = 0.0; 

	sys->DF2[4][3][0] = 0.0;
	sys->DF2[4][3][1] = 0.0;
	sys->DF2[4][3][2] = 0.0;
	sys->DF2[4][3][3] = 0.0;
	sys->DF2[4][3][4] = 0.0;
	sys->DF2[4][3][5] = 0.0;
	sys->DF2[4][3][6] = 0.0;
	sys->DF2[4][3][7] = 0.0;

	sys->DF2[4][4][0] = 0.0;
	sys->DF2[4][4][1] = 0.0;
	sys->DF2[4][4][2] =  qdd[4] * a * delta;
	sys->DF2[4][4][3] = 0.0;
	sys->DF2[4][4][4] =  qdd[4] * a * a;
	sys->DF2[4][4][5] = -qdd[4] * a * b;
	sys->DF2[4][4][6] =  qdd[4] * a * delta;
	sys->DF2[4][4][7] = 0.0; 

	sys->DF2[4][5][0] = 0.0;
	sys->DF2[4][5][1] = 0.0;
	sys->DF2[4][5][2] = -qdd[4] * b * delta;
	sys->DF2[4][5][3] = 0.0;
	sys->DF2[4][5][4] = -qdd[4] * a * b;
	sys->DF2[4][5][5] =  qdd[4] * b * b;
	sys->DF2[4][5][6] = -qdd[4] * b * delta;
	sys->DF2[4][5][7] = 0.0; 

	sys->DF2[4][6][0] = 0.0;
	sys->DF2[4][6][1] = 0.0;
	sys->DF2[4][6][2] =  qdd[4] * delta * delta;
	sys->DF2[4][6][3] = 0.0;
	sys->DF2[4][6][4] =  qdd[4] * a * delta;
	sys->DF2[4][6][5] = -qdd[4] * b * delta;
	sys->DF2[4][6][6] =  qdd[4] * delta * delta;
	sys->DF2[4][6][7] = 0.0; 

	sys->DF2[4][7][0] = 0.0;
	sys->DF2[4][7][1] = 0.0;
	sys->DF2[4][7][2] = 0.0;
	sys->DF2[4][7][3] = 0.0;
	sys->DF2[4][7][4] = 0.0;
	sys->DF2[4][7][5] = 0.0;
	sys->DF2[4][7][6] = 0.0;
	sys->DF2[4][7][7] = 0.0;

	sys->DF2[5][0][0] = 0.0;
	sys->DF2[5][0][1] = 0.0;
	sys->DF2[5][0][2] = 0.0;
	sys->DF2[5][0][3] = 0.0;
	sys->DF2[5][0][4] = 0.0;
	sys->DF2[5][0][5] = 0.0;
	sys->DF2[5][0][6] = 0.0;
	sys->DF2[5][0][7] = 0.0;

	sys->DF2[5][1][0] = 0.0;
	sys->DF2[5][1][1] = 0.0;
	sys->DF2[5][1][2] = 0.0;
	sys->DF2[5][1][3] = 0.0;
	sys->DF2[5][1][4] = 0.0;
	sys->DF2[5][1][5] = 0.0;
	sys->DF2[5][1][6] = 0.0;
	sys->DF2[5][1][7] = 0.0;

	sys->DF2[5][2][0] = 0.0;
	sys->DF2[5][2][1] = 0.0;
	sys->DF2[5][2][2] = 0.0;
	sys->DF2[5][2][3] = 0.0;
	sys->DF2[5][2][4] = 0.0;
	sys->DF2[5][2][5] = 0.0;
	sys->DF2[5][2][6] = 0.0;
	sys->DF2[5][2][7] = 0.0;

	sys->DF2[5][3][0] = 0.0;
	sys->DF2[5][3][1] = 0.0;
	sys->DF2[5][3][2] = 0.0;
	sys->DF2[5][3][3] =  qdd[5] * delta * delta;
	sys->DF2[5][3][4] = -qdd[5] * c * delta;
	sys->DF2[5][3][5] =  qdd[5] * d * delta;
	sys->DF2[5][3][6] = 0.0;
	sys->DF2[5][3][7] =  qdd[5] * delta * delta;

	sys->DF2[5][4][0] = 0.0;
	sys->DF2[5][4][1] = 0.0;
	sys->DF2[5][4][2] = 0.0;
	sys->DF2[5][4][3] = -qdd[5] * c * delta;
	sys->DF2[5][4][4] =  qdd[5] * c * c;
	sys->DF2[5][4][5] = -qdd[5] * c * d;
	sys->DF2[5][4][6] = 0.0;
	sys->DF2[5][4][7] = -qdd[5] * c * delta;

	sys->DF2[5][5][0] = 0.0;
	sys->DF2[5][5][1] = 0.0;
	sys->DF2[5][5][2] = 0.0;
	sys->DF2[5][5][3] =  qdd[5] * d * delta;
	sys->DF2[5][5][4] = -qdd[5] * c * d;
	sys->DF2[5][5][5] =  qdd[5] * d * d;
	sys->DF2[5][5][6] = 0.0;
	sys->DF2[5][5][7] =  qdd[5] * d * delta;

	sys->DF2[5][6][0] = 0.0;
	sys->DF2[5][6][1] = 0.0;
	sys->DF2[5][6][2] = 0.0;
	sys->DF2[5][6][3] = 0.0;
	sys->DF2[5][6][4] = 0.0;
	sys->DF2[5][6][5] = 0.0;
	sys->DF2[5][6][6] = 0.0;
	sys->DF2[5][6][7] = 0.0;

	sys->DF2[5][7][0] = 0.0;
	sys->DF2[5][7][1] = 0.0;
	sys->DF2[5][7][2] = 0.0;
	sys->DF2[5][7][3] =  qdd[5] * delta * delta;
	sys->DF2[5][7][4] = -qdd[5] * c * delta;
	sys->DF2[5][7][5] =  qdd[5] * d * delta;
	sys->DF2[5][7][6] = 0.0;
	sys->DF2[5][7][7] =  qdd[5] * delta * delta;

	sys->DF2[6][0][0] =  qdd[6] * delta * delta;
	sys->DF2[6][0][1] = 0.0;
	sys->DF2[6][0][2] = 0.0;
	sys->DF2[6][0][3] = 0.0;
	sys->DF2[6][0][4] =  qdd[6] * delta * delta;
	sys->DF2[6][0][5] = 0.0;
	sys->DF2[6][0][6] =  qdd[6] * a * delta;
	sys->DF2[6][0][7] = -qdd[6] * b * delta;

	sys->DF2[6][1][0] = 0.0;
	sys->DF2[6][1][1] = 0.0;
	sys->DF2[6][1][2] = 0.0;
	sys->DF2[6][1][3] = 0.0;
	sys->DF2[6][1][4] = 0.0;
	sys->DF2[6][1][5] = 0.0;
	sys->DF2[6][1][6] = 0.0;
	sys->DF2[6][1][7] = 0.0;

	sys->DF2[6][2][0] = 0.0;
	sys->DF2[6][2][1] = 0.0;
	sys->DF2[6][2][2] = 0.0;
	sys->DF2[6][2][3] = 0.0;
	sys->DF2[6][2][4] = 0.0;
	sys->DF2[6][2][5] = 0.0;
	sys->DF2[6][2][6] = 0.0;
	sys->DF2[6][2][7] = 0.0;

	sys->DF2[6][3][0] = 0.0;
	sys->DF2[6][3][1] = 0.0;
	sys->DF2[6][3][2] = 0.0;
	sys->DF2[6][3][3] = 0.0;
	sys->DF2[6][3][4] = 0.0;
	sys->DF2[6][3][5] = 0.0;
	sys->DF2[6][3][6] = 0.0;
	sys->DF2[6][3][7] = 0.0;

	sys->DF2[6][4][0] =  qdd[6] * delta * delta;
	sys->DF2[6][4][1] = 0.0;
	sys->DF2[6][4][2] = 0.0;
	sys->DF2[6][4][3] = 0.0;
	sys->DF2[6][4][4] =  qdd[6] * delta * delta;
	sys->DF2[6][4][5] = 0.0;
	sys->DF2[6][4][6] =  qdd[6] * a * delta;
	sys->DF2[6][4][7] = -qdd[6] * b * delta;

	sys->DF2[6][5][0] = 0.0;
	sys->DF2[6][5][1] = 0.0;
	sys->DF2[6][5][2] = 0.0;
	sys->DF2[6][5][3] = 0.0;
	sys->DF2[6][5][4] = 0.0;
	sys->DF2[6][5][5] = 0.0;
	sys->DF2[6][5][6] = 0.0;
	sys->DF2[6][5][7] = 0.0;

	sys->DF2[6][6][0] =  qdd[6] * a * delta;
	sys->DF2[6][6][1] = 0.0;
	sys->DF2[6][6][2] = 0.0;
	sys->DF2[6][6][3] = 0.0;
	sys->DF2[6][6][4] =  qdd[6] * a * delta;
	sys->DF2[6][6][5] = 0.0;
	sys->DF2[6][6][6] =  qdd[6] * a * a;
	sys->DF2[6][6][7] = -qdd[6] * a * b;

	sys->DF2[6][7][0] = -qdd[6] * b * delta;
	sys->DF2[6][7][1] = 0.0;
	sys->DF2[6][7][2] = 0.0;
	sys->DF2[6][7][3] = 0.0;
	sys->DF2[6][7][4] = -qdd[6] * b * delta;
	sys->DF2[6][7][5] = 0.0;
	sys->DF2[6][7][6] = -qdd[6] * a * b;
	sys->DF2[6][7][7] =  qdd[6] * b * b;

	sys->DF2[7][0][0] = 0.0;
	sys->DF2[7][0][1] = 0.0;
	sys->DF2[7][0][2] = 0.0;
	sys->DF2[7][0][3] = 0.0;
	sys->DF2[7][0][4] = 0.0;
	sys->DF2[7][0][5] = 0.0;
	sys->DF2[7][0][6] = 0.0;
	sys->DF2[7][0][7] = 0.0;

	sys->DF2[7][1][0] = 0.0;
	sys->DF2[7][1][1] =  qdd[7] * delta * delta;
	sys->DF2[7][1][2] = 0.0;
	sys->DF2[7][1][3] = 0.0;
	sys->DF2[7][1][4] = 0.0;
	sys->DF2[7][1][5] =  qdd[7] * delta * delta;
	sys->DF2[7][1][6] = -qdd[7] * c * delta;
	sys->DF2[7][1][7] =  qdd[7] * d * delta;

	sys->DF2[7][2][0] = 0.0;
	sys->DF2[7][2][1] = 0.0;
	sys->DF2[7][2][2] = 0.0;
	sys->DF2[7][2][3] = 0.0;
	sys->DF2[7][2][4] = 0.0;
	sys->DF2[7][2][5] = 0.0;
	sys->DF2[7][2][6] = 0.0;
	sys->DF2[7][2][7] = 0.0;

	sys->DF2[7][3][0] = 0.0;
	sys->DF2[7][3][1] = 0.0;
	sys->DF2[7][3][2] = 0.0;
	sys->DF2[7][3][3] = 0.0;
	sys->DF2[7][3][4] = 0.0;
	sys->DF2[7][3][5] = 0.0;
	sys->DF2[7][3][6] = 0.0;
	sys->DF2[7][3][7] = 0.0;

	sys->DF2[7][4][0] = 0.0;
	sys->DF2[7][4][1] = 0.0;
	sys->DF2[7][4][2] = 0.0;
	sys->DF2[7][4][3] = 0.0;
	sys->DF2[7][4][4] = 0.0;
	sys->DF2[7][4][5] = 0.0;
	sys->DF2[7][4][6] = 0.0;
	sys->DF2[7][4][7] = 0.0;

	sys->DF2[7][5][0] = 0.0;
	sys->DF2[7][5][1] =  qdd[7] * delta * delta;
	sys->DF2[7][5][2] = 0.0;
	sys->DF2[7][5][3] = 0.0;
	sys->DF2[7][5][4] = 0.0;
	sys->DF2[7][5][5] =  qdd[7] * delta * delta;
	sys->DF2[7][5][6] = -qdd[7] * c * delta;
	sys->DF2[7][5][7] =  qdd[7] * d * delta;

	sys->DF2[7][6][0] = 0.0;
	sys->DF2[7][6][1] = -qdd[7] * c * delta;
	sys->DF2[7][6][2] = 0.0;
	sys->DF2[7][6][3] = 0.0;
	sys->DF2[7][6][4] = 0.0;
	sys->DF2[7][6][5] = -qdd[7] * c * delta;
	sys->DF2[7][6][6] =  qdd[7] * c * c;
	sys->DF2[7][6][7] = -qdd[7] * c * d;

	sys->DF2[7][7][0] = 0.0;
	sys->DF2[7][7][1] =  qdd[7] * d * delta;
	sys->DF2[7][7][2] = 0.0;
	sys->DF2[7][7][3] = 0.0;
	sys->DF2[7][7][4] = 0.0;
	sys->DF2[7][7][5] =  qdd[7] * d * delta;
	sys->DF2[7][7][6] = -qdd[7] * c * d;
	sys->DF2[7][7][7] =  qdd[7] * d * d;


	switch(sys->variable){
		case DELTA:
		sys->Dl[0] =  qd[0] * (x[6] + x[2]);
		sys->Dl[1] = -qd[1] * (x[7] + x[3]);
		sys->Dl[2] =  qd[2] * (x[0] + x[4]);
		sys->Dl[3] = -qd[3] * (x[1] + x[5]);
		sys->Dl[4] =  qd[4] * (x[2] + x[6]);
		sys->Dl[5] = -qd[5] * (x[3] + x[7]);
		sys->Dl[6] =  qd[6] * (x[4] + x[0]);
		sys->Dl[7] = -qd[7] * (x[5] + x[1]);

		sys->DFL[0][0] =  qdd[0] * (x[6] + x[2]) * a;
		sys->DFL[0][1] = -qdd[0] * (x[6] + x[2]) * b;
		sys->DFL[0][2] =  qdd[0] * (x[6] + x[2]) * delta + q[0];
		sys->DFL[0][3] = 0.0;
		sys->DFL[0][4] = 0.0;
		sys->DFL[0][5] = 0.0;
		sys->DFL[0][6] =  qdd[0] * (x[6] + x[2]) * delta + q[0];
		sys->DFL[0][7] = 0.0;

		sys->DFL[1][0] = -qdd[0] * (x[7] + x[3]) * c;
		sys->DFL[1][1] =  qdd[0] * (x[7] + x[3]) * d;
		sys->DFL[1][2] =  0.0;
		sys->DFL[1][3] =  qdd[0] * (x[7] + x[3]) * delta - q[1];
		sys->DFL[1][4] = 0.0;
		sys->DFL[1][5] = 0.0;
		sys->DFL[1][6] = 0.0;
		sys->DFL[1][7] =  qdd[0] * (x[7] + x[3]) * delta - q[1];

		sys->DFL[2][0] =  qdd[2] * (x[0] + x[4]) * delta + q[2];
		sys->DFL[2][1] = 0.0;
		sys->DFL[2][2] =  qdd[2] * (x[0] + x[4]) * a;
		sys->DFL[2][3] = -qdd[2] * (x[0] + x[4]) * b;
		sys->DFL[2][4] =  qdd[2] * (x[0] + x[4]) * delta + q[2];
		sys->DFL[2][5] = 0.0;
		sys->DFL[2][6] = 0.0;
		sys->DFL[2][7] = 0.0;

		sys->DFL[3][0] = 0.0;
		sys->DFL[3][1] =  qdd[3] * (x[1] + x[5]) * delta - q[3];
		sys->DFL[3][2] = -qdd[3] * (x[1] + x[5]) * c;
		sys->DFL[3][3] =  qdd[3] * (x[1] + x[5]) * d;
		sys->DFL[3][4] = 0.0;
		sys->DFL[3][5] =  qdd[3] * (x[1] + x[5]) * delta - q[3];
		sys->DFL[3][6] = 0.0;
		sys->DFL[3][7] = 0.0;

		sys->DFL[4][0] = 0.0;
		sys->DFL[4][1] = 0.0;
		sys->DFL[4][2] =  qdd[4] * (x[2] + x[6]) * delta + q[4];
		sys->DFL[4][3] = 0.0;
		sys->DFL[4][4] =  qdd[4] * (x[2] + x[6]) * a;
		sys->DFL[4][5] = -qdd[4] * (x[2] + x[6]) * b;
		sys->DFL[4][6] =  qdd[4] * (x[2] + x[6]) * delta + q[4];
		sys->DFL[4][7] = 0.0;
		
		sys->DFL[5][0] = 0.0;
		sys->DFL[5][1] = 0.0;
		sys->DFL[5][2] = 0.0;
		sys->DFL[5][3] =  qdd[5] * (x[3] + x[7]) * delta - q[5];
		sys->DFL[5][4] = -qdd[5] * (x[3] + x[7]) * c;
		sys->DFL[5][5] =  qdd[5] * (x[3] + x[7]) * d;
		sys->DFL[5][6] = 0.0;
		sys->DFL[5][7] =  qdd[5] * (x[3] + x[7]) * delta - q[5];

		sys->DFL[6][0] =  qdd[6] * (x[4] + x[0]) * delta + q[6];
		sys->DFL[6][1] = 0.0;
		sys->DFL[6][2] = 0.0;
		sys->DFL[6][3] = 0.0;
		sys->DFL[6][4] =  qdd[6] * (x[4] + x[0]) * delta + q[6];
		sys->DFL[6][5] = 0.0;
		sys->DFL[6][6] =  qdd[6] * (x[4] + x[0]) * a;
		sys->DFL[6][7] = -qdd[6] * (x[4] + x[0]) * b;
		
		sys->DFL[7][0] = 0.0;
		sys->DFL[7][1] =  qdd[7] * (x[5] + x[1]) * delta - q[7]; 
		sys->DFL[7][2] = 0.0;
		sys->DFL[7][3] = 0.0;
		sys->DFL[7][4] = 0.0;
		sys->DFL[7][5] =  qdd[7] * (x[5] + x[1]) * delta - q[7]; 
		sys->DFL[7][6] = -qdd[7] * (x[5] + x[1]) * c;
		sys->DFL[7][7] =  qdd[7] * (x[5] + x[1]) * d;

		break;
	}

	/* complete differential eqs. */

	for (i = 0; i < NDE; i++){ 
		f[i] = -x[i] + q[i];
	}



	/* ############################################################ */

	/* DO NOT EDIT BELOW */

	for (i = 0; i < NDE; i++){
		for (j = 0; j < NDE; j++){
			for (k = j + 1; k < NDE; k++){
				double diff;
				diff = fabs(sys->DF2[i][j][k] - sys->DF2[i][k][j]);
				if (diff > 1e-8 && sys->mode != MODE_F){
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

    for (j = 0; j < NDE; j++){ f[NDE + NDE * NDE + j] += sys->Dl[j]; }

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
                Dfx2[i][j][k] = inner(sys->DF2[j][k], buf, NDE);
            }
        }
    }

    for (i = 0; i < NDE; i++){
        for (j = 0; j < NDE; j++){
            buf[j] = x[NDE + NDE * i + j]; 
        }
        for (j = 0; j < NDE; j++){
            Dfl[i][j] = inner(sys->DFL[j], buf, NDE);
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

    /* comprete 2nd variations for parameters (1)*/

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
}

	/* ############################################################ */



void StoreDf(double x[], double q[], SysData *sys)
{
	int i;
	double a, b, c, d, rx, ry, delta;
	double sigmoid(double);
	double qd[NDE];

	a = sys->param[A];
	b = sys->param[B];
	c = sys->param[C];
	d = sys->param[D];
	rx = sys->param[RX];
	ry = sys->param[RY];
	delta = sys->param[DELTA];

    q[0] = sigmoid(a * x[0] - b * x[1] + rx + delta * (x[6] + x[2]));
    q[1] = sigmoid(c * x[0] - d * x[1] + ry - delta * (x[7] + x[3]));
    q[2] = sigmoid(a * x[2] - b * x[3] + rx + delta * (x[0] + x[4]));
    q[3] = sigmoid(c * x[2] - d * x[3] + ry - delta * (x[1] + x[5]));
    q[4] = sigmoid(a * x[4] - b * x[5] + rx + delta * (x[2] + x[6]));
    q[5] = sigmoid(c * x[4] - d * x[5] + ry - delta * (x[3] + x[7]));
    q[6] = sigmoid(a * x[6] - b * x[7] + rx + delta * (x[4] + x[0]));
    q[7] = sigmoid(c * x[6] - d * x[7] + ry - delta * (x[5] + x[1]));


	for (i = 0; i < NDE; i++){ qd[i] = (1.0 - q[i]) * q[i]; }

    sys->Df[0][0] = -1.0 + qd[0] * a;
    sys->Df[0][1] = -qd[0] * b;
    sys->Df[0][2] = qd[0] * delta;
    sys->Df[0][3] = 0.0;
    sys->Df[0][4] = 0.0;
    sys->Df[0][5] = 0.0;
    sys->Df[0][6] = qd[0] * delta;
    sys->Df[0][7] = 0.0;
    
    sys->Df[1][0] = qd[1]* c;
    sys->Df[1][1] = -1.0 - qd[1]* d;
    sys->Df[1][2] = 0.0;
    sys->Df[1][3] = -qd[1] * delta;
    sys->Df[1][4] = 0.0;
    sys->Df[1][5] = 0.0;
    sys->Df[1][6] = 0.0;
    sys->Df[1][7] = -qd[1] * delta;

    sys->Df[2][0] = qd[2] * delta;
    sys->Df[2][1] = 0.0;
    sys->Df[2][2] = -1 + qd[2] * a;
    sys->Df[2][3] = -qd[2] * b;
    sys->Df[2][4] =  qd[2] * delta;
    sys->Df[2][5] =  0.0;
    sys->Df[2][6] =  0.0;
    sys->Df[2][7] =  0.0;

    sys->Df[3][0] =  0.0;
    sys->Df[3][1] = -qd[3] * delta;
    sys->Df[3][2] =  qd[3] * c;
    sys->Df[3][3] =  -1.0 - qd[3] * d;
    sys->Df[3][4] =  0.0;
    sys->Df[3][5] =  -qd[3] * delta;
    sys->Df[3][6] =  0.0;
    sys->Df[3][7] =  0.0;

    sys->Df[4][0] =  0.0;
    sys->Df[4][1] =  0.0;
    sys->Df[4][2] = qd[4] * delta;
    sys->Df[4][3] =  0.0;
    sys->Df[4][4] = -1.0+ qd[4]* a;
    sys->Df[4][5] =  -qd[4] * b;
    sys->Df[4][6] =  qd[4] * delta;
    sys->Df[4][7] =  0.0;

    sys->Df[5][0] =  0.0;
    sys->Df[5][1] =  0.0;
    sys->Df[5][2] =  0.0;
    sys->Df[5][3] =  -qd[5] * delta;
    sys->Df[5][4] =  qd[5] * c;
    sys->Df[5][5] =  -1.0 - qd[5] * d;
    sys->Df[5][6] =  0.0;
    sys->Df[5][7] =  -qd[5] * delta;

    sys->Df[6][0] =  qd[6] * delta;
    sys->Df[6][1] =  0.0;
    sys->Df[6][2] =  0.0;
    sys->Df[6][3] =  0.0;
    sys->Df[6][4] =  qd[6] * delta;
    sys->Df[6][5] =  0.0;
    sys->Df[6][6] = -1.0 + qd[6] * a;
    sys->Df[6][7] = -qd[6] * b;

    sys->Df[7][0] =  0.0;
    sys->Df[7][1] = -qd[7] * delta;
    sys->Df[7][2] =  0.0;
    sys->Df[7][3] =  0.0;
    sys->Df[7][4] =  0.0;
    sys->Df[7][5] = -qd[7] * delta;
    sys->Df[7][6] =  qd[7] * c;
    sys->Df[7][7] = -1.0 - qd[7] * d;

}

double sigmoid(double x)
{
    return 1.0 / (1.0 + exp(-x));
}

