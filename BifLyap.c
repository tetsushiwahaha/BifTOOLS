#include "BifSysDepend.h"

static char rcsid[] = "$Id: BifLyap.c,v 1.1 2006/08/25 03:39:08 tetsushi Exp $";

int lyapunov(Sysdata *sys)
{

	double *x, h;
	double xstart[NDE];
	double xbak[NDVE1];

	x = sys->x;
	h = sys->h;

	for (i = 0; i < NDE; i++) xstart[i] = xbak[i] = x[i];

	while (1){

		for (i = 0; i < NDE; i++){
			for (j = 0; j < NDE; j++){
				if (i == j){ x[NDE + NDE  * j + i] = 1.0; }
				else { x[NDE + NDE  * j + i] = 0.0;}
			}
		}
		for (loop = 0; loop < LE_IDLING; loop++){
			while(1){
				for (i = 0; i < NDVE1; i++) xbak[i] = x[i];
				runge(NDVE1, h, x, t, sys);
				sign = x[sys->section] - xstart[sys->section];

					? 1.0 : -1.0;

				if (x[sys->section] - xstart[sys->section] > 0.0


			
		}

		for (loop = 0; loop < LE_ITERATION; loop++){
		}




	
}


void sysvar(SysData *sys)
{
	int i,j;
	int loop0, loop1;
	double t, h;
	double tau = 0.0;
	double *x;
	double init[NDE][NDE];
	double dummy[NVALL];

	x = sys->x;

	for (i = NDE; i < NVALL; i++){ x[i] = 0.0; }

	for (i = 0; i < NDE; i++){
		for (j = 0; j < NDE; j++){
			if (i == j) x[NDE + j * NDE + i] = 1.0;
			else x[NDE + j * NDE + i] = 0.0; 
		}
	}
	/* to set P matrix */
	function(sys->x, dummy, 0.0, sys); 		

	if (sys->mode == MODE_P){
		/* pitchfork bif. */
		double p[NDE][NDE];
		double z[NDE], res[NDE];

		inverse(sys->P, &p[0][0], NDE);

		for (i = 0; i < NDE; i++){ z[i] = x[i]; }
		matprovec(p, z, res, NDE);
		for (i = 0; i < NDE; i++){ x[i] = res[i]; }

		for (i = 0; i < NDE; i++){
			for (j = 0; j < NDE; j++){
				x[NDE + j * NDE + i] = p[i][j];
			}
		}
	}


	h = sys->tau / sys->m;

	for (loop0 = 0; loop0 < sys->l; loop0++){
		for (loop1 = 0; loop1 < sys->m; loop1++){
			t = loop1 * h;
			runge(NVALL, h, sys->x, t, sys);
		}
		for (j = 0; j < NDE; j++) sys->xn[loop0][j] = sys->x[j];
	}

#ifdef DEBUG_SYSVAR
	printf("tau = %lf \n", sys->tau);
	for (i = 0; i < NDE; i++){
		printf("x[%d] = %lf\n", i, x[i]);
	}
	for (i = 0; i < NDE; i++){
		for (j = 0; j < NDE; j++){
			printf("%lf ", x[NDE + j * NDE + i] );
		}
	}
#endif
}

