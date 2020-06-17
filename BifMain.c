#include "BifSysDepend.h"
#include <stdlib.h>

static char rcsid[] = "$Id: BifMain.c,v 1.3 2019/11/22 09:18:42 tetsushi Exp $";

extern int optind;

int main(int argc, char **argv)
{
	SysData *sys;
	char message[BUFSIZ];
	void print_message0(char *, int);
	void print_result(SysData *);
	void print_result_lyap(SysData *);

	int i, j, k, iteration;

	sys = BifInit(argc, argv);

	while (1){
		fprintf(stdout, "------------------\n");
	  if (sys->mode == MODE_L){
			strcpy(message, "Lyapunov Exponents");
			print_message0(message, 0);
			/*
			laypunov(sys);
			*/
	  } else {
		iteration = fixed(sys);
		if (iteration > sys->ite_max) { 
			printf("iteration over %d\n", iteration);
			exit(-1);
		}
		switch(sys->mode){
			case MODE_T:
				strcpy(message, "tangent");
				break;
			case MODE_D:
				strcpy(message, "period-doubling");
				break;
			case MODE_F:
				strcpy(message, "tracing fixed point");
				break;
			case MODE_N:
				strcpy(message, "Neimark-Sacker");
				break;
			case MODE_P:
				strcpy(message, "pitchfork");
				break;
			case MODE_E:
				strcpy(message, "tracing equilibrium");
				break;
			case MODE_H:
				strcpy(message, "Hopf bifurcation");
				break;
			case MODE_ET:
				strcpy(message, "tangent bif. of equilibrium");
				break;
			case MODE_EP:
				strcpy(message, "pitchfork of equilibrium");
				break;
		}
		print_message0(message, iteration);

		if (sys->mode == MODE_L){
			print_result_lyap(sys);
		} else {
			print_result(sys);
		}

		printf("param(%d) : %lf  -> ", sys->increment, 
			sys->param[sys->increment]);
		sys->param[sys->increment] += sys->dparam[sys->increment];
		printf("%lf\n", sys->param[sys->increment]);
	}
  }
}


void print_message0(char *message, int iteration)
{
	fprintf(stdout, "\x1b[34m");
	fprintf(stdout, "(%s : %d)\n", message, iteration);
	fprintf(stdout, "\x1b[0m");
	fprintf(stdout, "\x1b[31m");
}




void print_result(SysData *sys){
	int i;
	double zr[NDE], zi[NDE];

	for (i = 0; i < NP; i++){
		if (i == sys->variable){
			fprintf(stdout, "%.14f ", sys->param[i]);
		} else {
			fprintf(stdout, "%.10f ", sys->param[i]);
		}
	}

	fprintf(stdout, "NDE: %d\n", NDE);
	fprintf(stdout, "x: ");
	for (i = 0; i < NDE; i++){ fprintf(stdout, "%.12f ", sys->x[i]); } 
	fprintf(stdout, "%lf\n", sys->tau);

	if (sys->outflag){
		for (i = 0; i < NP; i++){
			fprintf(sys->fpout, "%lf ", sys->param[i]);
		}
		for (i = 0; i < NP; i++){
			if (i == sys->variable){
				fprintf(sys->fpout0, "%.14f ", sys->param[i]);
			} else {
				fprintf(sys->fpout0, "%.12f ", sys->param[i]);
			}
		}
		fprintf(sys->fpout0, "\n");
		for (i = 0; i < NDE; i++){ 
			fprintf(sys->fpout0, "%.14f ", sys->x[i]); 
		} 
		fprintf(sys->fpout0, "\n");
		fprintf(sys->fpout0, "%.14f\n", sys->tau);
		fprintf(sys->fpout, "\n");
	}
	if (sys->mode == MODE_E || sys->mode == MODE_H || sys->mode == MODE_ET){
		DQR(sys->Df, zr, zi, NDE);
	} else {
		DQR(sys->Dx, zr, zi, NDE);
	}
	for (i = 0; i < NDE; i++){
		fprintf(stdout, "(%f%ci%f)\n", zr[i], 
			(zi[i] < 0.0) ? '-' : '+', fabs(zi[i]));
		if (fabs(zi[i]) > EPS){
			fprintf(stdout, "{ %lf }\n", 
				sqrt(zr[i] * zr[i] + zi[i] * zi[i]));
		}
	}
	if (sys->outflag){
		for (i = 0; i < NDE; i++){
			fprintf(sys->fpout0, "(%.12f%ci%.12f)\n", zr[i], 
				(zi[i] < 0.0) ? '-' : '+', fabs(zi[i]));
			if (fabs(zi[i]) > EPS){
				fprintf(sys->fpout0, "{ %lf }\n", 
					sqrt(zr[i] * zr[i] + zi[i] * zi[i]));
			}
		}
		fflush(sys->fpout0);
		fflush(sys->fpout);
	}
	fprintf(stdout, "\x1b[0m");
	fflush(stdout);
}

void print_result_lyap(SysData *a)
{
	;
}
