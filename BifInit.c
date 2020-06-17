#include "BifSysDepend.h"
#include <unistd.h>
#include <string.h>

static char rcsid[] = "$Id: BifInit.c,v 1.4 2019/11/23 01:27:22 tetsushi Exp $";

extern int optind;

SysData *BifInit(int argc, char **argv)
{
	char modeflag;
	char cmd[BUFSIZ];
	char command[BUFSIZ];
	char fout[BUFSIZ], fout0[BUFSIZ];
	int i, j, k, ch, n;
	int iteration;
	int outflag = 0;
	double *x;
	double u[NDE];
	double period;
	SysData *sys;
	FILE *fpin;

	sys = (SysData *)calloc(sizeof (SysData), 1);

	while ((ch = getopt(argc, argv, "o")) != -1){
		switch (ch){
			case 'o':
				outflag = 1;
			break;
			default:
			;
		}
	}
	argc -= optind;
	argv += optind;

	sprintf(cmd, "/usr/bin/cpp -E %s | sed 's/^#.*$//g'", *argv);
	if ((fpin = popen(cmd, "r")) == NULL){
		fprintf(stderr, "cannot open %s\n", cmd);
		exit(-1);
	}

	if (outflag){
		sprintf(fout, "%s.out", *argv);
		sprintf(fout0, "%s.log", *argv);

		if((sys->fpout = fopen(fout,"w"))==NULL){
			fprintf(stderr,"cannot open %s.out\n",argv[1]);
			exit(-1);
		}
		if((sys->fpout0 = fopen(fout0,"w"))==NULL){
			fprintf(stderr,"cannot open %s.log\n",argv[1]);
			exit(-1);
		}
		sys->outflag = 1;
	}

	fscanf(fpin, "%s", command);

	sys->mode = -1;

	if (!strcmp(command, "F") || !strcmp(command, "FIX")){
		modeflag = MODE_F;
	} else if (!strcmp(command, "D") || !strcmp(command, "PD") || 
		!strcmp(command, "I") || !strcmp(command, "PERIODDOUBLING")){
		modeflag = MODE_D;
	} else if (!strcmp(command, "T") || !strcmp(command, "G") || 
		!strcmp(command, "TANGENT") ){
		modeflag = MODE_T;
	} else if (!strcmp(command, "N") || !strcmp(command, "NEIMARK")){
		modeflag = MODE_N;
	} else if (!strcmp(command, "P") || !strcmp(command, "PITCHFORK")){
		modeflag = MODE_P;
	} else if (!strcmp(command, "E")){
		modeflag = MODE_E;
	} else if (!strcmp(command, "L")){
		modeflag = MODE_L;
	} else if (!strcmp(command, "H") || !strcmp(command, "EH")){
		modeflag = MODE_H;
	} else if (!strcmp(command, "EG") || !strcmp(command, "ET")){
		modeflag = MODE_ET;
	}

	sys->mode = modeflag;
	if (modeflag == -1){
		fprintf(stderr, "specified mode is not defined\n");
		exit(-1);
	}
	
	for (i = 0; i < NP; i++) { fscanf(fpin,"%lf", &sys->param[i]);}

	for (i = 0; i < NDE; i++){ fscanf(fpin,"%lf",&sys->x[i]); }

	fscanf(fpin,"%lf", &sys->tau);
	fscanf(fpin, "%lf", &sys->sigma);
	fscanf(fpin, "%lf", &sys->omega);

	for (i = 0; i < NP; i++) { fscanf(fpin,"%lf", &sys->dparam[i]);}

	fscanf(fpin,"%d",&sys->l);
	fscanf(fpin,"%d",&sys->section);
	fscanf(fpin,"%d",&sys->m);
	fscanf(fpin,"%lf",&sys->eps0);
	fscanf(fpin,"%lf",&sys->emax);
	fscanf(fpin,"%d",&sys->ite_max);
	fscanf(fpin,"%d",&sys->variable);
	fscanf(fpin,"%d",&sys->increment);

#ifdef DEBUG
	printf("l:%d ",sys->l);
	printf("sec:%d ",sys->section);
	printf("m:%d ",sys->m);
	printf("eps0:%lf ",sys->eps0);
	printf("emax:%lf ",sys->emax);
	printf("itemax:%d ",sys->ite_max);
	printf("val:%d ",sys->variable);
	printf("inc:%d ",sys->increment);
#endif

	pclose(fpin);

	x = sys->x;
	sys->xe = x[sys->section]; 

	/* we assume that omega(=sin(theta)) > 0 */

	sys->omega = fabs(sys->omega);

	if (sys->mode == MODE_N) {	/* assume NS bif. */
		sys->on_real_axis = 0;
		sys->theta = asin(sys->omega);
		if (sys->sigma < 0.0) { 
			sys->theta = M_PI - sys->theta; 
		}
	} else {
		sys->on_real_axis = 1;
	}

#ifdef DEBUG
	printf("sys->l = %d\n", sys->l);
	printf("command = %s\n", command);
#endif

	return sys;
	/* supposed to be not null */

}
