#include "BifSysDepend.h"

static char rcsid[] = "$Id: BifNewton.c,v 1.8 2019/12/20 03:44:50 tetsushi Exp $";

int fixed(SysData *sys)
{
	int i,j;
	double sum;
	double prev[NNM];
	double dlt0, dlt1, dlt2, dlt3;
	double s;
	double ferror = 0.0;
	void newtonequil(SysData *);
	void newtonequilbif(SysData *);

	i = 0;
	sum = 0.0;

	while (1){

		for (j = 0; j < NDE; j++){ prev[j] = sys->x[j];}


		dlt1 = dlt2 = dlt3 = sum = 0.0;
		switch (sys->mode){
			case MODE_N:
				prev[NDE + 2] = sys->theta;
			case MODE_D:
			case MODE_P:
			case MODE_T:
				prev[NDE + 1] = sys->tau;
			case MODE_H:
			case MODE_EP:
			case MODE_ET:
				prev[NDE] = sys->param[sys->variable];
			break;
		}
		if (sys->mode == MODE_E){ 
			newtonequil(sys); 
		} else if (sys->mode == MODE_H || sys->mode == MODE_ET
			|| sys->mode == MODE_EP){ 
			newtonequilbif(sys); 
		} else { 
			newton(sys); 
		}

		switch (sys->mode){
			case MODE_N:
				dlt3 = fabs(prev[NDE + 2] - sys->theta);
			case MODE_D:
			case MODE_P:
			case MODE_T:
				dlt2 = fabs(prev[NDE + 1] - sys->tau);
			case MODE_H:
			case MODE_ET:
				dlt1 = fabs(prev[NDE] - sys->param[sys->variable]);
			break;
		}

		for ( j = 0; j < NDE; j++){	/* one of them should be zero */
			s = sys->x[j] - prev[j]; 
			s *= s;
			sum += s;
		}
		dlt0 = sqrt(sum) / NDE;

#ifdef SHOWDELTA
		switch (sys->mode){
			case MODE_T:
			case MODE_P:
			case MODE_D:
				printf("delta x, l, t = %.10f %.10f %.10f\n", 
					dlt0, dlt1, dlt2);
			break;
			case MODE_N:
				printf("delta x, l, t, q = %.10f %.10f %.10f %.10f\n", 
					dlt0, dlt1, dlt2, dlt3);
			break;
			case MODE_E:
			case MODE_F:
				printf("delta = %.10f\n", dlt0); 
			break;	
			case MODE_H:
			case MODE_ET:
				printf("delta x, l = %.10f %.10f\n", dlt0, dlt1);
		}	
		fflush(stdout);
#endif
		if (dlt0 > sys->emax || dlt1 > sys->emax){
			fprintf(stderr, "divergence. (%lf) (%lf) (%lf)\n", dlt0, dlt1, sys->emax);
			exit(1);
		}
		i++;

		if (sys->mode == MODE_D || sys->mode == MODE_T || sys->mode == MODE_P){
			if (dlt0 + dlt1 + dlt2  < sys->eps0) break;
		} else if (sys->mode == MODE_N){
			if (dlt0 + dlt1 + dlt3  < sys->eps0) break;
		} else if (sys->mode == MODE_F || sys->mode == MODE_E){
			if (dlt0 < sys->eps0) break;
		} else if (sys->mode == MODE_H || sys->mode == MODE_ET
				|| sys->mode == MODE_EP ){
			if (dlt0 + dlt1 < sys->eps0) break;
		}


		if (i > sys->ite_max){
			printf("iteration over %d\n",sys->ite_max);
			exit(1);
		}
	}
	return i;
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

#ifdef DEBUG_SYSVAR
	printf("BEFOREBEFORE\n");
	printf("tau = %lf \n", sys->tau);
	printf("m = %d \n", sys->m);
	printf("h = %lf \n", h);
	for (i = 0; i < NDE; i++){
		printf("x[%d] = %.15f\n", i, x[i]);
	}
	for (i = 0; i < NDE; i++){
		for (j = 0; j < NDE; j++){
			printf("x(%d) = %lf\n", NDE + i * NDE + j, x[NDE + i * NDE + j] );
		}
	}
	for (i = 0; i < NVALL; i++){
		for (j = 0; j < NDE; j++){
			printf(":x(%d) = %.15f\n", NDE + i * NDE + j, x[NDE + i * NDE +j] );
		}
	}
#endif

	for (loop0 = 0; loop0 < sys->l; loop0++){
		for (loop1 = 0; loop1 < sys->m; loop1++){
			t = loop1 * h;
			runge(NVALL, h, sys->x, t, sys);
		}
		for (j = 0; j < NDE; j++) sys->xn[loop0][j] = sys->x[j];
	}

#ifdef DEBUG_SYSVAR
	printf("AFTERAFTER\n");
	printf("tau = %lf \n", sys->tau);
	printf("m = %d \n", sys->m);
	printf("h = %lf \n", h);
	for (i = 0; i < NDE; i++){
		printf("x[%d] = %.15f\n", i, x[i]);
	}
	for (i = 0; i < NDE; i++){
		for (j = 0; j < NDE; j++){
			printf("x(%d) = %lf\n", NDE + i * NDE + j, x[NDE + i * NDE + j] );
		}
	}
	for (i = 0; i < NVALL; i++){
		for (j = 0; j < NDE; j++){
			printf(":x(%d) = %.15f\n", NDE + i * NDE + j, x[NDE + i * NDE +j] );
		}
	}
#endif
}

void newton(SysData *sys)
{
	int i, j, k, loop;
	int Fpos;
	double t;
	double prev[NNM];
	double h[NNM];

	double Dx[NDE][NDE];
	double NDx[NNM][NNM];
	double Dq[NDE];

	double P[NDE][NDE];
	double Pbak[NDE][NDE];
	double PBAK[NDE][NDE][NDE];
	double Dh[NDE][NDE];
	double dphidl[NDE];
	double Dhinv[NDE][NDE];
	double dphidx[NDE][NDE];
	double buf[NDE][NDE];
	double bufv[NDE];
	double chi;
	double re, im;
	double f0[NDE];
	double f[NVALL];
	DecompMat *dets;

	double F[NDE], F2[NDE];
	double dchidx[NDE], dchidu[NDE];
	double dchidl, dchidt;


	switch (sys->mode){
		case MODE_P:
			/* pitchfork bifurcation */
			Fpos = NDE + 1;
			sys->tau /= 2.0;
			break;
		case MODE_N:
			/* Neimark-Sacker bifurcation */
			prev[NDE + 1] = sys->theta;
			sys->omega = sin(sys->theta);		
			sys->sigma = cos(sys->theta);
			Fpos = NDE + 2;
			break;
		case MODE_D:
		case MODE_T:
			/* PD or tangent bif. */
			Fpos = NDE + 1;
			break;
		case MODE_F:
			Fpos = NDE;
	}

	for (i = 0, j = 0; i < NDE; i++){
		if (i == sys->section) continue;
		prev[j] = sys->x[i];
		j++;
	}
	prev[CLM_T] = sys->tau;
	prev[CLM_L] = sys->param[sys->variable];


	/* make Dh */
	for (i = 0; i < NDE - 1; i++){
		for (j = 0, k = 0; j < NDE; j++){
			Dh[i][j] = 0.0; 
			if (j == sys->section) continue;
			if (k == i){ Dh[i][j] = 1.0; }
			k++;
			
		}
	}

	/* make Dhinv by transpose of Dh */
	for (i = 0; i < NDE - 1; i++){
		for (j = 0; j < NDE; j++){ 
			Dhinv[j][i] = Dh[i][j];
		}
	}

	/* make Dq */
	for (i = 0; i < NDE; i++){
		if (i == sys->section) Dq[i] = 1.0;
		else Dq[i] = 0.0;
	}

	/* assuming */
	sys->x[sys->section] = sys->xe;

	/* set variables */
	sysvar(sys);

	/* check */
	for (i = 0; i < NDE; i++){
		for (j = i; j < NDE; j++){
			for (k = 0 ; k < NDE; k++){
				double diff;
				diff = fabs(sys->x[NE1 + i * NDE * NDE + j * NDE + k]
				 - sys->x[NE1 + j * NDE * NDE + i * NDE  + k]);
				if (diff > 1e-7){
					printf("something wrong with 2nd variations\n");
					printf("i = %d, j = %d, k = %d : val = %.15f\n", 
						i, j, k, diff);
				}
			}
		}
	}

	k = NDE;
	for (i = 0; i < NDE; i++){
		for (j = 0; j < NDE; j++){
			dphidx[j][i] = sys->x[k];
			k++;
		}
	}

#ifdef EXMAT
	printf("DphiDx:\n");
	for (i = 0; i < NDE; i++){
		for (j = 0; j < NDE; j++){
			printf("%lf\t", dphidx[i][j]);
		}
		printf("\n");
	}
	printf("detA = %.12f\n", det(dphidx, NDE));
#endif

/* possibly for mode_f, below codes should be avoided.*/

	for (i = 0; i < NDE; i++){
		dphidl[i] = sys->x[k];
		k++;
	}

	for (i = 0; i < NDE; i++){
		for (j = 0; j < NDE ; j++){
			sys->Dx[i][j] = dphidx[i][j];
		}
	}

	matpromatmn(Dh, dphidx, buf, NDE-1, NDE, NDE);
	matpromatmn(buf, Dhinv, Dx, NDE-1, NDE, NDE-1);

	for (i = 0; i < NDE; i++){
		for (j = 0; j < NDE; j++){
			NDx[i][j] = Dx[i][j];
		}
	}

	/* make F */
	function(sys->x, f, 0.0, sys); 		
	for (i = 0; i < NDE; i++ ) f0[i] = f[i];

#ifdef DEBUG
	for (i = 0; i < NDE; i++ ) printf("f[%d] = %.12f\n", i, f0[i]);
#endif

	/* set variables to Jacobian matrix of Newton's method */
	/* DT_l - I */
	for (i = 0; i < NDE-1; i++){
		for (j = 0; j < NDE-1; j++){
			if (i == j) NDx[i][j] -= 1.0;
		}
	}


	/* set the rest of Ndx[ROW_Q] */
	vecpromat(Dq, dphidx, bufv, NDE, NDE);
	vecpromat(bufv, Dhinv, NDx[ROW_Q], NDE, NDE-1);

	/* complete column-lambda */
	matprovec(Dh, dphidl, bufv, NDE);
	for (i = 0; i < NDE - 1; i++){ NDx[i][CLM_L] = bufv[i]; }
	NDx[ROW_Q][CLM_L] = inner(Dq, dphidl, NDE);

	matprovec(Dh, f0, bufv, NDE);
	for (i = 0; i < NDE - 1; i++){ NDx[i][CLM_T] = bufv[i]; }
	NDx[ROW_Q][CLM_T] = inner(Dq, f0, NDE);

	/* -F() */ 
	for (i = 0, j = 0; i < NDE; i++){
		if (i != sys->section){ 
			NDx[j][Fpos] = prev[j] - sys->x[i];
#ifdef DEBUG
			printf("prev[%d] - sys->x[%d]\n = %lf - %lf = %lf\n",
				j, i, prev[j], sys->x[i], prev[j] - sys->x[i]);
			printf("NDx[%d][%d] = %lf\n", j, Fpos, NDx[j][Fpos]);
#endif
			j++;
		} 
	}

	/* q(x) = 0 */
	NDx[ROW_Q][Fpos] = -(sys->x[sys->section] - sys->xe);

#ifdef DEBUG
			printf("prev[%d] - sys->xe = %lf - %lf = %lf\n",
				ROW_Q, sys->x[sys->section], sys->xe, 
				sys->x[sys->section] - sys->xe );
			printf("NDx[%d][%d] = %lf\n", ROW_Q, Fpos, NDx[ROW_Q][Fpos]);
#endif

	if (sys->on_real_axis){
		double mu;

		if (sys->mode == MODE_D || (sys->mode == MODE_P && 
			fabs(sys->sigma + 1.0) < EPS)) {	

			mu = -1.0;
		/* ======= CHI ==========*/
		matcopy(dphidx, P, NDE);

		for (i = 0; i < NDE; i++){        /* make chi */
		    for (j = 0; j < NDE; j++){
			if (i == j) { P[i][j] -= mu; }
		    }
		}
		matcopy(P, Pbak, NDE);            /* save chi matrix*/
		chi = det(P, NDE);

		/* ======= DCHIDU ==========*/

		for (i = 0; i < NDE; i++){
	   		dchidx[i] = 0.0;
		    for (j = 0 ; j < NDE; j++){
				matcopy(Pbak, P, NDE);
				for (k = 0; k < NDE; k++){
				    P[k][j] = sys->x[NE1 + NDE * NDE * j + NDE * i + k];    
				}
				dchidx[i] += det(P, NDE);
				matcopy(Pbak, P, NDE);    /* <- do not forget this! */
#ifdef DEBUG
				{	int i, j;
					for (i = 0; i < NDE; i++){	
						for (j = 0; j < NDE; j++){	
							printf("%lf ", P[i][j]);
						}
						printf("\n");
					}
				}
#endif
		    }
		}

#ifdef DEBUG
		for (i = 0; i < NDE; i++){
			printf(">>dchidx[%d] = %.12f\n", i, dchidx[i]);
		}
#endif

		/* matpromatmn(dchidx, Dhinv, dchidu, 1, NDE, NDE-1);*/

		vecpromat(dchidx, Dhinv, dchidu, NDE, NDE-1);

		for (i = 0; i < NDE-1; i++){
		    NDx[NDE][i] = dchidu[i];
		}

		/* ======= DCHIDL ==========*/

		dchidl = 0.0;
		for (j = 0 ; j < NDE; j++){
		    for (k = 0; k < NDE; k++){
				P[k][j] = sys->x[NE1 + NVE2 + NDE * j + k];    
		    }
		    dchidl += det(P, NDE);
		    matcopy(Pbak, P, NDE);
		}
		NDx[NDE][CLM_L] = dchidl;

		/* ======= DCHIDT ==========*/
		/* using values of f() */

		dchidt = 0.0;
		for (j = 0 ; j < NDE; j++){
			for (k = 0; k < NDE; k++){ P[k][j] = f[NDE + NDE * j + k]; }
			dchidt += det(P, NDE);
			matcopy(Pbak, P, NDE);
		}
		NDx[NDE][CLM_T] = dchidt;

	} else if (sys->mode == MODE_T || (sys->mode == MODE_P && 
			fabs(sys->sigma - 1.0) < EPS)) {	
		/* tangent bifurcation */

/************************** TANGENT BIF **************************/
#ifdef DEBUG
	printf("TANGENT!\n");
#endif

		mu = 1.0;

		/* ======= CHI ==========*/

		matcopy(dphidx, P, NDE);
		chi = 0.0;

		for (loop = 0; loop < NDE; loop++){
			matcopy(dphidx, P, NDE);
		    for (i = 0; i < NDE; i++){		/* make chi */
				for (j = 0; j < NDE; j++){
				    if (i == j) { P[i][j] -= mu; }
				}
		    }
		    for (j = 0; j < NDE; j++){
				if (loop == j){ P[j][loop] = -1.0;}
				else P[j][loop] = 0.0;
		    }
		    matcopy(P, PBAK[loop], NDE);			/* save chi matrix*/
			chi += det(P, NDE);
		}
#ifdef DEBUG
		printf("chi(tan) = %lf\n", chi);
#endif

		/* ======= DCHIDU ==========*/

		for (loop = 0; loop < NDE; loop++){
		    dchidx[loop] = 0.0;
		    for (i = 0; i < NDE; i++){
				for (j = 0 ; j < NDE; j++){
				    matcopy(PBAK[i], P, NDE);
				    if (i == j){ 
						for (k = 0; k < NDE; k++){
						    if (k == j){ P[k][j] = 0.0; } 
						}
				    } else {
						for (k = 0; k < NDE; k++){
						    P[k][j] = 
						    	sys->x[NE1 + NDE * NDE * j 
						    		+ NDE * loop + k];	
						}
#ifdef XXX
		printf("--> P\n");
		{	int i, j;
			for (i = 0; i < NDE; i++){	
				for (j = 0; j < NDE; j++){	
					printf("%lf ", P[i][j]);
				}
				printf("\n");
			}
		}
	
#endif
						dchidx[loop] += det(P, NDE);
				    }
				}
		    }
		}

#ifdef DEBUG
	for (i = 0; i < NDE; i++)
		printf("dchidx[%d] = %lf\n", i, dchidx[i]);
#endif
	
		/* matpromatmn(dchidx, Dhinv, dchidu, 1, NDE, NDE-1);*/

		vecpromat(dchidx, Dhinv, dchidu, NDE, NDE-1);
	
		for (i = 0; i < NDE-1; i++){
		    NDx[NDE][i] = dchidu[i];
		}

		/* ======= DCHIDL ==========*/
	
		dchidl = 0.0;
		for (i = 0; i < NDE; i++){
			for (j = 0 ; j < NDE; j++){
				matcopy(PBAK[i], P, NDE);
			
				if (i == j){
		   			for (k = 0; k < NDE; k++){
						if (k == j){ P[k][j] = 0.0;}
					}
				} else {
				    for (k = 0; k < NDE; k++){
						P[k][j] = 
						sys->x[NE1 + NVE2 + NDE * j + k];	
		    		}
#ifdef DEBUG
	{	int ii, jj;
		for (ii = 0; ii < NDE; ii++){
			for (jj = 0; jj < NDE; jj++){
				printf("---------> %.5f ", P[ii][jj]);
			}
			printf("\n");
		}
	}
#endif		
				    dchidl += det(P, NDE);
				}
	    	}
		}
#ifdef DEBUG
	printf("dchidl_tan = %lf\n", dchidl);
#endif
		NDx[NDE][CLM_L] = dchidl;

		/* ======= DCHIDT ==========*/
		/* using values of f() */

			dchidt = 0.0;
			for (i = 0; i < NDE; i++){
				for (j = 0 ; j < NDE; j++){
					matcopy(PBAK[i], P, NDE);
					if (i == j){ 
			   			P[k][j] = 0.0;
					} else {
						for (k = 0; k < NDE; k++){ 
							P[k][j] = f[NDE + NDE * j + k]; 
						}
			   		 	dchidt += det(P, NDE);
					}
				 }
			}	
			NDx[NDE][CLM_T] = dchidt;
		}
		NDx[NDE][NDE+1] = -chi;

/****************************** NS BIF *******************************/
		
	} else if (sys->mode == MODE_N){

		/* NS bif. */
		dets = ExpandDet(sys->Dx, sys->sigma, sys->omega);
		ComplexDet(sys->Dx, dets, &re, &im);

		NDx[ROW_CHIR][Fpos] = -re;
		NDx[ROW_CHII][Fpos] = -im;

		/* assume that phi is not related to the argument of mu */ 
		for (i = 0; i < NDE; i++){ 
			NDx[i][Fpos-1] = 0.0;
		}

		for (i = 0, j = 0; i < NDE; i++){
			if (i == sys->section) continue;
			DerivDetByX(dets, i, &re, &im, sys);
			NDx[ROW_CHIR][j] = re; 
			NDx[ROW_CHII][j] = im;
			++j;
		}
		DerivDetByL(dets, &re, &im, sys);
		NDx[ROW_CHIR][CLM_L] = re; NDx[ROW_CHII][CLM_L] = im;
	
		DerivDetByTau(dets, f, &re, &im, sys);
		NDx[ROW_CHIR][CLM_T] = re; NDx[ROW_CHII][CLM_T] = im;
	
		DerivDetByTheta(dets, &re, &im, sys);
		NDx[ROW_CHIR][CLM_TH] = re; NDx[ROW_CHII][CLM_TH] = im;
	}

#ifdef EXMAT
	printf("\n");
	printf("extended NDx:\n");
	for (i = 0; i < Fpos; i++){
		for (j = 0; j < Fpos + 1 ; j++){ printf("%.8f ",NDx[i][j]); }
		printf("\n");
	}
	fflush(stdout);
#endif
	
	gauss(Fpos, NDx, h);

#ifdef DEBUG_H
	for (i = 0; i < Fpos; i++) printf("h[%d] = %.12f\n", i, h[i]);
	if (!sys->on_real_axis){
		printf("theta: %.10f -> %.10f\n", sys->theta, prev[5]+h[5]);
	}
#endif

	/* poincare section */
	sys->x[sys->section] = sys->xe; 
	for (i = 0, j = 0; i < NDE; i++){
		if (i == sys->section) continue;
#ifdef DEBUG_H
		printf("x[%d]: %.12f -> %.12f\n", i, prev[j], prev[j] + h[j]);
#endif
		sys->x[i] = prev[j] + h[j];
		++j;
	}
	sys->tau  = prev[CLM_T] + h[CLM_T];
#ifdef DEBUG_H
		printf("tau: %.12f -> %.12f\n", prev[CLM_T], prev[CLM_T] + h[CLM_T]);
		printf("para: %.12f -> %.12f\n", sys->param[sys->variable],  
			prev[CLM_L] + h[CLM_L]);
#endif



	if (sys->mode == MODE_P) sys->tau *= 2.0;

/*
	if (sys->mode != MODE_F || sys->mode != MODE_E){
		sys->param[sys->variable] =  prev[CLM_L] + h[CLM_L];
	}
*/

	switch (sys->mode){
		case MODE_F:
		case MODE_E:
			break;
		default:
			sys->param[sys->variable] =  prev[CLM_L] + h[CLM_L];
	}
			

	if (!sys->on_real_axis){
		sys->theta = prev[NDE+1] + h[NDE+1];
	    FreeDecompMat(dets);
	}
}



void newtonequil(SysData *sys)
{
	/* very simple! */

	int i, j;
	double f[NVALL];
	double prev[NNM], h[NNM];
	double NDx[NNM][NNM];
	double *x;

	x = sys->x;

	for (i = 0; i < NDE; i++){ prev[i] = x[i]; }

	function(sys->x, f, 0.0, sys); 		

	for (i = 0; i < NDE; i++){
		for (j = 0; j < NDE; j++){
			NDx[i][j] = sys->Df[i][j];	
		}
	}
			
	for (i = 0; i < NDE; i++){ NDx[i][NDE] = -f[i]; }

#ifdef EXMAT
	printf("DF\n");
	for (i = 0; i < NDE; i++){
		for (j = 0; j < NDE+ 1; j++){
			printf(" %lf ", NDx[i][j]);
		}
		printf("\n");
	}
#endif

	gauss(NDE, NDx, h);

	for (i = 0; i < NDE; i++){ sys->x[i] = prev[i] + h[i]; }
}			



void newtonequilbif(SysData *sys)
{
	double f[NVALL], prev[NNM], h[NNM];
	double NDx[NNM][NNM];
	double *x;
	double re, im;
	double P[NDE][NDE];
	double Pbak[NDE][NDE];
	DecompMat *dets;

	int i, j, k;

	x = sys->x;

	for (i = 0; i < NDE; i++){
		prev[i]	 = x[i];
	}
	prev[NDE] = sys->param[sys->variable];
	if (sys->mode == MODE_H){
		prev[NDE+1] = sys->omega;
	}

	function(sys->x, f, 0.0, sys); 		

	for (i = 0; i < NDE; i++){
		for (j = 0; j < NDE; j++){
			NDx[i][j] = sys->Df[i][j];	
		}
	}
	for (i = 0; i < NDE; i++){
		NDx[i][NDE] = sys->Dl[i];
	}

	for (i = 0; i < NDE; i++){
		NDx[i][NDE+1] = 0.0;
	}
	
	if (sys->mode == MODE_H){

		for (i = 0; i < NDE; i++){
			NDx[i][NDE+2] = -f[i];
		}
		dets = ExpandDet(sys->Df, sys->sigma, sys->omega);
		ComplexDet(sys->Df, dets, &re, &im);

	
		NDx[NDE][NDE+2] = -re;
		NDx[NDE+1][NDE+2] = -im;
	
		for (i = 0; i < NDE; i++){
			DerivDetEqByX(dets, i, &re, &im, sys);
			NDx[NDE][i] = re; 
			NDx[NDE+1][i] = im; 
		}

		DerivDetEqByL(dets, &re, &im, sys);
		NDx[NDE][NDE] = re; 
		NDx[NDE+1][NDE] = im; 
	
		DerivDetEqByOmega(dets, &re, &im, sys);
		NDx[NDE][NDE+1] = re; 
		NDx[NDE+1][NDE+1] = im; 

	} else {

		/* tangent bifurcation of an equilibrium */
		/* MODE_ET */

		double chi, mu;
		double dchidx[NDE];
		double dchidl;

		for (i = 0; i < NDE; i++){
			NDx[i][NDE+1] = -f[i];
		}

		mu = 0.0; /* tangent bif */

		/* ======= CHI ==========*/
		matcopy(sys->Df, P, NDE);

		for (i = 0; i < NDE; i++){        /* make chi */
		    for (j = 0; j < NDE; j++){
			if (i == j) { P[i][j] -= mu; }
		    }
		}
		matcopy(P, Pbak, NDE);            /* save chi matrix*/
		NDx[NDE][NDE+1] = det(P, NDE);
	
		/* ======= DCHIDU ==========*/

		for (i = 0; i < NDE; i++){
	   		dchidx[i] = 0.0;
		    for (j = 0 ; j < NDE; j++){
				for (k = 0; k < NDE; k++){
					P[k][j] = sys->DF2[i][k][j];
				}
				dchidx[i] += det(P, NDE);
				matcopy(Pbak, P, NDE);
			}
		}
		for (i = 0; i < NDE; i++){
		    NDx[NDE][i] = dchidx[i];
		}

		/* ======= DCHIDL ==========*/

		dchidl = 0.0;
		for (j = 0 ; j < NDE; j++){
			for (k = 0; k < NDE; k++){
				P[k][j] = sys->DFL[k][j];
			}
		    dchidl += det(P, NDE);
		    matcopy(Pbak, P, NDE);
		}
		NDx[NDE][NDE] = dchidl;
	}

#ifdef EXMAT
		printf("\nDX:\n");
		for (i = 0; i < NDE + 2; i++){
			for (j = 0; j < NDE + 2; j++){
				printf("%.15f ", NDx[i][j]);
			}
			printf("\n");
		}
#endif

	switch (sys->mode){
		case MODE_ET:
			gauss(NDE+1, NDx, h);
		break;
		case MODE_H:
			gauss(NDE+2, NDx, h);
		break;
	}
	
	for (i = 0; i < NDE; i++){
		sys->x[i] = prev[i] + h[i];
	}
#ifdef DEBUG
	for (i = 0; i < NDE+1; i++){
		printf("h[%d] = %.15f\n", i, h[i]);
	}
#endif
	sys->param[sys->variable]= prev[NDE] + h[NDE];

	if (sys->mode == MODE_H){
		sys->omega = prev[NDE+1] + h[NDE+1];
	}
}				
