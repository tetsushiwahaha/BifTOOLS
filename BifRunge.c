#include "BifSysDepend.h"

static char rcsid[] = "$Id: BifRunge.c,v 1.1 2001/01/18 04:46:07 tetsushi Exp $";

void runge(int n, double h, double x[], double t,  SysData *Dat)
{
    double k1[NVALL], k2[NVALL], k3[NVALL], k4[NVALL];
    double xtemp[NVALL];
    int i;

    function(x,k1,t, Dat);
    for (i=0; i<n; i++) xtemp[i]=x[i]+k1[i]/2.0*h;
    t+=h/2.0;

    function(xtemp,k2,t, Dat);
    for (i=0; i<n; i++) xtemp[i]=x[i]+k2[i]/2.0*h;

    function(xtemp,k3,t, Dat);
    for (i=0; i<n; i++) xtemp[i]=x[i]+k3[i]*h;
    t+=h/2.0;

    function(xtemp,k4,t, Dat);
    for (i=0; i<n; i++) x[i]+=(k1[i]+2.0*(k2[i]+k3[i])+k4[i])/6.0*h;

}

