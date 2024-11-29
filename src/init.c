#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <Rdefines.h>
#include <Rmath.h>


/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */

extern void F77_NAME(idmlikelihood)(double *,int *, int *, double *,int*, double *,double *,double *,
              int *,int *,int *,int *, int *, double *,double*, double *,
              int *,int *,int *,int *,int *,int *, double *,double*, double *, double *,int *,int *,double *);

extern void F77_NAME(idmlikelihoodweib)(double *,int *, int *, double *,int*,
              int *,int *, double *,double*, double *,
              int *,int *,int *,int *,int *,int *, double *,double*, double *, double *,int *,int *,int *,double *);


extern void F77_NAME(derivaweib)(double *,int *, int *, double *,int*,
              int *,int *, double *,double*, double *,
              int *,int *,int *,int *,int *,int *, double *,double*, double *, double *,int *,int *,double *);

extern void F77_NAME(derivaweibdiag)(double *,int *, int *, double *,int*,
              int *,int *, double *,double*, double *,
              int *,int *,int *,int *,int *,int *, double *,double*, double *, double *,int *,int *,double *);

extern void F77_NAME(firstderivaweib)(double *,int *, int *, double *,int*,
              int *,int *, double *,double*, double *,
              int *,int *,int *,int *,int *,int *, double *,double*, double *, double *,int *,int *,double *);


extern void F77_NAME(derivaspline)(double *,int *, int *, double *,int*, double *,double *,double *,
                     int *,int *,int *,int *, int *, double *,double*, double *,
                     int *,int *,int *,int *,int *,int *, double *,double*, double *, double *,int *,double *);

extern void F77_NAME(derivasplinediag)(double *,int *, int *, double *,int*, double *,double *,double *,
                     int *,int *,int *,int *, int *, double *,double*, double *,
                     int *,int *,int *,int *,int *,int *, double *,double*, double *, double *,int *,double *);

extern void F77_NAME(firstderivaspline)(double *,int *, int *, double *,int*, double *,double *,double *,
                     int *,int *,int *,int *, int *, double *,double*, double *,
                     int *,int *,int *,int *,int *,int *, double *,double*, double *, double *,int *,double *);


static const R_FortranMethodDef FortranEntries[] = {
    {"idmlikelihood",(DL_FUNC) &F77_NAME(idmlikelihood),    29},
    {"idmlikelihoodweib",(DL_FUNC) &F77_NAME(idmlikelihoodweib),    24},
    {"derivaweib",(DL_FUNC) &F77_NAME(derivaweib),    23},
    {"derivaweibdiag",(DL_FUNC) &F77_NAME(derivaweib),    23},
    {"firstderivaweib",(DL_FUNC) &F77_NAME(derivaweib),    23},
    {"derivaspline",(DL_FUNC) &F77_NAME(derivaspline),    28},
    {"derivasplinediag",(DL_FUNC) &F77_NAME(derivaspline),    28},
    {"firstderivaspline",(DL_FUNC) &F77_NAME(derivaspline),    28},
    {NULL, NULL, 0}
};


void R_init_SmoothHazardoptim9(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

