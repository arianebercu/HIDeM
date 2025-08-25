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
              int *,int *,int *,int *,int *,int *, double *,double*, double *, double *,int *,int *,double *);

extern void F77_NAME(idmlikelihoodweibtimedep)(double *,int *, int *, double *,int*,
              int *,int *, double *,double*, double *, double *,double*, double *,
              int *,int *,int *,int *,int *,int *,int *,
              int *,int *,int *,int *,int *,int *, double *,double*, double *, double *,int *,double *);


extern void F77_NAME(idmlikelihoodweibtimedepgrid)(double *,int *, int *, double *,int*,
              int *,int *, double *,double*, double *, double *,double*, double *,
              int *,int *,int *,int *,int *,int *,int *, double *,
              int *,int *,int *,int *,int *,int *, double *,double*, double *, double *,int *,double *);




extern void F77_NAME(derivaweib)(double *,int *, int *, double *,int*,
              int *,int *, double *,double*, double *,
              int *,int *,int *,int *,int *,int *, double *,double*, double *, double *,int *,double *);

extern void F77_NAME(derivaweiballpara)(double *,int *, int *, double *,int*,
              int *,int *, double *,double*, double *,
              int *,int *,int *,int *,int *,int *, double *,double*, double *, double *,int *,double *);

extern void F77_NAME(derivaweiballparadiag)(double *,int *, int *, double *,int*,
              int *,int *, double *,double*, double *,
              int *,int *,int *,int *,int *,int *, double *,double*, double *, double *,int *,double *);

extern void F77_NAME(derivaweiballparafirstderiv)(double *,int *, int *, double *,int*,
              int *,int *, double *,double*, double *,
              int *,int *,int *,int *,int *,int *, double *,double*, double *, double *,int *,double *);

extern void F77_NAME(derivaweibdiag)(double *,int *, int *, double *,int*,
              int *,int *, double *,double*, double *,
              int *,int *,int *,int *,int *,int *, double *,double*, double *, double *,int *,double *);

extern void F77_NAME(firstderivaweib)(double *,int *, int *, double *,int*,
              int *,int *, double *,double*, double *,
              int *,int *,int *,int *,int *,int *, double *,double*, double *, double *,int *,double *);


extern void F77_NAME(derivaspline)(double *,int *, int *, double *,int*, double *,double *,double *,
                     int *,int *,int *,int *, int *, double *,double*, double *,
                     int *,int *,int *,int *,int *,int *, double *,double*, double *, double *,int *,double *);

extern void F77_NAME(derivasplinediag)(double *,int *, int *, double *,int*, double *,double *,double *,
                     int *,int *,int *,int *, int *, double *,double*, double *,
                     int *,int *,int *,int *,int *,int *, double *,double*, double *, double *,int *,double *);

extern void F77_NAME(firstderivaspline)(double *,int *, int *, double *,int*, double *,double *,double *,
                     int *,int *,int *,int *, int *, double *,double*, double *,
                     int *,int *,int *,int *,int *,int *, double *,double*, double *, double *,int *,double *);


extern void F77_NAME(derivaweibfirstderiv)(double *,int *, int *, double *,int*,
                int *,int *, double *,double*, double *,
                int *,int *,int *,int *,int *,int *, double *,double*, double *, double *,int *,double *);

 extern void F77_NAME(derivasplinesfirstderiv)(double *,int *, int *, double *,int*, double *,double *,double *,
                        int *,int *,int *,int *, int *, double *,double*, double *,
                        int *,int *,int *,int *,int *,int *, double *,double*, double *, double *,int *,double *);

extern void F77_NAME(derivasplinessecondderiv)(double *,int *, int *, double *,int*, double *,double *,double *,
                        int *,int *,int *,int *, int *, double *,double*, double *,
                        int *,int *,int *,int *,int *,int *, double *,double*, double *, double *,int *,double *);

 extern void F77_NAME(causalidmlikelihoodweib)(double *,int *, int *, double *,int*,
                            int *,int *, double *,double*, double *,
                            int *,int *,int *, int *,int*, int *,int *,int *, double *,double*, double *, double *,int *,double *);

 extern void F77_NAME(derivaweibsecondderiv)(double *,int *, int *, double *,int*,
                int *,int *, double *,double*, double *,
                int *,int *,int *,int *,int *,int *, double *,double*, double *, double *,int *,double *);
static const R_FortranMethodDef FortranEntries[] = {
    {"idmlikelihood",(DL_FUNC) &F77_NAME(idmlikelihood),    29},
    {"idmlikelihoodweib",(DL_FUNC) &F77_NAME(idmlikelihoodweib),    23},
    {"idmlikelihoodweibtimedep",(DL_FUNC) &F77_NAME(idmlikelihoodweibtimedep),    32},
     {"idmlikelihoodweibtimedepgrid",(DL_FUNC) &F77_NAME(idmlikelihoodweibtimedepgrid),    33},
 {"causalidmlikelihoodweib",(DL_FUNC) &F77_NAME(causalidmlikelihoodweib),    24},
    {"derivaweib",(DL_FUNC) &F77_NAME(derivaweib),    22},
    {"derivaweiballpara",(DL_FUNC) &F77_NAME(derivaweiballpara),    22},
 {"derivaweibsecondderiv",(DL_FUNC) &F77_NAME(derivaweibsecondderiv),    22},
    {"derivaweiballparadiag",(DL_FUNC) &F77_NAME(derivaweiballparadiag),    22},
    {"derivaweiballparafirstderiv",(DL_FUNC) &F77_NAME(derivaweiballparafirstderiv),    22},
 {"derivaweibfirstderiv",(DL_FUNC) &F77_NAME(derivaweibfirstderiv),    22},
    {"derivaweibdiag",(DL_FUNC) &F77_NAME(derivaweibdiag),    22},
    {"firstderivaweib",(DL_FUNC) &F77_NAME(firstderivaweib),    22},
    {"derivaspline",(DL_FUNC) &F77_NAME(derivaspline),    28},
    {"derivasplinediag",(DL_FUNC) &F77_NAME(derivasplinediag),    28},
    {"firstderivaspline",(DL_FUNC) &F77_NAME(firstderivaspline),    28},
  {"derivasplinesfirstderiv",(DL_FUNC) &F77_NAME(derivasplinesfirstderiv),    28},
{"derivasplinessecondderiv",(DL_FUNC) &F77_NAME(derivasplinessecondderiv),    28},
    {NULL, NULL, 0}
};


void R_init_HIDeM(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

