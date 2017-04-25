// Automatically generated, editing not advised.
#ifndef R_GAM_H
#define R_GAM_H
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("gam", String)
#else
#define _(String) (String)
#endif

#define FDEF(name)  {#name, (DL_FUNC) &F77_SUB(name), sizeof(name ## _t)/sizeof(name ## _t[0]), name ##_t}
void F77_SUB(lo0)(
		  double *x,
		  double *y,
		  double *z,
		  int *n,
		  int *d,
		  int *p,
		  int *nvmax,
		  double *span,
		  int *degree,
		  int *match,
		  int *nef,
		  double *dof,
		  double *s,
		  double *var,
		  double *beta,
		  int *iv,
		  int *liv,
		  int *lv,
		  double *v,
		  int *iwork,
		  double *work
		  );

static R_NativePrimitiveArgType lo0_t[] = {
  REALSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  INTSXP,
  REALSXP
};
void F77_SUB(lowese)(
		     int *iv,
		     int *liv,
		     int *lv,
		     double *wv,
		     int *m,
		     double *z,
		     double *s
		     );

static R_NativePrimitiveArgType lowese_t[] = {
  INTSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  INTSXP,
  REALSXP,
  REALSXP
};
void F77_SUB(sknotl)(
		     double *x,
		     int *n,
		     double *knot,
		     int *k
		     );

static R_NativePrimitiveArgType sknotl_t[] = {
  REALSXP,
  INTSXP,
  REALSXP,
  INTSXP
};
void F77_SUB(splsm)(
		    double *x,
		    double *y,
		    double *w,
		    int *n,
		    int *match,
		    int *nef,
		    double *spar,
		    double *dof,
		    double *smo,
		    double *s0,
		    double *cov,
		    int *ifcov,
		    double *work
		    );

static R_NativePrimitiveArgType splsm_t[] = {
  REALSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  LGLSXP,
  REALSXP
};
void F77_SUB(bvalus)(
		     int *n,
		     double *knot,
		     double *coef,
		     int *nk,
		     double *x,
		     double *s,
		     int *order
		     );

static R_NativePrimitiveArgType bvalus_t[] = {
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP
};
void F77_SUB(baklo)(
		    double *x,
		    double *y,
		    double *w,
		    int *npetc,
		    int *wddnfl,
		    double *spatol,
		    int *match,
		    double *etal,
		    double *s,
		    double *eta,
		    double *beta,
		    double *var,
		    double *dof,
		    double *qr,
		    double *qraux,
		    int *qpivot,
		    double *effect,
		    int *iv,
		    double *v,
		    int *iwork,
		    double *work
		    );

static R_NativePrimitiveArgType baklo_t[] = {
  REALSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  REALSXP,
  INTSXP,
  REALSXP,
  INTSXP,
  REALSXP
};
void F77_SUB(bakfit)(
		     double *x,
		     int *npetc,
		     double *y,
		     double *w,
		     int *which,
		     double *spar,
		     double *dof,
		     int *match,
		     int *nef,
		     double *etal,
		     double *s,
		     double *eta,
		     double *beta,
		     double *var,
		     double *tol,
		     double *qr,
		     double *qraux,
		     int *qpivot,
		     double *effect,
		     double *work
		     );

static R_NativePrimitiveArgType bakfit_t[] = {
  REALSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  REALSXP,
  REALSXP
};

static R_FortranMethodDef fMethods[] = {
  FDEF(lo0) ,
  FDEF(lowese) ,
  FDEF(sknotl) ,
  FDEF(splsm) ,
  FDEF(bvalus) ,
  FDEF(baklo) ,
  FDEF(bakfit) ,
  {NULL, NULL, 0}
};

void R_init_gam(DllInfo *dll){
  R_registerRoutines(dll, NULL, NULL, fMethods, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
#endif
