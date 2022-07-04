/*

 *  R : A Computer Language for Statistical Data Analysis

 *  Copyright (C) 2001-2   The R Development Core Team.

 *  Copyright (C) 2003     The R Foundation

 *

 *  This program is free software; you can redistribute it and/or modify

 *  it under the terms of the GNU General Public License as published by

 *  the Free Software Foundation; either version 2 of the License, or

 *  (at your option) any later version.

 *

 *  This program is distributed in the hope that it will be useful,

 *  but WITHOUT ANY WARRANTY; without even the implied warranty of

 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the

 *  GNU General Public License for more details.

 *

 *  You should have received a copy of the GNU General Public License

 *  along with this program; if not, write to the Free Software

 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

 */



#ifndef R_MODREG_H

#define R_MODREG_H



#include <R.h>
#include <Rinternals.h>



SEXP R_isoreg(SEXP y);



void BDRksmooth(double *x, double *y, int *n,

		double *xp, double *yp, int *np,

		int *kern, double *bandwidth);



void

loess_raw(double *y, double *x, double *weights, double *robust, int *d,

	  int *n, double *span, int *degree, int *nonparametric,

	  int *drop_square, int *sum_drop_sqr, double *cell,

	  char **surf_stat, double *surface, int *parameter,

	  int *a, double *xi, double *vert, double *vval, double *diagonal,

	  double *trL, double *one_delta, double *two_delta, int *setLf);

void

loess_dfit(double *y, double *x, double *x_evaluate, double *weights,

	   double *span, int *degree, int *nonparametric,

	   int *drop_square, int *sum_drop_sqr,

	   int *d, int *n, int *m, double *fit);

void

loess_dfitse(double *y, double *x, double *x_evaluate, double *weights,

	     double *robust, int *family, double *span, int *degree,

	     int *nonparametric, int *drop_square,

	     int *sum_drop_sqr,

	     int *d, int *n, int *m, double *fit, double *L);

void

loess_ifit(int *parameter, int *a, double *xi, double *vert,

	   double *vval, int *m, double *x_evaluate, double *fit);

void

loess_ise(double *y, double *x, double *x_evaluate, double *weights,

	  double *span, int *degree, int *nonparametric,

	  int *drop_square, int *sum_drop_sqr, double *cell,

	  int *d, int *n, int *m, double *fit, double *L);



void Srunmed(double *y, double *smo,

	     int *n, int *band, int *end_rule, int *debug);



void Trunmed(int *nn,/* = length(data) */

	     int *kk,/* is odd <= nn */

	     const double *data,

	     double *median, /* (n) */

	     int   *outlist,/* (k+1) */

	     int   *nrlist,/* (2k+1) */

	     double *window,/* (2k+1) */

	     int   *end_rule,

	     int   *print_level);



/* Fortran : */



void F77_SUB(lowesw)(double *res, int *n, double *rw, int *pi);

void F77_SUB(lowesp)(int *n, double *y, double *yhat, double *pwgts,

		     double *rwgts, int *pi, double *ytilde);

void F77_SUB(setppr)(double *span1, double *alpha1,

	int *optlevel, int *ism, double *df1, double *gcvpen1);

void F77_SUB(smart)(int *m, int *mu, int *p, int * q, int *n,

		    double *w, double *x, double *y,

		    double *ww, double *smod, int *nsmod, double *sp,

		    int *nsp, double *dp, int *ndp, double *edf);

void F77_SUB(pppred)(int *np, double *x, double *smod,

		     double *y, double *sc);

void F77_SUB(qsbart)(double *penalt, double *dofoff,

		     double *xs, double *ys, double *ws, double *ssw,

		     int *n, double *knot, int *nk, double *coef,

		     double *sz, double *lev, double *crit, int *iparms,

		     double *spar, double *parms, int *isetup,

		     double *scrtch, int *ld4, int *ldnk, int *ier);



void F77_NAME(sbart)

    (double *penalt, double *dofoff,

     double *xs, double *ys, double *ws, double *ssw,

     int *n, double *knot, int *nk, double *coef,

     double *sz, double *lev, double *crit, int *icrit,

     double *spar, int *ispar, int *iter, double *lspar,

     double *uspar, double *tol, double *eps, int *isetup,

     double *xwy, double *hs0, double *hs1, double *hs2,

     double *hs3, double *sg0, double *sg1, double *sg2,

     double *sg3, double *abd, double *p1ip, double *p2ip,

     int *ld4, int *ldnk, int *ier);



void F77_NAME(sgram)(double *sg0, double *sg1, double *sg2, double *sg3,

		     double *tb, int *nb);

void F77_NAME(stxwx)(double *x, double *z, double *w,

		     int *k, double *xknot, int *n, double *y,

		     double *hs0, double *hs1, double *hs2, double *hs3);

void F77_NAME(sslvrg)(double *penalt, double *dofoff,

		      double *x, double *y, double *w, double *ssw, int *n,

		      double *knot, int *nk, double *coef, double *sz,

		      double *lev, double *crit, int *icrit, double *lambda,

		      double *xwy,

		      double *hs0, double *hs1, double *hs2, double *hs3,

		      double *sg0, double *sg1, double *sg2, double *sg3,

		      double *abd, double *p1ip, double *p2ip,

		      int *ld4, int *ldnk, int *info);



void F77_SUB(bvalus)(int *n, double *knot, double *coef,

		     int *nk, double *x, double *s, int *order);

void F77_SUB(supsmu)(int *n, double *x, double *y,

		     double *w, int *iper, double *span, double *alpha,

		     double *smo, double *sc, double *edf);

#endif

