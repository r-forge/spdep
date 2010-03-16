/* Copyright 2010 by Roger S. Bivand. */

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>

/** */
static int c__1 = 1;


/**
 * Calculate the sum of squared errors term for spatial regression
 * using an environment to hold data
 *
 * @param env pointer to an SEXP environment
 * @param coef current value of coefficient being optimzed
 * 
 * @return double, value of SSE for current coef
 *
 */
SEXP R_ml_sse_env(SEXP env, SEXP coef) {

  SEXP res;
  SEXP y, x, wy, WX;
  int i, j, k, n, p;
  double *yl, *xlq, *xlqyl, *qy;
  double tol=1e-7, cyl, cxlqyl, sse;
  int *jpvt;
  double *work, *qraux;
  char *trans = "T";
  double one = 1.0, zero = 0.0;
  double lambda = NUMERIC_POINTER(coef)[0];
  int pc=0;

  n = INTEGER_POINTER(findVarInFrame(env, install("n")))[0];
  y = findVarInFrame(env, install("y"));
  x = findVarInFrame(env, install("x"));
  p = INTEGER_POINTER(findVarInFrame(env, install("p")))[0];
  wy = findVarInFrame(env, install("wy"));
  WX = findVarInFrame(env, install("WX"));

  yl = (double *) R_alloc(n, sizeof(double));
  xlq = (double *) R_alloc(n*p, sizeof(double));
  qy = (double *) R_alloc(n*p, sizeof(double));
  xlqyl = (double *) R_alloc(p, sizeof(double));
  jpvt = (int *) R_alloc(p, sizeof(int));
  work = (double *) R_alloc(p*2, sizeof(double));
  qraux = (double *) R_alloc(p, sizeof(double));

  for (i=0; i<n; i++) yl[i] = NUMERIC_POINTER(y)[i] - 
    lambda * NUMERIC_POINTER(wy)[i];
  for (i=0; i<n*p; i++) xlq[i] = NUMERIC_POINTER(x)[i] - 
    lambda * NUMERIC_POINTER(WX)[i];

  F77_CALL(dqrdc2)(xlq, &n, &n, &p, &tol, &k, qraux, jpvt, work);
  if (p != k) warning("Q looses full rank");

  for (i=0; i<n*k; i++) qy[i] = 0.0;
  for (i=0; i<k; i++) qy[(i +(n*i))] = 1.0;

  F77_CALL(dqrqy)(xlq, &n, &k, qraux, qy, &k, qy);

  F77_CALL(dgemv)(trans, &n, &k, &one, qy, &n, yl, &c__1, &zero, xlqyl, &c__1);

  cyl = F77_CALL(ddot)(&n, yl, &c__1, yl, &c__1);

  cxlqyl = F77_CALL(ddot)(&k, xlqyl, &c__1, xlqyl, &c__1);

  sse = cyl - cxlqyl;

  PROTECT(res=NEW_NUMERIC(1)); pc++;
  NUMERIC_POINTER(res)[0] = sse;
  UNPROTECT(pc);

  return(res);

}


