/* Copyright 2010 by Roger S. Bivand. */

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>
#define ROFFSET 1

SEXP opt_error_free(SEXP ptr);
SEXP hess_error_free(SEXP ptr);
SEXP hess_lag_free(SEXP ptr);
SEXP opt_error_init();
SEXP hess_error_init();
SEXP hess_lag_init();
SEXP R_ml_sse_env(SEXP env, SEXP coef);
SEXP R_ml1_sse_env(SEXP env, SEXP lambda, SEXP beta);
SEXP R_ml2_sse_env(SEXP env, SEXP rho, SEXP beta);

void opt_error_set(SEXP env);
void hess_error_set(SEXP env);
void hess_lag_set(SEXP env);

SEXP card(SEXP nb);
SEXP listw2dsT(SEXP nbs, SEXP wts, SEXP card, SEXP ncard2);
SEXP listw2dgR(SEXP nbs, SEXP wts, SEXP card, SEXP ncard);
SEXP listw2sn(SEXP nbs, SEXP wts, SEXP card, SEXP ncard);
SEXP dnearneigh(SEXP din1, SEXP din2, SEXP pnte, SEXP p, SEXP test, SEXP lonlat);
SEXP gearyw(SEXP nb, SEXP weights, SEXP x, SEXP card, SEXP zeropolicy, SEXP ftype);
SEXP gsymtest(SEXP nb, SEXP glist, SEXP card);
SEXP spInsiders(SEXP bbbi, SEXP bbbj);
SEXP jcintern(SEXP nb, SEXP weights, SEXP dum, SEXP card);
SEXP lagw(SEXP nb, SEXP weights, SEXP x, SEXP card, SEXP zeropolicy, SEXP naok);
SEXP nbdists(SEXP nb, SEXP x, SEXP np, SEXP dim, SEXP lonlat);
SEXP polypoly(SEXP p1, SEXP n01, SEXP p2, SEXP n02, SEXP snap);
SEXP symtest(SEXP nb, SEXP card, SEXP verbose);
SEXP g_components(SEXP nblst, SEXP cmpnm);

void dfs(SEXP nblst, SEXP cmpnm, SEXP visited, int curcmp, int nodeid);
void compute_gabriel(int *no_nodes, int *g1, int *g2, int *nogab, int *ngaballoc,  double *nodes_xd, double *nodes_yd);
void compute_relative(int *no_nodes, int *g1, int *g2, int *nogab, int *ngaballoc, double *nodes_xd, double *nodes_yd);
void prunemst(int *e1, int *e2, int *ne, int *gr);
void compute_soi(int *no_nodes, int *g1, int *g2, int *noedges, int *noneigh, int *neigh, int *nearneigh, double *rad, double *nodes_xd, double *nodes_yd);

void gcdist(double *lon1, double *lon2, double *lat1, double *lat2, double *dist);
void knearneigh(int *kin, int *pnte, int *p, double *test, int *res, double *dists, int *lonlat);







