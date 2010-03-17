/* Copyright 2010 by Roger S. Bivand. */

#include "spdep.h"

#include <R_ext/Rdynload.h>

static const R_CMethodDef CEntries[]  = {
    {"dfs", &dfs, 5},
    {"compute_gabriel", &compute_gabriel, 7},
    {"compute_relative", &compute_relative, 7},
    {"prunemst", &prunemst, 4},
    {"compute_soi", &compute_soi, 10},
    {"gcdist", &gcdist, 5},
    {"knearneigh", &knearneigh, 7},
    {NULL, NULL, 0}
};

static R_CallMethodDef CallEntries[] = {
    {"opt_error_free", &opt_error_free, 1},
    {"hess_error_free", &hess_error_free, 1},
    {"hess_lag_free", &hess_lag_free, 1},
    {"opt_error_init", &opt_error_init, 0},
    {"hess_error_init", &hess_error_init, 0},
    {"hess_lag_init", &hess_lag_init, 0},
    {"R_ml_sse_env", &R_ml_sse_env, 2},
    {"R_ml1_sse_env", &R_ml1_sse_env, 3},
    {"R_ml2_sse_env", &R_ml2_sse_env, 3},
    {"card", &card, 1},
    {"listw2dsT", &listw2dsT, 4},
    {"listw2dgR", &listw2dgR, 4},
    {"listw2sn", &listw2sn, 4},
    {"dnearneigh", &dnearneigh, 6},
    {"gearyw", &gearyw, 6},
    {"gsymtest", &gsymtest, 3},
    {"spInsiders", &spInsiders, 2},
    {"jcintern", &jcintern, 4},
    {"lagw", &lagw, 6},
    {"nbdists", &nbdists, 5},
    {"polypoly", &polypoly, 5},
    {"symtest", &symtest, 3},
    {"g_components", &g_components, 2},
    {NULL, NULL, 0}
};


void 
#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
R_init_spdep(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);

}



