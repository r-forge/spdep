/* Copyright 2001 by Roger S. Bivand. 
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
**/

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>
#define ROFFSET 1

SEXP symtest(SEXP nb, SEXP card, SEXP verbose)
{
	int i, icard, j, k, k1, flag, fstop, n=length(nb), pc=0;
	SEXP ans;
	PROTECT(ans = NEW_LOGICAL(1)); pc++;
	LOGICAL_POINTER(ans)[0] = TRUE;

	fstop = 0;
	for (i=0; i < n; i++) {
	    icard = INTEGER_POINTER(card)[i];
	    for (j=0; j<icard; j++) {
		flag = 0;
		k = INTEGER_POINTER(VECTOR_ELT(nb, i))[j];
		if (k > 0 && k <= n) {
		    for (k1=0; k1<INTEGER_POINTER(card)[k-ROFFSET]; k1++) {
			if (i+ROFFSET == INTEGER_POINTER(VECTOR_ELT(nb,
			    k-ROFFSET))[k1]) flag += 1;
		    }
		    if (flag != 1) {
			fstop++;
			if (LOGICAL_POINTER(verbose)[0] == TRUE)
			    Rprintf("Non matching contiguities: %d and %d\n",
				i+ROFFSET, k);
		    }
		}
	    }
	}
	if (fstop > 0) LOGICAL_POINTER(ans)[0] = FALSE;

	UNPROTECT(pc); /* ans */
	return(ans);
}

