/* Copyright 2000 by Roger S. Bivand. 
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

SEXP card(SEXP nb)
{
	int i, n=length(nb), pc=0;
	SEXP ans;
	PROTECT(ans = NEW_INTEGER(n)); pc++;

	for (i=0; i < n; i++) {
	    if (INTEGER_POINTER(VECTOR_ELT(nb, i))[0] == 0) 
		INTEGER_POINTER(ans)[i] = 0;
	    else
		INTEGER_POINTER(ans)[i] = length(VECTOR_ELT(nb, i));
	}

	UNPROTECT(pc); /* ans */
	return(ans);
}

