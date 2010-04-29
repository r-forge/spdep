/* Copyright 2004 by Roger S. Bivand. */

#include "spdep.h"

SEXP polypoly(SEXP p1, SEXP n01, SEXP p2, SEXP n02, SEXP snap)
{
	int n1=INTEGER_POINTER(n01)[0], n2=INTEGER_POINTER(n02)[0], pc=0;
	int i, j, k=0;
	double sn=NUMERIC_POINTER(snap)[0], dist;
	double x1, x2, y1, y2, xd, yd, sn2 =sn*sn;

	SEXP ans;
	PROTECT(ans = NEW_INTEGER(1)); pc++;

	for (i=0; (i < n1) && (k < 2); i++) {
		x1 = NUMERIC_POINTER(p1)[i];
		y1 = NUMERIC_POINTER(p1)[n1 + i];
		for (j=0; (j < n2) && (k < 2); j++) {
			x2 = NUMERIC_POINTER(p2)[j];
			y2 = NUMERIC_POINTER(p2)[n2 + j];
/*			dist = pythag((x1-x2), (y1-y2));
			if (dist < sn) k++;
			if (k > 1) break;*/
/* following lines Micah Altman 2010 */
			xd = x1-x2;
			if (fabs(xd)>sn) { continue; }
			yd = y1-y2;
			if (fabs(yd)>sn) { continue; }
			dist = xd*xd + yd*yd;
			if (dist <= sn2) k++;
		}
	}
	
	INTEGER_POINTER(ans)[0] = k;

	UNPROTECT(pc); /* ans */
	return(ans);
}

/* function by Micah Altman */

SEXP spOverlap(SEXP bbbi, SEXP bbbj) {

	int pc=0,overlap=1;
	double bbi[4], bbj[4];
	SEXP ans;

	PROTECT(ans = NEW_INTEGER(1)); pc++;
	bbi[0] = NUMERIC_POINTER(bbbi)[0];
	bbi[1] = NUMERIC_POINTER(bbbi)[1];
	bbi[2] = NUMERIC_POINTER(bbbi)[2];
	bbi[3] = NUMERIC_POINTER(bbbi)[3];
	bbj[0] = NUMERIC_POINTER(bbbj)[0];
	bbj[1] = NUMERIC_POINTER(bbbj)[1];
	bbj[2] = NUMERIC_POINTER(bbbj)[2];
	bbj[3] = NUMERIC_POINTER(bbbj)[3];

        if ((bbi[0]>bbj[2]) | (bbi[1]>bbj[3]) | 
		(bbi[2]<bbj[0]) | (bbi[3]<bbj[1]) ) {
		overlap=0;
	}

	INTEGER_POINTER(ans)[0] = overlap;		
	UNPROTECT(pc); /* ans */
	return(ans);
}


