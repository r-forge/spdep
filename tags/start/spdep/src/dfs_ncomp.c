/* Copyright 2001 by Nicholas Lewin-Koh. 
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

#include<R.h>
#include<Rinternals.h>

#define BLACK 1
#define WHITE 0

void dfs(SEXP nblst, SEXP cmpnm, SEXP visited, int curcmp, int nodeid){
  int n,i;

  INTEGER(cmpnm)[nodeid]=curcmp;
  INTEGER(visited)[nodeid]=BLACK;
  n=length(VECTOR_ELT(nblst,nodeid));

  for(i=0;i<n;i++){
    if(INTEGER(visited)[(INTEGER(VECTOR_ELT(nblst,nodeid))[i]-1)]==WHITE){ 
      dfs(nblst, cmpnm, visited, curcmp,
	  INTEGER(VECTOR_ELT(nblst,nodeid))[i]-1 );
    }
  }
}


SEXP g_components(SEXP nblst, SEXP cmpnm){
  int i, curcmp=1, nvert;
  SEXP visited;
  
  nvert=length(nblst);
  PROTECT(visited=allocVector(INTSXP,nvert));
  
  for(i=0; i < nvert; i++){
    INTEGER(visited)[i]=WHITE;
  }

  for(i=0; i < nvert; i++){
    if(INTEGER(visited)[i]==WHITE){ 
      INTEGER(visited)[i]=BLACK;
      if(INTEGER(VECTOR_ELT(nblst,i))[0]==0){
	INTEGER(cmpnm)[i]=curcmp;
	curcmp++;
      }
      else{
	dfs(nblst, cmpnm, visited, curcmp, i);
	curcmp++;
      }
    }    
  }
  UNPROTECT(1);
  return(cmpnm);
}
