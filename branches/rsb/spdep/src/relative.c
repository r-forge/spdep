/* Copyright 2001 by Nicholas Lewin-Koh. 
* *
* * This program is free software; you can redistribute it and/or modify
* * it under the terms of the GNU General Public License as published by
* * the Free Software Foundation; either version 2 of the License, or
* * (at your option) any later version.
* *
* * This program is distributed in the hope that it will be useful,
* * but WITHOUT ANY WARRANTY; without even the implied warranty of
* * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* * GNU General Public License for more details.
* **/


#include <R.h>
#include <Rmath.h>


static double distance(double x1, double y1, double x2, double y2){

  return(pythag(x1-x2,y1-y2));

}

void compute_relative(int *no_nodes, int *g1, int *g2, int *nogab,
                      double *nodes_xd, double *nodes_yd)
{
  int i,j,l, no_gab=0;
  double rad;

  for(i=0;i<*no_nodes;i++)
    {
      for(j=i+1;j<*no_nodes;j++)
	{
	  rad=distance(nodes_xd[i],nodes_yd[i],nodes_xd[j],nodes_yd[j]);
          /*Rprintf("hi \n");*/
	  for(l=0;l<*no_nodes;l++)
	    {
	      if((l!=i)&&(l!=j)&&
	      (distance(nodes_xd[i],nodes_yd[i],nodes_xd[l],nodes_yd[l])<rad)&&
              (distance(nodes_xd[j],nodes_yd[j],nodes_xd[l],nodes_yd[l])<rad))
		break;
	    }

	  if(l==*no_nodes)
	    {
	      g1[no_gab]=i+1; g2[no_gab++]=j+1;
	    }
	}
    }
  *nogab=no_gab;
  return;

}
