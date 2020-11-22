/* kmp.c : longest extension functions */

#include <stdlib.h>
#include "defs.h"

int *computeLIx (char *pattern, int patLen, int dir, int maxPerBnd) 
     /* Roman's modification of computeLIe */
{
  int *LIx = malloc (maxPerBnd*sizeof(int));

  int i=1 ;

  while(i<maxPerBnd)
    {
      int j = 0 ;
      
      while ( i+j < patLen && 
	      pattern[dir*j] == pattern[dir*(i+j)] )
	j++ ;
      
      LIx[i]=j ;

      i++;
    }

  return LIx ;
}


int *computeLOx (char *pattern, int patLen, int dir, int maxPerBnd) 
     /* Roman's modification of computeLOe */
{
  int *LOx = malloc (maxPerBnd*sizeof(int));

  int i=1 ;

  while(i<maxPerBnd)
    {
      int j = 0 ;
      
      while ( j < patLen && 
	      pattern[dir*j] == pattern[dir*(j-i)] )
	j++ ;
      
      LOx[i]=j ;

      i++;
    }

  return LOx ;
}


