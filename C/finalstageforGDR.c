/* Search for GDRs of the second type and assembling GDRs into runs */

#include <stdlib.h>
#include "defs.h"

#define FNDRUNS
void EmergSortFirstTypeRepeats(void);
void FindSecondTypeRepeats(void);
void FindSecondTypeExactRepeats(void);
void CopyFactorRepeats(void);
void CopyFactorExactRepeats(void);
void NewSecondTypeRepeat(int repeatPer, int repeatEnd);
void SortAllRepeats(void);
void FindAllRuns(void);
void FindAllStepRuns(void);
void AssemblRuns(void);
void AssemblStepRuns(void);
void FindRepEnd(int *endPointer);
int FindRepStart(int repInitPosition);

int prevHeadPstn, nextHeadPstn, minGapLen, curFact, startRep, maxInitPstn;
listreps *presortLists;    // repetitions of the first type sorted by end positions
listreps *sortedLists;     // repetitions of both types sorted by start positions
listreps *nextRepPtr, *perLists, *maxIntPtr, *initLists, curListRep;
transreps maxTransRep=0, nxtTransRep, *newTransRep;
int  actTransPrd, dblPrd, repPrd;
extern int wordLength, dblMinPer, minPer, maxPer, actPerBnd, tooBigReps;
/* extern int *numRep; */
extern int *copyEnds, *factorEnds, *copyBgns, *factorStarts, factorNum;
extern char *seq_original, prevLetterCheck, nextLetterCheck;
extern int step, LastWindow, maxNumErr;
char *seq_original2;



void EmergSortFirstTypeRepeats(void)
     /* analogous to SortFirstTypeRepeats, for the case when
	there is no enough memory for 'sortedLists' */
{  
  int endPstn, initPstn, maxEndPstn;
  listreps curListRep, nxtListRep, prvListRep, *initPtr;

  sortedLists=presortLists;
  maxEndPstn= wordLength + 1;
  for(endPstn=dblMinPer;endPstn<=maxEndPstn;endPstn++)
    {  
      if ((curListRep=presortLists[endPstn]))
	{  
	  presortLists[endPstn]=NULL;
	  do
	    {  
	      initPstn=curListRep->rep.initpos;
	      nxtListRep=curListRep->next;
	      curListRep->next=sortedLists[initPstn];
	      sortedLists[initPstn]=curListRep;
	    }  
	  while((curListRep=nxtListRep));
	}
    }
  maxIntPtr=sortedLists+(wordLength - dblMinPer);
  for(initPtr=sortedLists;initPtr<=maxIntPtr;initPtr++)
    {  
      if ((prvListRep=(*initPtr)))
	{  
	  if ((curListRep=prvListRep->next))
	    {  
	      prvListRep->next=NULL;
	      while((nxtListRep=curListRep->next))
		{  
		  curListRep->next=prvListRep;
		  prvListRep=curListRep;
		  curListRep=nxtListRep;
		} 
	      curListRep->next=prvListRep;
	      *initPtr=curListRep;
	    }
	}
    }
}


void FindSecondTypeRepeats(void)
     /* searching for approximate reps of the second type */
{  
  minGapLen=dblMinPer+2;
  nextHeadPstn=1;
  for(curFact=2; curFact <= factorNum; curFact++)
    {  
      prevHeadPstn=nextHeadPstn;
      nextHeadPstn=factorEnds[curFact];
      if (nextHeadPstn-prevHeadPstn>minGapLen)
	CopyFactorRepeats();
    }
}


void FindSecondTypeExactRepeats(void)
     /* searching for exact reps of the second type */
{  
  minGapLen=dblMinPer+1;
  nextHeadPstn=1;
  for(curFact=1; curFact < factorNum; curFact++)
    {  
/*       if ( (prevHeadPstn=nextHeadPstn) >= step ) */
      prevHeadPstn=nextHeadPstn;
      if ( step > 0 && prevHeadPstn >= VIRTUAL_STEP )  // GK
	break;
      nextHeadPstn=factorStarts[curFact + 1];
      if (nextHeadPstn-prevHeadPstn>minGapLen)
	CopyFactorExactRepeats();
    }
  free(copyBgns);
  free(factorStarts);
}


void CopyFactorRepeats(void)
     /* computing all approximate reps of the second type, 
	contained in the current LZ-factor in the current window */
{
  listreps copyListRep, initFstRep;
  repetition *copyRep;
  int copyRepEnd, curCopyEnd, shift, maxStart;

  curCopyEnd=copyEnds[curFact];
  shift=nextHeadPstn-curCopyEnd;
  maxStart=nextHeadPstn-minGapLen;
  for(startRep=prevHeadPstn+1;startRep<=maxStart;startRep++)
    {  
      nextRepPtr=sortedLists+startRep;
      if ( (copyListRep=(*(nextRepPtr-shift))) )
	{  
	  initFstRep=(*nextRepPtr);
	  do
	    {  
	      copyRep=&(copyListRep->rep);
	      if ( (copyRepEnd=copyRep->endpos)<curCopyEnd )
		NewSecondTypeRepeat(copyRep->period, copyRepEnd+shift);
	      else
		break;
	    }  
	  while( (copyListRep=copyListRep->next) );
	  *nextRepPtr=initFstRep;
	}
    }
}


void CopyFactorExactRepeats(void)
     /* computing all exact reps of the second type,
	contained in the current s-factor in the current window */
{
  listreps copyListRep, initFstRep;
  repetition *copyRep;
  int copyRepEnd, curCopyEnd, shift, maxStart;

  shift=prevHeadPstn-copyBgns[curFact];
  curCopyEnd=nextHeadPstn-shift;
  maxStart=nextHeadPstn-minGapLen;
  for(startRep=prevHeadPstn+1;startRep<=maxStart;startRep++)
    {  
      nextRepPtr=sortedLists+startRep;
      if ( (copyListRep=(*(nextRepPtr-shift))) )
	{  
	  initFstRep=(*nextRepPtr);
	  do
	    {  
	      copyRep=&(copyListRep->rep);
	      if ( (copyRepEnd=copyRep->endpos)<curCopyEnd )
		NewSecondTypeRepeat(copyRep->period, copyRepEnd+shift);
	      else
		break;
	    }  
	  while( (copyListRep=copyListRep->next) );
	  *nextRepPtr=initFstRep;
	}
    }
}


void NewSecondTypeRepeat(int repeatPer, int repeatEnd)
     /* called from CopyFactorRepeats and CopyFactorExactRepeats;
	processing reps of the second type. 
	Creates for a found rep a corresponding structure 'listreps'
	and inserts it into the corresponding list from 'sortedLists' */
{  
  repetition *newRep;

   if ( (*nextRepPtr=(listreps)calloc(1, sizeof(struct s_listreps)))==NULL )
   {  
     printf("The output is too big.\n");
      printf("Narrow the period bounds.\n");
      exit(15);
   }
   newRep=&((*nextRepPtr)->rep);
   newRep->initpos=startRep;
   newRep->endpos=repeatEnd;
   newRep->period=repeatPer;
   nextRepPtr=&((*nextRepPtr)->next);
}


void SortAllRepeats(void) 
     /* sorting approximate reps taken from 'sortedLists'
	into 'perLists' (array of lists of found reps 
	sorted by their periods) */
{
  listreps curListRep, nextListRep, *initPtr;
  if (actPerBnd>maxPer+1)
    actPerBnd=maxPer+1;
  for(initPtr=maxIntPtr;initPtr>=sortedLists;initPtr--)
    {  
      if ((curListRep=(*initPtr)))
      {  
	do
	{
	  nextListRep=curListRep->next;
	  curListRep->next=perLists[curListRep->rep.period];
	  perLists[curListRep->rep.period]=curListRep;
	}  
	while((curListRep=nextListRep));
      }
    }
  free(sortedLists);  
}


void EmergSortAllRepeats(void)
     /* analogous to 'SortAllRepeats' in the case when
	there is no enough memory for 'perLists' */
{  
  int repPer;
  listreps curListRep, nextListRep, *initPtr;
  if (actPerBnd>maxPer+1)
    actPerBnd=maxPer+1;
  perLists= sortedLists + wordLength + 1;
  for(initPtr=maxIntPtr;initPtr>=sortedLists;initPtr--)
    {  
      if ((curListRep=(*initPtr)))
	{  
	  *initPtr=NULL;
	  do
	    {  
	      repPer=curListRep->rep.period;
	      nextListRep=curListRep->next;
	      curListRep->next=(*(perLists-repPer));
	      *(perLists-repPer)=curListRep;
	    }  
	  while((curListRep=nextListRep));
	}
    }
  for(repPer=minPer;repPer<actPerBnd;repPer++)
    sortedLists[repPer]=(*(perLists-repPer));
  perLists=(listreps *)realloc(sortedLists, actPerBnd*sizeof(listreps));
}


void FindAllRuns(void)
     /* assembling GDRs stored in 'perLists' into runs,
	and storing these runs into array of lists 'initLists'
	in the case when the current window is the last one
	('initLists' - array of lists of runs of the current window,
	sorted by start positions) */

     /* case where investigated fragment is completely inside the window */
{ transreps remTransRep;
  listreps remListRep;
  int remRepEnd;
  maxInitPstn= wordLength - minPer;
  if ( (initLists=(listreps *)calloc(maxInitPstn, sizeof(listreps)))==NULL )
    {  
      printf("Not enough memory for checking the repeats\n");
      exit(4);
    }
  seq_original2=seq_original-2;
  if (maxTransRep)
  {  actTransPrd=maxTransRep->period;
     nxtTransRep=maxTransRep->next;
     free(maxTransRep);
  }
  else
     actTransPrd=0;
  for(repPrd=actPerBnd-1; repPrd>=minPer; repPrd--)
   {  if ((curListRep=perLists[repPrd]))
      {  dblPrd=2*repPrd;
         if (repPrd==actTransPrd)
         {  if (nxtTransRep)
            {  actTransPrd=nxtTransRep->period;
               nxtTransRep=(remTransRep=nxtTransRep)->next;
               free(remTransRep);
            }
            else
               actTransPrd=0;
            remRepEnd=curListRep->rep.endpos;
            curListRep=(remListRep=curListRep)->next;
            free(remListRep);
            while(curListRep)
            {  if ( curListRep->rep.initpos<=remRepEnd-dblPrd )
               {  remRepEnd=curListRep->rep.endpos;
                  curListRep=(remListRep=curListRep)->next;
	          free(remListRep);
               }
               else
               {  AssemblRuns();
                  break;
               }
            }
         }
         else
            AssemblRuns();
      }
   }
  free(perLists);
}

void AssemblRuns()
{  int  repIntPos, *curRepEnd;
   listreps nxtListRep;
   repIntPos=FindRepStart(curListRep->rep.initpos);
   nxtListRep=curListRep->next;
   curListRep->next=initLists[repIntPos];
   initLists[repIntPos]=curListRep;
   curRepEnd=&(curListRep->rep.endpos);
   while((curListRep=nxtListRep))
   {  if ( curListRep->rep.initpos<=(*curRepEnd)-dblPrd )
      {  *curRepEnd=curListRep->rep.endpos;
	 nxtListRep=curListRep->next;
	 free(curListRep);
      }
      else
      {  FindRepEnd(curRepEnd);
         repIntPos=FindRepStart(curListRep->rep.initpos);
	 nxtListRep=curListRep->next;
         curListRep->next=initLists[repIntPos];
	 initLists[repIntPos]=curListRep;
         curRepEnd=&(curListRep->rep.endpos);
      }
   }
   FindRepEnd(curRepEnd);
}

void FindAllStepRuns(void)
     /* analogous to FindAllRuns for the case when 
	the current window is not the last one */

     /* case where investigated fragment extends beyond the window */
     /* called only if step is specified and smaller than wordLength */
{ transreps remTransRep;
  listreps remListRep;
  int remRepEnd;
  maxInitPstn=step+actPerBnd;
  if ( (initLists=(listreps *)calloc(maxInitPstn, sizeof(listreps)))==NULL )
    {  
      printf("Not enough memory for checking the repeats\n");
      exit(4);
    }
  seq_original2=seq_original-2;
  if (maxTransRep)
  {  actTransPrd=maxTransRep->period;
     nxtTransRep=maxTransRep->next;
     free(maxTransRep);
     maxTransRep=0;
  }
  else
     actTransPrd=0;
  newTransRep=&maxTransRep;
  for(repPrd=actPerBnd-1; repPrd>=minPer; repPrd--)
   {  if ((curListRep=perLists[repPrd]))
      {  dblPrd=2*repPrd;
	 if (repPrd==actTransPrd)
	 {  if (nxtTransRep)
            {  actTransPrd=nxtTransRep->period;
               nxtTransRep=(remTransRep=nxtTransRep)->next;
               free(remTransRep);
            }
            else
               actTransPrd=0;
            remRepEnd=curListRep->rep.endpos;
            curListRep=(remListRep=curListRep)->next;
            free(remListRep);
            while(curListRep)
            {  if ( curListRep->rep.initpos<=remRepEnd-dblPrd )
               {  remRepEnd=curListRep->rep.endpos;
                  curListRep=(remListRep=curListRep)->next;
	          free(remListRep);
               }
               else
               {  AssemblStepRuns();
                  break;
               }
            }
         }
         else
            AssemblStepRuns();
      }
   }
  free(perLists);
}


void AssemblStepRuns(void)
     /* assemble runs with a given period */
{  int  repIntPos, *curRepEnd;
 listreps nxtListRep;
 if ( (repIntPos=curListRep->rep.initpos)<step )
   {  
     repIntPos=FindRepStart(repIntPos);  // find the beginning of the run
     nxtListRep=curListRep->next;
     curListRep->next=initLists[repIntPos];
     initLists[repIntPos]=curListRep;
     curRepEnd=&(curListRep->rep.endpos);
     while((curListRep=nxtListRep))
       {  
	 if ( curListRep->rep.initpos<=(*curRepEnd)-dblPrd )
	   {  
	     *curRepEnd=curListRep->rep.endpos;
	     nxtListRep=curListRep->next;
	     free(curListRep);
	   }
	 else
	   {  
	     if ( (repIntPos=curListRep->rep.initpos)<step )
	       {  
		 FindRepEnd(curRepEnd);  // find the end of the run
		 repIntPos=FindRepStart(repIntPos);  // find the beginning of the following one
		 nxtListRep=curListRep->next;
		 curListRep->next=initLists[repIntPos];
		 initLists[repIntPos]=curListRep;
		 curRepEnd=&(curListRep->rep.endpos);
	       }
	     else
	       {  
		 if (*curRepEnd>step+dblPrd)
		   {  
		     if ( (*newTransRep=(transreps)calloc(1, sizeof(struct s_transreps)))==NULL )
		       {  
			 printf("Not enough memory for this step.\n");
			 exit(115);
		       }
		     (*newTransRep)->period=repPrd;
		     newTransRep=&((*newTransRep)->next);
		   }
		 FindRepEnd(curRepEnd);
		 break;
	       }
	   }
       }
     if (!curListRep)
       {  
	 if (*curRepEnd>step+dblPrd)
	   {  
	     if ( (*newTransRep=(transreps)calloc(1, sizeof(struct s_transreps)))==NULL )
	       {  
		 printf("Not enough memory for this step.\n");
		 exit(116);
	       }
	     (*newTransRep)->period=repPrd;
	     newTransRep=&((*newTransRep)->next);
	     if (*curRepEnd>2*step)
	       {  
		 initLists[repIntPos]=(nxtListRep=initLists[repIntPos])->next;
		 free(nxtListRep);
		 tooBigReps=YES;
		 return;
	       }
	   }
	 FindRepEnd(curRepEnd);
	 return;
       }         
   }
 do
   {  
     nxtListRep=curListRep->next;
     free(curListRep);
   }  
 while((curListRep=nxtListRep));
}


int FindRepStart(int repInitPos)
     /* process the beginning of the run and trim it out if it is highly erroneous */
{  int  curInitPos, actInitPos, initPosBnd, numBgnErr=0;
 curInitPos=repInitPos;
 initPosBnd=repInitPos+repPrd;
 while( seq_original[curInitPos]!= seq_original[curInitPos+repPrd])
   {  
     curInitPos++;
     numBgnErr++;
   }
 actInitPos=curInitPos;
 do
   {  
     do
       {  
	 curInitPos++;
       }  
     while( seq_original[curInitPos]== seq_original[curInitPos+repPrd]);
     do
       {  
	 if ( (++curInitPos)>=initPosBnd )
	   return(actInitPos);
         numBgnErr++;
       }  
     while( seq_original[curInitPos]!= seq_original[curInitPos+repPrd]);
     /* check the "formula" and set the "real beginning" */
     if ( 2*numBgnErr*repPrd > (curInitPos-repInitPos)*(repPrd+maxNumErr) )
       actInitPos=curInitPos;
   }  
 while(YES);
}


void FindRepEnd(int *endPtr)
     /* process the end of the run and trim it out if it is highly erroneous */
{  
  int  EndPosBnd, curEndPos, actEndPos, numEndErr=0;
  curEndPos=(*endPtr);
  EndPosBnd=curEndPos-repPrd;
  while( seq_original2[curEndPos]!= seq_original2[curEndPos-repPrd] )
    {  
      curEndPos--;
      numEndErr++;
    }
  actEndPos=curEndPos;
  do
    {  
      do
	{  
	  curEndPos--;
	}  
      while( seq_original2[curEndPos]==seq_original2[curEndPos-repPrd] );
      do
	{  
	  if ( (--curEndPos)<=EndPosBnd )
	    {  
	      *endPtr=actEndPos;
	      return;
	    }
	  numEndErr++;
	}  
      while( seq_original2[curEndPos]!=seq_original2[curEndPos-repPrd] );
      if ( 2*numEndErr*repPrd >
	   ((*endPtr)-curEndPos)*(repPrd+maxNumErr) )
	actEndPos=curEndPos;
    }  
  while(YES);
}
