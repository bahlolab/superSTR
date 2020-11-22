/* Search for GDRs of the first type */

#include <stdlib.h>
#include "defs.h"

void FindLocalHeadGDR(void);
void FindFinalLocalHeadGDR(void);
void FindHeadGDR(void);
void FindHeadGDR1(void);
void FindHeadGDR2(void);
void FindHeadGDR3(void);
void FindStartGlobalHeadGDR(void);
void FindStartHeadGDR(void);
void FindFinalGlobalHeadGDR(void);
void FindFinalHeadGDR0(void);
void FindFinalHeadGDR1(void);
void FindGlobalHeadGDR0(void);
void FindAllGlobalHeadGDR(void);
void FreeFindAllGlobalHeadGDR(void);
void FastFindAllGlobalHeadGDR(void);
void EasyFindAllGlobalHeadGDR(void);
void FastFreeFindAllGlobalHeadGDR(void);
void HardFindAllGlobalHeadGDR(void);
void FindOnlyLeftGlobalHeadGDR(void);
void ComputeRadFact(void);

void FreeFindAllHeadGDR(void);
void RightFreeFindAllHeadGDR(void);
void FindAllHeadGDR(void);
void FindOnlyRightHeadGDR(void);
void BriefFindAllHeadGDR(void);
void LeftFreeFindAllHeadGDR(void);
void FindOnlyLeftHeadGDR(void);
void NewFirstTypeRepeat(int,int);
void EasyFindRightHeadRuns(void);

extern int Period;
int curMaxPer, curBlock, *curBlockEnds, *radFact;
int prevDist, nextDist, hlfPrvDist, hlfNxtDist;

extern  char *seq_working_copy;
extern int *factorEnds, factorNum, numGlobFact, numFullBlock;
extern int blockSize, finBlockSize, *maxSumFact, actPerBnd;
extern int *subRepEnd, *subRepBgn, *subPrfEnd, *subSufBgn;
extern int curHeadPstn, actHeadPstn, minRepBnd, maxRepBnd;
extern int numPrfErr, numSufErr, maxNumErr;
extern int prfEnd, lastPrfEnd, sufBgn, lastSufBgn, maxPrfEnd, minSufBgn;
extern int wordLength, minPer, maxPer, actPerBnd;
extern listreps *presortLists;


void mainstage(void)
{  subRepBgn=subRepEnd=(int *)calloc(maxNumErr+1, sizeof(int));
   if( ( presortLists=(listreps *) calloc ((2 + wordLength), sizeof(listreps))
	) == NULL )
   {  printf("Not enough memory for saving the repetitions\n");
      exit(100);  }
   ComputeRadFact();
   curBlock=1;
   FindStartGlobalHeadGDR();
   curBlockEnds=factorEnds;
   FindLocalHeadGDR();
   minRepBnd=-1;
   actHeadPstn=curHeadPstn=(*curBlockEnds);
   if (numFullBlock>1)
   {  maxRepBnd=curBlockEnds[blockSize];
      FindHeadGDR();
      while(++curBlock<numFullBlock)
      {  FindLocalHeadGDR();
         curHeadPstn=(*curBlockEnds);
         maxRepBnd=curBlockEnds[blockSize];
         if (curHeadPstn>=2*maxSumFact[curBlock])
            FindGlobalHeadGDR0();
         else     
         {  minRepBnd=-1;
            actHeadPstn=curHeadPstn;
            FindHeadGDR();
         }
      }
      FindLocalHeadGDR();
      curHeadPstn=(*curBlockEnds);
      maxRepBnd=curBlockEnds[finBlockSize];
      if (curHeadPstn>=2*maxSumFact[curBlock])
         FindGlobalHeadGDR0();
      else     
      {  minRepBnd=-1;
         actHeadPstn=curHeadPstn;
         FindHeadGDR();
      }
   }
   else
   {  maxRepBnd=curBlockEnds[finBlockSize];
      FindHeadGDR();
   }
   if (finBlockSize>1)
      FindFinalLocalHeadGDR();
   FindFinalGlobalHeadGDR();
   free(maxSumFact);
}


void FindLocalHeadGDR(void)
{  int fact;
   for(fact=1;fact<blockSize;fact++)
   {  minRepBnd=curBlockEnds[fact-radFact[fact]];
      actHeadPstn=curHeadPstn=curBlockEnds[fact];
      maxRepBnd=curBlockEnds[fact+radFact[fact]];
      FindHeadGDR();
   }
   curBlockEnds+=blockSize;
}


void FindFinalLocalHeadGDR(void)
{  int fact;
   for(fact=1;fact<finBlockSize;fact++)
   {  minRepBnd=curBlockEnds[fact-radFact[fact]];
      actHeadPstn=curHeadPstn=curBlockEnds[fact];
      if (fact+radFact[fact]<finBlockSize)
         maxRepBnd=curBlockEnds[fact+radFact[fact]];
      else
         maxRepBnd=curBlockEnds[finBlockSize];
      FindHeadGDR();
   }
}


void FindHeadGDR(void)
{  prevDist=(curHeadPstn-minRepBnd)-2; 
   nextDist=(maxRepBnd-curHeadPstn)-2;
   if ( (curMaxPer=(prevDist+nextDist+1)/2)>=minPer )
   {  hlfPrvDist=prevDist/2;
      hlfNxtDist=nextDist/2;
      Period=minPer; 
      if (prevDist<nextDist)
         FindHeadGDR1();
      else
      {  if (prevDist>nextDist)
            FindHeadGDR2();
         else
            FindHeadGDR3();
      }
   }
}


void FindHeadGDR1(void)
{  int thrdDist;
   if (hlfPrvDist<maxPer)
   {  while(Period<=hlfPrvDist)
      {  FreeFindAllHeadGDR();
         Period++;  }
   }
   else
   {  while(Period<=maxPer)
      {  FreeFindAllHeadGDR();
         Period++;  }
      return;
   }
   if (hlfNxtDist<=prevDist)
   {  if (hlfNxtDist<maxPer)
      {  while(Period<=hlfNxtDist)
         {  RightFreeFindAllHeadGDR();
            Period++;  }
      }
      else
      {  while(Period<=maxPer)
         {  RightFreeFindAllHeadGDR();
            Period++;  }
         return;
      }
      if ( (thrdDist=(prevDist+nextDist)/3)<maxPer )
      {  while(Period<=thrdDist)
         {  FindAllHeadGDR();
            Period++;  }
      }
      else
      {  while(Period<=maxPer)
         {  FindAllHeadGDR();
            Period++;  }
         return;
      }
   }
   else
   {  if (prevDist<maxPer)
      {  while(Period<=prevDist)
         {  RightFreeFindAllHeadGDR();
            Period++;  }
      }
      else
      {  while(Period<=maxPer)
         {  RightFreeFindAllHeadGDR();
            Period++;  }
         return;
      }
      if (hlfNxtDist<maxPer)
      {  while(Period<=hlfNxtDist)
         {  FindOnlyRightHeadGDR();
            Period++;  }
      }
      else
      {  while(Period<=maxPer)
         {  FindOnlyRightHeadGDR();
            Period++;  }
         return;
      }
   }
   if (curMaxPer<maxPer)
   {  while(Period<=curMaxPer)
      {  BriefFindAllHeadGDR();
         Period++;  }
   }
   else
   {  while(Period<=maxPer)
      {  BriefFindAllHeadGDR();
         Period++;  }
   }
}


void FindHeadGDR2(void)
{  int thrdDist;
   if (hlfNxtDist<maxPer)
   {  while(Period<=hlfNxtDist)
      {  FreeFindAllHeadGDR();
         Period++;  }
   }
   else
   {  while(Period<=maxPer)
      {  FreeFindAllHeadGDR();
         Period++;  }
      return;
   }
   if (hlfPrvDist<=nextDist)
   {  if (hlfPrvDist<maxPer)
      {  while(Period<=hlfPrvDist)
         {  LeftFreeFindAllHeadGDR();
            Period++;  }
      }
      else
      {  while(Period<=maxPer)
         {  LeftFreeFindAllHeadGDR();
            Period++;  }
         return;
      }
      if ( (thrdDist=(prevDist+nextDist)/3)<maxPer )
      {  while(Period<=thrdDist)
         {  FindAllHeadGDR();
            Period++;  }
      }
      else
      {  while(Period<=maxPer)
         {  FindAllHeadGDR();
            Period++;  }
         return;
      }
   }
   else
   {  if (nextDist<maxPer)
      {  while(Period<=nextDist)
         {  LeftFreeFindAllHeadGDR();
            Period++;  }
      }
      else
      {  while(Period<=maxPer)
         {  LeftFreeFindAllHeadGDR();
            Period++;  }
         return;
      }
      if (hlfPrvDist<maxPer)
      {  while(Period<=hlfPrvDist)
         {  FindOnlyLeftHeadGDR();
            Period++;  }
      }
      else
      {  while(Period<=maxPer)
         {  FindOnlyLeftHeadGDR();
            Period++;  }
         return;
      }
   }
   if (curMaxPer<maxPer)
   {  while(Period<=curMaxPer)
      {  BriefFindAllHeadGDR();
         Period++;  }
   }
   else
   {  while(Period<=maxPer)
      {  BriefFindAllHeadGDR();
         Period++;  }
   }
}


void FindHeadGDR3(void)
{  int thrdDist;
   if (hlfPrvDist<maxPer)
   {  while(Period<=hlfPrvDist)
      {  FreeFindAllHeadGDR();
         Period++;  }
   }
   else
   {  while(Period<=maxPer)
      {  FreeFindAllHeadGDR();
         Period++;  }
      return;
   }
   if ( (thrdDist=(prevDist+nextDist)/3)<maxPer )
   {  while(Period<=thrdDist)
      {  FindAllHeadGDR();
         Period++;  }
   }
   else
   {  while(Period<=maxPer)
      {  FindAllHeadGDR();
         Period++;  }
      return;
   }
   if (curMaxPer<maxPer)
   {  while(Period<=curMaxPer)
      {  BriefFindAllHeadGDR();
         Period++;  }
   }
   else
   {  while(Period<=maxPer)
      {  BriefFindAllHeadGDR();
         Period++;  }
   }
}


void FindStartGlobalHeadGDR(void)
{  maxRepBnd=factorEnds[blockSize];
   curMaxPer=maxRepBnd/2-1;
   Period=minPer; 
   if (curMaxPer<maxPer)
   {  while(Period<=curMaxPer)
      {  FindStartHeadGDR();
         Period++;  }
   }
   else
   {  while(Period<=maxPer)
      {  FindStartHeadGDR();
         Period++;  }
   }
}


void FindStartHeadGDR(void)
{  lastPrfEnd=0;
   prfEnd=Period;
   for(numPrfErr=0;numPrfErr<=maxNumErr;numPrfErr++)
   {  while(seq_working_copy[++prfEnd]==seq_working_copy[++lastPrfEnd]);
      if (prfEnd>=maxRepBnd)
         return;
   }
   if (Period<lastPrfEnd)
      NewFirstTypeRepeat(0, prfEnd);
}


void FindFinalGlobalHeadGDR(void)
{  curHeadPstn= wordLength + 1;
   curMaxPer=maxSumFact[curBlock]-2;
   Period=minPer; 
   if (wordLength > 2 * curMaxPer)
   {  if (curMaxPer<maxPer)
      {  while(Period<=curMaxPer)
         {  FindFinalHeadGDR0();
            Period++;  }
      }
      else
      {  while(Period<=maxPer)
         {  FindFinalHeadGDR0();
            Period++;  }
      }
   }
   else
   {  if (curMaxPer<maxPer)
      {  while(Period<=curMaxPer)
         {  FindFinalHeadGDR1();
            Period++;  }
      }
      else
      {  while(Period<=maxPer)
         {  FindFinalHeadGDR1();
            Period++;  }
      }
   }
}


void FindFinalHeadGDR0(void)
{  sufBgn=minSufBgn=curHeadPstn-Period;
   lastSufBgn=curHeadPstn;
   for(numSufErr=0;numSufErr<=maxNumErr;numSufErr++)
      while(seq_working_copy[--sufBgn]==seq_working_copy[--lastSufBgn]);
   if (lastSufBgn<minSufBgn)
      NewFirstTypeRepeat(sufBgn, curHeadPstn);
}

   
void FindFinalHeadGDR1(void)
{  sufBgn=minSufBgn=curHeadPstn-Period;
   lastSufBgn=curHeadPstn;
   for(numSufErr=0;numSufErr<=maxNumErr;numSufErr++)
   {  while(seq_working_copy[--sufBgn]==seq_working_copy[--lastSufBgn]);
      if (sufBgn<=0)
      {  printf("The input sequence is a globally-defined repeat ");
         printf("of period %d with %d errors!!!\n", Period, numSufErr);
         NewFirstTypeRepeat(0, curHeadPstn);
         return;
      }
   }
   if (lastSufBgn<minSufBgn)
      NewFirstTypeRepeat(sufBgn, curHeadPstn);
}


void FindGlobalHeadGDR0(void)
{  int maxLeftRoot;
   maxLeftRoot=maxSumFact[curBlock]-2;
   nextDist=(maxRepBnd-curHeadPstn)-2;
   hlfNxtDist=nextDist/2;
   Period=minPer; 
   if (hlfNxtDist<maxPer)
   {  while(Period<=hlfNxtDist)
      {  FreeFindAllGlobalHeadGDR();
         Period++;  }
   }
   else
   {  while(Period<=maxPer)
      {  FreeFindAllGlobalHeadGDR();
         Period++;  }
      return;
   }
   if (nextDist<maxPer)
   {  while(Period<=nextDist)
      {  FindAllGlobalHeadGDR();
         Period++;  }
   }
   else
   {  while(Period<=maxPer)
      {  FindAllGlobalHeadGDR();
         Period++;  }
      return;
   }
   if (maxLeftRoot<maxPer)
   {  while(Period<=maxLeftRoot)
      {  FindOnlyLeftGlobalHeadGDR();
         Period++;  } 
   }
   else
   {  while(Period<=maxPer)
      {  FindOnlyLeftGlobalHeadGDR();
         Period++;  } 
   }
}


void FindAllGlobalHeadGDR(void)
{  maxPrfEnd=maxRepBnd-(Period+1);
   subPrfEnd=subRepEnd;
   *subPrfEnd=sufBgn=curHeadPstn-Period;
   prfEnd=lastSufBgn=curHeadPstn;
   while(seq_working_copy[prfEnd]==seq_working_copy[*subPrfEnd])
   {  prfEnd++;
      (*subPrfEnd)++;
   }
   if (prfEnd>=maxRepBnd)
      return;
   while(seq_working_copy[--sufBgn]==seq_working_copy[--lastSufBgn]);
   for(numPrfErr=maxNumErr;numPrfErr>0;numPrfErr--)
   {  *(subPrfEnd+1)=(*subPrfEnd);
      subPrfEnd++;
      while(seq_working_copy[++prfEnd]==seq_working_copy[++(*subPrfEnd)]);
      if (prfEnd>=maxRepBnd)
      {  EasyFindAllGlobalHeadGDR();
	 return;
      }
   }
   if (prfEnd>=maxPrfEnd)
   {  if (*subPrfEnd>lastSufBgn)
         NewFirstTypeRepeat(sufBgn, prfEnd);
      FastFindAllGlobalHeadGDR();
   }
   else
      HardFindAllGlobalHeadGDR();
}


void FreeFindAllGlobalHeadGDR(void)
{  maxPrfEnd=curHeadPstn+Period;
   subPrfEnd=subRepEnd;
   *subPrfEnd=sufBgn=curHeadPstn-Period;
   prfEnd=lastSufBgn=curHeadPstn;
   while(seq_working_copy[prfEnd]==seq_working_copy[*subPrfEnd])
   {  prfEnd++;
      (*subPrfEnd)++;
   }
   if (prfEnd>=maxRepBnd)
      return;
   while(seq_working_copy[--sufBgn]==seq_working_copy[--lastSufBgn]);
   for(numPrfErr=maxNumErr;numPrfErr>0;numPrfErr--)
   {  *(subPrfEnd+1)=(*subPrfEnd);
      subPrfEnd++;
      while(seq_working_copy[++prfEnd]==seq_working_copy[++(*subPrfEnd)]);
      if (prfEnd>=maxRepBnd)
      {  EasyFindAllGlobalHeadGDR();
	 return;
      }
   }
   if (prfEnd>=maxPrfEnd)
   {  NewFirstTypeRepeat(sufBgn, prfEnd);
      FastFreeFindAllGlobalHeadGDR();
   }
   else
      HardFindAllGlobalHeadGDR();
}


void FastFindAllGlobalHeadGDR(void)
{  int maxSubPrfEnd;
   for(numSufErr=1;numSufErr<=maxNumErr;numSufErr++)
   {  while(seq_working_copy[--sufBgn]==seq_working_copy[--lastSufBgn]);
      if (*(--subPrfEnd)>lastSufBgn)
         NewFirstTypeRepeat(sufBgn, (*subPrfEnd)+Period);
   }
   maxSubPrfEnd=maxPrfEnd-Period;
   subPrfEnd=subRepEnd;
   lastPrfEnd=subRepEnd[maxNumErr];
   while(*subPrfEnd<maxSubPrfEnd)
   {  while(seq_working_copy[++prfEnd]==seq_working_copy[++lastPrfEnd]);
      if (lastPrfEnd-(*subPrfEnd)>Period)
      {  if (prfEnd<maxRepBnd)
            NewFirstTypeRepeat(*subPrfEnd, prfEnd);
         else
            return;
      }
      subPrfEnd++;
   }
}


void EasyFindAllGlobalHeadGDR(void)
{  numSufErr=1;
   while(numSufErr<numPrfErr)
   {  while(seq_working_copy[--sufBgn]==seq_working_copy[--lastSufBgn]);
      numSufErr++;
   }
   do
   {  while(seq_working_copy[--sufBgn]==seq_working_copy[--lastSufBgn]);
      if (*(--subPrfEnd)>lastSufBgn)
         NewFirstTypeRepeat(sufBgn, (*subPrfEnd)+Period);
   }  while(++numSufErr<=maxNumErr);
}


void FastFreeFindAllGlobalHeadGDR(void)
{  for(numSufErr=1;numSufErr<=maxNumErr;numSufErr++)
   {  while(seq_working_copy[--sufBgn]==seq_working_copy[--lastSufBgn]);
      if (*(--subPrfEnd)>lastSufBgn)
         NewFirstTypeRepeat(sufBgn, (*subPrfEnd)+Period);
   }
   subPrfEnd=subRepEnd;
   lastPrfEnd=subRepEnd[maxNumErr];
   while(*subPrfEnd<curHeadPstn)
   {  while(seq_working_copy[++prfEnd]==seq_working_copy[++lastPrfEnd]);
      if (lastPrfEnd-(*subPrfEnd)>Period)
      {  if (prfEnd<maxRepBnd)
            NewFirstTypeRepeat(*subPrfEnd, prfEnd);
         else
            break;
      }
      subPrfEnd++;
   }
   if (*subPrfEnd==curHeadPstn)
   {  while(seq_working_copy[++prfEnd]==seq_working_copy[++lastPrfEnd]);
      if ( (lastPrfEnd>maxPrfEnd) && (prfEnd<maxRepBnd) )
         NewFirstTypeRepeat(curHeadPstn, prfEnd);
   }
}      


void HardFindAllGlobalHeadGDR(void)
{  if (*subPrfEnd>lastSufBgn)
      NewFirstTypeRepeat(sufBgn, prfEnd);
   for(numSufErr=1;numSufErr<=maxNumErr;numSufErr++)
   {  while(seq_working_copy[--sufBgn]==seq_working_copy[--lastSufBgn]);
      if (*(--subPrfEnd)>lastSufBgn)
         NewFirstTypeRepeat(sufBgn, (*subPrfEnd)+Period);
   }
   subPrfEnd=subRepEnd;
   *subPrfEnd=sufBgn=curHeadPstn;
   prfEnd=lastSufBgn=curHeadPstn+Period;
   while(seq_working_copy[++prfEnd]==seq_working_copy[++(*subPrfEnd)]);
   if (prfEnd>=maxRepBnd)
      return;
   while(seq_working_copy[sufBgn]==seq_working_copy[lastSufBgn])
   {  sufBgn--;
      lastSufBgn--;
   } 
   for(numPrfErr=maxNumErr;numPrfErr>0;numPrfErr--)
   {  *(subPrfEnd+1)=(*subPrfEnd);
      subPrfEnd++;
      while(seq_working_copy[++prfEnd]==seq_working_copy[++(*subPrfEnd)]); 
      if (prfEnd>=maxRepBnd)
      {  EasyFindRightHeadRuns();
	 return;
      }
   }
   if (*subPrfEnd>lastSufBgn)
      NewFirstTypeRepeat(sufBgn, prfEnd);
   for(numSufErr=1;numSufErr<=maxNumErr;numSufErr++)
   {  while(seq_working_copy[--sufBgn]==seq_working_copy[--lastSufBgn]);
      if (*(--subPrfEnd)>lastSufBgn)
         NewFirstTypeRepeat(sufBgn, (*subPrfEnd)+Period);
   }
}


void FindOnlyLeftGlobalHeadGDR(void)
{  subPrfEnd=subRepEnd;
   *subPrfEnd=sufBgn=curHeadPstn-Period;
   prfEnd=lastSufBgn=curHeadPstn;
   while(seq_working_copy[prfEnd]==seq_working_copy[*subPrfEnd])
   {  prfEnd++;
      (*subPrfEnd)++;
   }
   if (prfEnd>=maxRepBnd)
      return;
   while(seq_working_copy[--sufBgn]==seq_working_copy[--lastSufBgn]);
   for(numPrfErr=maxNumErr;numPrfErr>0;numPrfErr--)
   {  *(subPrfEnd+1)=(*subPrfEnd);
      subPrfEnd++;
      while(seq_working_copy[++prfEnd]==seq_working_copy[++(*subPrfEnd)]);
      if (prfEnd>=maxRepBnd)
      {  EasyFindAllGlobalHeadGDR();
	 return;
      }
   }
   if (*subPrfEnd>lastSufBgn)
      NewFirstTypeRepeat(sufBgn, prfEnd);
   for(numSufErr=1;numSufErr<=maxNumErr;numSufErr++)
   {  while(seq_working_copy[--sufBgn]==seq_working_copy[--lastSufBgn]);
      if (*(--subPrfEnd)>lastSufBgn)
         NewFirstTypeRepeat(sufBgn, (*subPrfEnd)+Period);
   }
}


void ComputeRadFact(void)
{  int rad, fact;
   radFact=(int *)calloc(blockSize, sizeof(int));
   for(rad=1;rad<=blockSize/2;rad*=2)
   {  for(fact=rad;fact<blockSize;fact+=2*rad)
         radFact[fact]=rad;
   }
}

