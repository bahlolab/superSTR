#include <stdlib.h>
#include "defs.h"

void FindAllHeadGDR(void);
void RightFreeFindAllHeadGDR(void);
void LeftFreeFindAllHeadGDR(void);
void FreeFindAllHeadGDR(void);
void QuiteEasyFreeFindAllHeadGDR(void);
void QuiteEasyFindAllHeadGDR(void);
void QuiteEasyFreeFindRightHeadGDR(void);
void QuiteEasyFindRightHeadGDR(void);
void EasyFindAllHeadGDR(void);
void FastFindAllHeadGDR(void);
void FastFreeFindAllHeadGDR(void);
void HardFindAllHeadGDR(void);
void EasyFindRightHeadRuns(void);
void BriefFindAllHeadGDR(void);
void FindOnlyLeftHeadGDR(void);
void FindOnlyRightHeadGDR(void);
void EasyFindOnlyRightHeadGDR(void);
void NewFirstTypeRepeat(int newRepBgn, int newRepEnd);
char *seq_working_copy;
int curHeadPstn, actHeadPstn, minRepBnd, maxRepBnd, maxPrfEnd, minSufBgn;
int *subRepEnd, *subRepBgn, *subPrfEnd, *subSufBgn;
int prfEnd, lastPrfEnd, sufBgn, lastSufBgn, numPrfErr, numSufErr;

extern int Period, maxNumErr;
extern listreps *presortLists;


void FindAllHeadGDR(void)
{  maxPrfEnd=maxRepBnd-(Period+1);
   subPrfEnd=subRepEnd;
   *subPrfEnd=sufBgn=minRepBnd+(Period+1);
   prfEnd=lastSufBgn=actHeadPstn=sufBgn+Period;
   while(seq_working_copy[++prfEnd]==seq_working_copy[++(*subPrfEnd)]);
   if (prfEnd>=maxRepBnd)
      return;
   while(seq_working_copy[sufBgn]==seq_working_copy[lastSufBgn])
   {  sufBgn--;
      lastSufBgn--;
   }
   if (sufBgn<=minRepBnd)
   {  if (prfEnd<maxPrfEnd)
         QuiteEasyFindAllHeadGDR();
      return;
   }
   for(numPrfErr=maxNumErr;numPrfErr>0;numPrfErr--)
   {  *(subPrfEnd+1)=(*subPrfEnd);
      subPrfEnd++;
      while(seq_working_copy[++prfEnd]==seq_working_copy[++(*subPrfEnd)]);
      if (prfEnd>=maxRepBnd)
      {  EasyFindAllHeadGDR();
	 return;
      }
   }
   if (prfEnd>=maxPrfEnd)
   {  if (*subPrfEnd>lastSufBgn)
         NewFirstTypeRepeat(sufBgn, prfEnd);
      FastFindAllHeadGDR();
   }
   else
      HardFindAllHeadGDR();
}


void RightFreeFindAllHeadGDR(void)
{  maxPrfEnd=curHeadPstn+Period;
   subPrfEnd=subRepEnd;
   *subPrfEnd=sufBgn=minRepBnd+(Period+1);
   prfEnd=lastSufBgn=sufBgn+Period;
   while(seq_working_copy[++prfEnd]==seq_working_copy[++(*subPrfEnd)]);
   if (prfEnd>=maxRepBnd)
      return;
   while(seq_working_copy[sufBgn]==seq_working_copy[lastSufBgn])
   {  sufBgn--;
      lastSufBgn--;
   }
   if (sufBgn<=minRepBnd)
   {  if (prfEnd<=maxPrfEnd)
         QuiteEasyFreeFindAllHeadGDR();
      return;
   }
   for(numPrfErr=maxNumErr;numPrfErr>0;numPrfErr--)
   {  *(subPrfEnd+1)=(*subPrfEnd);
      subPrfEnd++;
      while(seq_working_copy[++prfEnd]==seq_working_copy[++(*subPrfEnd)]);
      if (prfEnd>=maxRepBnd)
      {  EasyFindAllHeadGDR();
	 return;
      }
   }
   if (prfEnd>=maxPrfEnd)
   {  if (*subPrfEnd>lastSufBgn)
         NewFirstTypeRepeat(sufBgn, prfEnd);
      FastFreeFindAllHeadGDR();
   }
   else
      HardFindAllHeadGDR();
}


void LeftFreeFindAllHeadGDR(void)
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
   if (sufBgn<=minRepBnd)
   {  if (prfEnd<maxPrfEnd)
         QuiteEasyFindAllHeadGDR();
      return;
   }
   for(numPrfErr=maxNumErr;numPrfErr>0;numPrfErr--)
   {  *(subPrfEnd+1)=(*subPrfEnd);
      subPrfEnd++;
      while(seq_working_copy[++prfEnd]==seq_working_copy[++(*subPrfEnd)]);
      if (prfEnd>=maxRepBnd)
      {  EasyFindAllHeadGDR();
	 return;
      }
   }
   if (prfEnd>=maxPrfEnd)
   {  if (*subPrfEnd>lastSufBgn)
         NewFirstTypeRepeat(sufBgn, prfEnd);
      FastFindAllHeadGDR();
   }
   else
      HardFindAllHeadGDR();
}

   
void FreeFindAllHeadGDR(void)
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
   if (sufBgn<=minRepBnd)
   {  if (prfEnd<=maxPrfEnd)
         QuiteEasyFreeFindAllHeadGDR();
      return;
   }
   for(numPrfErr=maxNumErr;numPrfErr>0;numPrfErr--)
   {  *(subPrfEnd+1)=(*subPrfEnd);
      subPrfEnd++;
      while(seq_working_copy[++prfEnd]==seq_working_copy[++(*subPrfEnd)]);
      if (prfEnd>=maxRepBnd)
      {  EasyFindAllHeadGDR();
	 return;
      }
   }
   if (prfEnd>=maxPrfEnd)
   {  NewFirstTypeRepeat(sufBgn, prfEnd);
      FastFreeFindAllHeadGDR();
   }
   else
      HardFindAllHeadGDR();
}


void QuiteEasyFreeFindAllHeadGDR(void)
{  minSufBgn=prfEnd;
   subSufBgn=subRepBgn;
   sufBgn=lastPrfEnd=curHeadPstn;
   *subSufBgn=prfEnd=curHeadPstn+Period;
   while(seq_working_copy[++prfEnd]==seq_working_copy[++lastPrfEnd]);
   if (prfEnd>=maxRepBnd)
      return;
   while(seq_working_copy[sufBgn]==seq_working_copy[*subSufBgn])
   {  sufBgn--;
      (*subSufBgn)--;
   }
   if (*subSufBgn==minSufBgn)
   {  for(numPrfErr=maxNumErr;numPrfErr>0;numPrfErr--)
      {  while(seq_working_copy[++prfEnd]==seq_working_copy[++lastPrfEnd]);
         if (prfEnd>=maxRepBnd)
            return;
      }
      if (lastPrfEnd>minSufBgn)
         NewFirstTypeRepeat(sufBgn, prfEnd);
      return;
   }
   for(numSufErr=maxNumErr;numSufErr>0;numSufErr--)
   {  *(subSufBgn+1)=*subSufBgn;
      subSufBgn++;
      while(seq_working_copy[--sufBgn]==seq_working_copy[--(*subSufBgn)]);
      if (*subSufBgn==minSufBgn)
      {  QuiteEasyFreeFindRightHeadGDR();
         return;
      }
   }
   if (*subSufBgn<lastPrfEnd)
      NewFirstTypeRepeat(sufBgn, prfEnd);
   for(numPrfErr=1;numPrfErr<=maxNumErr;numPrfErr++)
   {  while(seq_working_copy[++prfEnd]==seq_working_copy[++lastPrfEnd]);
      if (*(--subSufBgn)<lastPrfEnd)
      {  if (prfEnd<maxRepBnd)
            NewFirstTypeRepeat((*subSufBgn)-Period, prfEnd);
         else
            return;
      }
   }
}


void QuiteEasyFindAllHeadGDR(void)
{  *subPrfEnd=sufBgn=minSufBgn=prfEnd;
   prfEnd=lastSufBgn=minSufBgn+Period;
   while(seq_working_copy[++prfEnd]==seq_working_copy[++(*subPrfEnd)]);
   if (prfEnd>=maxRepBnd)
      return;
   while(seq_working_copy[sufBgn]==seq_working_copy[lastSufBgn])
   {  sufBgn--;
      lastSufBgn--;
   } 
   if (lastSufBgn==minSufBgn)
   {  lastPrfEnd=(*subPrfEnd);
      for(numPrfErr=maxNumErr;numPrfErr>0;numPrfErr--)
      {  while(seq_working_copy[++prfEnd]==seq_working_copy[++lastPrfEnd]);
         if (prfEnd>=maxRepBnd)
            return;
      }
      NewFirstTypeRepeat(sufBgn, prfEnd);
      return;
   } 
   for(numPrfErr=maxNumErr;numPrfErr>0;numPrfErr--)
   {  *(subPrfEnd+1)=(*subPrfEnd);
      subPrfEnd++;
      while(seq_working_copy[++prfEnd]==seq_working_copy[++(*subPrfEnd)]); 
      if (prfEnd>=maxRepBnd)
      {  QuiteEasyFindRightHeadGDR();
         return;
      }
   }
   if (*subPrfEnd>lastSufBgn)
      NewFirstTypeRepeat(sufBgn, prfEnd);
   for(numSufErr=1;numSufErr<=maxNumErr;numSufErr++)
   {  while(seq_working_copy[--sufBgn]==seq_working_copy[--lastSufBgn]);
      if (*(--subPrfEnd)>lastSufBgn)
      {  NewFirstTypeRepeat(sufBgn, *subPrfEnd+Period);
         if (lastSufBgn==minSufBgn)
            return;
      }
   }
}


void QuiteEasyFreeFindRightHeadGDR(void)
{  numPrfErr=1;
   while(numPrfErr<numSufErr)
   {  while(seq_working_copy[++prfEnd]==seq_working_copy[++lastPrfEnd]);
      if (prfEnd>=maxRepBnd)
         return;
      numPrfErr++;
   }
   if (minSufBgn<lastPrfEnd)
      NewFirstTypeRepeat(sufBgn, prfEnd);
   do
   {  while(seq_working_copy[++prfEnd]==seq_working_copy[++lastPrfEnd]);
      if (*(--subSufBgn)<lastPrfEnd)
      {  if (prfEnd<maxRepBnd)
            NewFirstTypeRepeat((*subSufBgn)-Period, prfEnd);
         else
            return;
      }
   }  while(++numPrfErr<=maxNumErr);
}


void QuiteEasyFindRightHeadGDR(void)
{  numSufErr=1;
   while(numSufErr<numPrfErr)
   {  while(seq_working_copy[--sufBgn]==seq_working_copy[--lastSufBgn]);
      if (lastSufBgn==minSufBgn)
         return;
      numSufErr++;
   }
   do
   {  while(seq_working_copy[--sufBgn]==seq_working_copy[--lastSufBgn]);
      if (*(--subPrfEnd)>lastSufBgn)
      {  NewFirstTypeRepeat(sufBgn, *subPrfEnd+Period);
         if (lastSufBgn==minSufBgn)
            return;
      }
   }  while(++numSufErr<=maxNumErr);
}


void EasyFindAllHeadGDR(void)
{  numSufErr=1;
   while(numSufErr<numPrfErr)
   {  while(seq_working_copy[--sufBgn]==seq_working_copy[--lastSufBgn]);
      if (sufBgn<=minRepBnd)
         return;
      numSufErr++;
   }
   do
   {  while(seq_working_copy[--sufBgn]==seq_working_copy[--lastSufBgn]);
      if (*(--subPrfEnd)>lastSufBgn)
      {  if (sufBgn>minRepBnd)
            NewFirstTypeRepeat(sufBgn, *subPrfEnd+Period);
         else
            return;
      }
   }  while(++numSufErr<=maxNumErr);
}


void FastFindAllHeadGDR(void)
{  int maxSubPrfEnd;
   for(numSufErr=1;numSufErr<=maxNumErr;numSufErr++)
   {  while(seq_working_copy[--sufBgn]==seq_working_copy[--lastSufBgn]);
      if (*(--subPrfEnd)>lastSufBgn)
      {  if (sufBgn>minRepBnd)
            NewFirstTypeRepeat(sufBgn, *subPrfEnd+Period);
         else
            break;
      }
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


void FastFreeFindAllHeadGDR(void)
{  for(numSufErr=1;numSufErr<=maxNumErr;numSufErr++)
   {  while(seq_working_copy[--sufBgn]==seq_working_copy[--lastSufBgn]);
      if (*(--subPrfEnd)>lastSufBgn)
      {  if (sufBgn>minRepBnd)
            NewFirstTypeRepeat(sufBgn, *subPrfEnd+Period);
         else
            break;
      }
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


void HardFindAllHeadGDR(void)
{  if (*subPrfEnd>lastSufBgn)
      NewFirstTypeRepeat(sufBgn, prfEnd);
   for(numSufErr=1;numSufErr<=maxNumErr;numSufErr++)
   {  while(seq_working_copy[--sufBgn]==seq_working_copy[--lastSufBgn]);
      if (*(--subPrfEnd)>lastSufBgn)
      {  if (sufBgn>minRepBnd)
            NewFirstTypeRepeat(sufBgn, *subPrfEnd+Period);
         else
            break;
      }
   }
   subPrfEnd=subRepEnd;
   *subPrfEnd=sufBgn=actHeadPstn;
   prfEnd=lastSufBgn=actHeadPstn+Period;
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
         NewFirstTypeRepeat(sufBgn, *subPrfEnd+Period);
   }
}


void EasyFindRightHeadRuns(void)
{  numSufErr=1;
   while(numSufErr<numPrfErr)
   {  while(seq_working_copy[--sufBgn]==seq_working_copy[--lastSufBgn]);
      numSufErr++;
   }
   do
   {  while(seq_working_copy[--sufBgn]==seq_working_copy[--lastSufBgn]);
      if (*(--subPrfEnd)>lastSufBgn)
         NewFirstTypeRepeat(sufBgn, *subPrfEnd+Period);
   }  while(++numSufErr<=maxNumErr);
}


void BriefFindAllHeadGDR(void)
{  subPrfEnd=subRepEnd;
   *subPrfEnd=sufBgn=minRepBnd+(Period+1);
   prfEnd=lastSufBgn=sufBgn+Period;
   while(seq_working_copy[++prfEnd]==seq_working_copy[++(*subPrfEnd)]);
   if (prfEnd>=maxRepBnd)
      return;
   while(seq_working_copy[sufBgn]==seq_working_copy[lastSufBgn])
   {  sufBgn--;
      lastSufBgn--;
   }
   if (sufBgn<=minRepBnd)
      return;
   for(numPrfErr=maxNumErr;numPrfErr>0;numPrfErr--)
   {  *(subPrfEnd+1)=(*subPrfEnd);
      subPrfEnd++;
      while(seq_working_copy[++prfEnd]==seq_working_copy[++(*subPrfEnd)]);
      if (prfEnd>=maxRepBnd)
      {  EasyFindAllHeadGDR();
	 return;
      }
   }
   if (*subPrfEnd>lastSufBgn)
      NewFirstTypeRepeat(sufBgn, prfEnd);
   for(numSufErr=1;numSufErr<=maxNumErr;numSufErr++)
   {  while(seq_working_copy[--sufBgn]==seq_working_copy[--lastSufBgn]);
      if (*(--subPrfEnd)>lastSufBgn)
      {  if (sufBgn>minRepBnd)
            NewFirstTypeRepeat(sufBgn, *subPrfEnd+Period);
         else
            return;
      }
   }
}


void FindOnlyLeftHeadGDR(void)
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
   if (sufBgn<=minRepBnd)
      return;
   for(numPrfErr=maxNumErr;numPrfErr>0;numPrfErr--)
   {  *(subPrfEnd+1)=(*subPrfEnd);
      subPrfEnd++;
      while(seq_working_copy[++prfEnd]==seq_working_copy[++(*subPrfEnd)]);
      if (prfEnd>=maxRepBnd)
      {  EasyFindAllHeadGDR();
	 return;
      }
   }
   if (*subPrfEnd>lastSufBgn)
      NewFirstTypeRepeat(sufBgn, prfEnd);
   for(numSufErr=1;numSufErr<=maxNumErr;numSufErr++)
   {  while(seq_working_copy[--sufBgn]==seq_working_copy[--lastSufBgn]);
      if (*(--subPrfEnd)>lastSufBgn)
      {  if (sufBgn>minRepBnd)
            NewFirstTypeRepeat(sufBgn, *subPrfEnd+Period);
         else
            return;
      }
   }
}


void FindOnlyRightHeadGDR(void)
{  subSufBgn=subRepBgn;
   sufBgn=lastPrfEnd=curHeadPstn;
   *subSufBgn=prfEnd=curHeadPstn+Period;
   while(seq_working_copy[sufBgn]==seq_working_copy[*subSufBgn])
   {  sufBgn--;
      (*subSufBgn)--;
   }
   if (sufBgn<=minRepBnd)
      return;
   while(seq_working_copy[++prfEnd]==seq_working_copy[++lastPrfEnd]);
   if (prfEnd>=maxRepBnd)
      return;
   for(numSufErr=maxNumErr;numSufErr>0;numSufErr--)
   {  *(subSufBgn+1)=*subSufBgn;
      subSufBgn++;
      while(seq_working_copy[--sufBgn]==seq_working_copy[--(*subSufBgn)]);
      if (sufBgn<=minRepBnd)
      {  EasyFindOnlyRightHeadGDR();
         return;
      }
   }
   if (*subSufBgn<lastPrfEnd)
      NewFirstTypeRepeat(sufBgn, prfEnd);
   for(numPrfErr=1;numPrfErr<=maxNumErr;numPrfErr++)
   {  while(seq_working_copy[++prfEnd]==seq_working_copy[++lastPrfEnd]);
      if (*(--subSufBgn)<lastPrfEnd)
      {  if (prfEnd<maxRepBnd)
            NewFirstTypeRepeat((*subSufBgn)-Period, prfEnd);
         else
            return;
      }
   }
}


void EasyFindOnlyRightHeadGDR(void)
{  numPrfErr=1;
   while(numPrfErr<numSufErr)
   {  while(seq_working_copy[++prfEnd]==seq_working_copy[++lastPrfEnd]);
      if (prfEnd>=maxRepBnd)
         return;
      numPrfErr++;
   }
   do
   {  while(seq_working_copy[++prfEnd]==seq_working_copy[++lastPrfEnd]);
      if (*(--subSufBgn)<lastPrfEnd)
      {  if (prfEnd<maxRepBnd)
            NewFirstTypeRepeat((*subSufBgn)-Period, prfEnd);
         else
            return;
      }
   }  while(++numPrfErr<=maxNumErr);
}


void NewFirstTypeRepeat(int newRepBgn, int newRepEnd)
{  listreps newListRep;
   repetition *newRep;
   if ( (newListRep=(listreps)malloc(sizeof(struct s_listreps)))==NULL )
   {  printf("The output is too big.\n");
      printf("Narrow the Period bounds.\n");
      exit(14);
   }
   newRep=&(newListRep->rep);
   newRep->initpos=newRepBgn;
   newRep->endpos=newRepEnd;
   newRep->period=Period;
   newListRep->next=presortLists[newRepEnd];
   presortLists[newRepEnd]=newListRep;
}
