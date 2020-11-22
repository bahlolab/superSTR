/* Computing exact repetitions */

#include <stdlib.h>
#include <stdio.h>
#include "defs.h"

// From main.c
extern limits lim;
extern int step, LastWindow, maxPer, oldPerFactorBound;;
extern char nextLetterCheck, prevLetterCheck;

// From factorizeforGDR.c
extern int FACTYPE, factorNum, *factorStarts;
extern char *seq_factor_copy;
extern sfactorization FactorizeforGDR(void);
extern void guessDepth(int length);

// From mainRepeats.c
extern void mainRepeats(int prvFacBgn, int curFacBgn, int nxtFacBgn);
extern void finalRepeats(int prvFacBgn, int curFacBgn);

// From kmp.c
extern int *computeLIx(char *pattern, int patLen, int dir, int maxPerBnd);

// From finalstageforGDR.c
extern listreps *presortLists, *sortedLists, *maxIntPtr;
extern void EmergSortFirstTypeRepeats(void);
extern void FindSecondTypeExactRepeats(void);

// From findReps.c
extern listreps *allReps;
extern void sequence_sanitise(int chkLength);

// From printOutput.c
extern void print_score(head_t*, int, int, int , int , float );

// Locally declared variables
char *seq_original, *seq_original1;   // treated window of original sequence
int verbose = 1;         // defined as extern in defs.h
int verbose_long = 1;   // defined as extern in defs.h
int wordLength;    // declare it here!
int tooBigReps;

// Declaration of methods
void Savemreps(head_t*);
void Showmreps(head_t*);
void print_mrepeat(head_t*, int rinitpos, int rendpos, int rperiod);
void ComputePerBound(void);
int check_mrepeat(listreps checkRepeat);

void maxreps(head_t* headnode, char *seq_parameter, int length_parameter) { /* seaching for exact reps in the current window */
    int i;
    // Variable to hold s-factorisation
    sfactorization sfact;
    // Store original sequences in globals:
    seq_original = seq_parameter; // Sequence address
    wordLength = length_parameter;    // Sequence length
    // Rename to better variable
    seq_original1 = seq_original-1;
    // The maximum possible period is the length/2 or the user set limit.
    maxPer = (lim.max_period>length_parameter/2) ? length_parameter/2 :lim.max_period;
    /* Generate s-factorization */
    // Copy the sequence:
    seq_factor_copy = (char *)malloc(length_parameter+2);
    seq_factor_copy[0]='\0';
    // Sanitise the input - check Ns, capitalise everything.
    sequence_sanitise(length_parameter);
    FACTYPE=S_FACT_EXACT_MATCH;
    guessDepth(length_parameter);
    sfact = FactorizeforGDR();
    if (lim.err_number)
        ComputePerBound();
    else
        free(seq_factor_copy);
    // Debug - can print sfactorisation count here: printf("s-factorization : %d factors\n", sfact->nbfactors);
    // Detect repeats of type 1 (mainRepeats.c)
    if(( presortLists=(listreps *) calloc ((2 + wordLength), sizeof(listreps)) ) == NULL ) { // TODO: Exception
        fprintf(stderr,"Not enough memory for saving the repetitions\n");
        exit(20);
    }
    for (i=1; i< sfact -> nbfactors ; i++) {
        // Debug - can print sfact here with printf(" s-fact: [%d, %d]\n", sfact->factor[i], sfact->factor[i+1]);
        mainRepeats(sfact->factor[i - 1], sfact->factor[i], sfact->factor[i + 1]);
    }
    finalRepeats(sfact->factor[i - 1], sfact->factor[i]);
    // Sort repeats
    EmergSortFirstTypeRepeats();
    /* Detect repeats of type 2 */
    FindSecondTypeExactRepeats();
    free (sfact);
    /* Return results */
    if (lim.err_number==0) {
        Showmreps(headnode);
    } else {
        Savemreps(headnode);
    }
}


void Showmreps(head_t* repeat_head)
/* filtering out exact reps which start in the second half of the window (unless it is the last one),
or extend beyond the window; call print_mrepeat on those reps which go through filtering */
{
    int endPstn, initPos, repPer;
    listreps current_rep, *current_pointer, remListRep;
    current_pointer=sortedLists;
    if ( (current_rep=(*current_pointer)) )  // if there are reps starting at the first position treat them separately (because of prevLetterCheck)
        do {
            repPer=current_rep->rep.period;
            if (seq_original1[repPer]!=prevLetterCheck) { // repeat does not extend to the left
                if ((endPstn=current_rep->rep.endpos) < wordLength ) {// it does not touch the end
                    print_mrepeat(repeat_head, 0, endPstn, repPer);
                } else {
                    if (seq_original[wordLength - repPer] != nextLetterCheck) {
                        print_mrepeat(repeat_head, 0, wordLength, repPer);
                    } else {
                        tooBigReps = YES;
                    }
                }
            }
            current_rep=(remListRep=current_rep)->next;
            free(remListRep);
        }
        while(current_rep);

    while(++current_pointer<=maxIntPtr) {
        if ((current_rep=(*current_pointer))) {
            initPos=current_rep->rep.initpos;
            do {
                repPer=current_rep->rep.period;
                if ((endPstn=current_rep->rep.endpos) < wordLength ) // the repeat does not
                    // touch the end of the window
                    print_mrepeat(repeat_head, initPos, endPstn, repPer);
                else {
                    if (seq_original[wordLength - repPer] != nextLetterCheck)
                        print_mrepeat(repeat_head, initPos, wordLength, repPer);
                    else
                        tooBigReps=YES;
                }
                current_rep=(remListRep=current_rep)->next;
                free(remListRep);
            }
            while(current_rep);
        }
    }
    free(sortedLists);
}

void  Savemreps(head_t* dr)
/* filtering out exact reps which start in the second half of the window (unless it is the last one),
or extend beyond the window; save those reps which go through filtering */
{
    listreps current_rep, *current_pointer, next_rep;

    if ((allReps=(listreps *)calloc(oldPerFactorBound, sizeof(listreps))) == NULL )
    {
        printf("Error: Not enough memory for storing exact repeats\n");
        exit(20);
    }

    if (LastWindow==NO)
        maxIntPtr=sortedLists+step-1;                  // output repetition starting in the first half only

    for(current_pointer=maxIntPtr;current_pointer>sortedLists;current_pointer--)
    {
        if ((current_rep=(*current_pointer)))
        { while ( (next_rep=current_rep->next) )
            {  if ( check_mrepeat(current_rep) )
                {  current_rep->next=allReps[current_rep->rep.period];
                    allReps[current_rep->rep.period]=current_rep;
                }
                else
                    free(current_rep);
                current_rep=next_rep;
            }
            if ((current_rep->rep.endpos < wordLength) ||              // the repeat does not
                 (seq_original[wordLength - current_rep->rep.period]    // cross the end of the window
                  != nextLetterCheck) )
            {  if ( check_mrepeat(current_rep) )
                {  current_rep->next=allReps[current_rep->rep.period];
                    allReps[current_rep->rep.period]=current_rep;
                }
                else
                    free(current_rep);
            }
            else
            {  free(current_rep);
                tooBigReps=YES;
            }
        }
    }

    if ( (current_rep=(*current_pointer)) )  // if there are reps starting at the first position
        // treat them separately (because of prevLetterCheck)
    {  while ( (next_rep=current_rep->next) )
        {  if (seq_original1[current_rep->rep.period]!=prevLetterCheck)
                // repeat does not extend to the left
            {  if ( check_mrepeat(current_rep) )
                {  current_rep->next=allReps[current_rep->rep.period];
                    allReps[current_rep->rep.period]=current_rep;
                }
                else
                    free(current_rep);
            }
            else
                free(current_rep);
            current_rep=next_rep;
        }
        if (seq_original1[current_rep->rep.period]!=prevLetterCheck)
            // repeat does not extend to the left
        {  if ((current_rep->rep.endpos < wordLength) ||              // the repeat does not
                (seq_original[wordLength - current_rep->rep.period]    // cross the end of the window
                 != nextLetterCheck) )
            {  if ( check_mrepeat(current_rep) )
                {  current_rep->next=allReps[current_rep->rep.period];
                    allReps[current_rep->rep.period]=current_rep;
                }
                else
                    free(current_rep);
            }
            else
            {  free(current_rep);
                tooBigReps=YES;
            }
        }
        else
            free(current_rep);
    }

    free(sortedLists);
}


int check_mrepeat(listreps checkRep)
/* check if the period of the given exact rep is primitive and pass it further for
   saving in structure "allReps" */
{ int rinitpos, rendpos;
    int rlength, j, *LP;

    rinitpos=checkRep->rep.initpos++;
    rendpos=checkRep->rep.endpos++;
    rlength=rendpos-rinitpos;

    if (lim.min_period>1)
    {
/*   compute the minimal period */
        LP=computeLIx(seq_original+rinitpos, rlength, 1, lim.min_period);
        for (j = 1 ; j<lim.min_period ; j++)
        {
            if ( j+LP[j]>=rlength )
            { free(LP);       // throw away the repetition if its minimal period
                return(NO); }   // is smaller than the user-specified minimal period
        }
        free(LP);
    }

    return(YES);
}


void ComputePerBound() {
    int factorIndex, factorSum;
    if (factorNum < 3) { // TODO: Convert to exception
        printf("Error: number of errors is too big for this sequence\n");
        printf("       (or the sequence has too regular structure)\n");
        printf("       Run with a smaller number of errors\n");
        exit(71);
    }
    oldPerFactorBound=factorStarts[factorIndex=2];
    while(++factorIndex <= factorNum) {
        if ((factorSum= factorStarts[factorIndex] - factorStarts[factorIndex - 2]) > oldPerFactorBound)
            oldPerFactorBound=factorSum;
    }
}


void print_mrepeat(head_t* repeat_head, int rinitpos, int rendpos, int rperiod) {
/* check if the period of the given exact rep is primitive and pass it further for printing
to function 'print_rep' (for exe version) or 'printWeb' (for web version)
(both defined in file printOutput.c) */
    int rlength, j, *LP;
    float rexp;
    /* compute the minimal period */
    rlength=rendpos-rinitpos; // the length of the repeat
    if (lim.min_period>1) {
        LP=computeLIx(seq_original+rinitpos, rlength, 1, lim.min_period);
        for (j = 1 ; j<lim.min_period ; j++) {
            if ( j+LP[j]>=rlength ) {
                free(LP);   // throw away the repetition if its minimal period is smaller than user-specified minimal period
                return;
            }
        }
        free(LP);
    }
    rexp = (float)rlength/rperiod; // the exp of the repeat
    if ((rlength >= lim.min_size) && (rlength <= lim.max_size)    // LIMIT for size
         && (rexp >= lim.min_exponent)) {                 // LIMIT for exponent
        //This part is for the exe version of the program. Comment it out for web version.
        print_score(repeat_head, rinitpos + 1, rendpos, rlength, rperiod, 0.0);
    }
}

