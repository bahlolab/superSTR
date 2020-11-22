/* mainRepeats.c : computing repetitions of the first type (exact case) */

#include <stdlib.h>
#include <stdio.h>
#include "defs.h"

#define MINDIFLENPER       9
extern int allowsmall;

extern int *computeLIx(char *pattern, int patLen, int dir, int maxPerBnd);

extern int *computeLOx(char *pattern, int patLen, int dir, int maxPerBnd);

/* extern void NewFirstTypeRepeat(int newRepBgn, int newRepEnd); */
extern listreps *presortLists;
/* extern listreps *sortedLists; */
extern int maxPer;
extern char *seq_original;
/* extern float minExp; */
extern limits lim;
extern int step, LastWindow;

extern FILE *output_file;                 // output file (for xml output)
extern char *output_file_name;
extern int xmloutput;


void NewFirstTypeExactRepeat(int prRepBgn, int prRepEnd);

void startrepets(int curFacBgn, int nxtFacBgn);

void finalRepeats(int prvFacBgn, int curFacBgn);

int Period;                            // declare it here! (rather than in searchforHeadGDR.c)
int prvFacLen, maxLftLen, curFacLen, maxLftRoot;


void mainRepeats(int prvFacBgn, int curFacBgn, int nxtFacBgn)
/* searching for exact reps of the first type crossing the
current border between factors in Crochemore's s-factorisation,
and located within the current window;
called by maxreps */
{
    int rightLen;
    int *LPu, *LSt, *LPt, *LSu;

    curFacLen = nxtFacBgn - curFacBgn;
    prvFacLen = curFacBgn - prvFacBgn;
    maxLftRoot = curFacLen + prvFacLen;
    if ((maxLftLen = maxLftRoot + prvFacLen) > curFacBgn) {
        maxLftLen = curFacBgn;
        if (maxLftRoot > curFacBgn) {
            maxLftRoot = curFacBgn + 1;
            if (maxLftRoot <= curFacLen) {
                startrepets(curFacBgn, nxtFacBgn);
                return;
            }
        }

    }

    /* computing Main's extension functions */
    LPu = computeLIx(seq_original + curFacBgn, curFacLen, 1, curFacLen);
    LSt = computeLIx(seq_original + curFacBgn - 1, maxLftLen, -1, maxLftRoot);

    LPt = computeLOx(seq_original + curFacBgn, curFacLen, 1, maxLftRoot);
    LSu = computeLOx(seq_original + curFacBgn - 1, maxLftLen, -1, curFacLen);

    Period = lim.min_period;
    while (Period < curFacLen && Period <= maxPer)    // LIMIT
    {
        /* reps with a full period in previous factors */

        if (LPt[Period] + LSt[Period] >= Period) {
            if ((rightLen = LPt[Period]) >= Period) {
                if ((rightLen < curFacLen) ||
                    (seq_original[nxtFacBgn] != seq_original[nxtFacBgn - Period]))
                    NewFirstTypeExactRepeat(curFacBgn - (LSt[Period] + Period),
                                            curFacBgn + rightLen);
                Period++;
                continue;
            } else {
                if ((rightLen > 0) ||
                    (LSt[Period] + Period < prvFacLen))
                    NewFirstTypeExactRepeat(curFacBgn - (LSt[Period] + Period),
                                            curFacBgn + rightLen);
            }
        }

        /* m-repets with a full period in current factor */

        if (LSu[Period] + LPu[Period] >= Period) {
            if ((LPu[Period] + Period < curFacLen) ||
                (seq_original[nxtFacBgn] != seq_original[nxtFacBgn - Period]))
                NewFirstTypeExactRepeat(curFacBgn - LSu[Period],
                                        curFacBgn + LPu[Period] + Period);
        }
        Period++;
    }

    while (Period < maxLftRoot && Period <= maxPer)    // LIMIT
    {
        /* reps with a full period in previous factors only */

        if (LPt[Period] + LSt[Period] >= Period) {
            if ((rightLen = LPt[Period]) < curFacLen) {
                if ((rightLen > 0) ||
                    (LSt[Period] + Period < prvFacLen))
                    NewFirstTypeExactRepeat(curFacBgn - (LSt[Period] + Period),
                                            curFacBgn + rightLen);
            } else {
                if (seq_original[nxtFacBgn] != seq_original[nxtFacBgn - Period])
                    NewFirstTypeExactRepeat(curFacBgn - (LSt[Period] + Period),
                                            curFacBgn + rightLen);
            }
        }
        Period++;
    }
    free(LPu);
    free(LSu);
    free(LPt);
    free(LSt);
}


void NewFirstTypeExactRepeat(int newRepBgn, int newRepEnd)
/* called from mainRepeats, startrepets and finalRepeats;
processing an exact rep of the first type. For a found rep,
creates the corresponding structure 'listreps' and
inserts it into the corresponing list from 'presortLists'
(array of lists of reps of the first type found in the
current window and sorted by end positions) */
{
    listreps newListRep, nextRep;
    repetition *newRep;
    int rsize;

    rsize = newRepEnd - newRepBgn;         // the size of the rep

    if (((rsize - Period >= MINDIFLENPER) || allowsmall) &&
        (step < 0 || newRepBgn < VIRTUAL_STEP)                        // GK
            ) {
        if ((nextRep = presortLists[newRepEnd])) {
            if (nextRep->rep.initpos == newRepBgn) {
                return;
            }
        }
        if ((newListRep = (listreps) malloc(sizeof(struct s_listreps))) == NULL) {
            fprintf(stderr, "The output is too big.\n");
            fprintf(stderr, "Narrow period bounds.\n");
            exit(20);
        }
        newRep = &(newListRep->rep);
        newRep->initpos = newRepBgn;
        newRep->endpos = newRepEnd;
        newRep->period = Period;
        newListRep->next = nextRep;
        presortLists[newRepEnd] = newListRep;
    }
}


void startrepets(int curFacBgn, int nxtFacBgn)
/* called from mainrepets;
similar to 'mainRepeats' for the case when the reps can
contain the beginning of the current window */
{
    int rightLen;
    int *LPu, *LSt, *LPt, *LSu;


    /* computing longest extension functions */
    LPu = computeLIx(seq_original + curFacBgn, curFacLen, 1, curFacLen);
    LSt = computeLIx(seq_original + curFacBgn - 1, maxLftLen, -1, maxLftRoot);

    LPt = computeLOx(seq_original + curFacBgn, curFacLen, 1, maxLftRoot);
    LSu = computeLOx(seq_original + curFacBgn - 1, maxLftLen, -1, curFacLen);

    Period = lim.min_period;
    while (Period < maxLftRoot && Period <= maxPer)    // LIMIT
    {
        /* reps with a full period in previous factors */

        if (LPt[Period] + LSt[Period] >= Period) {
            if ((rightLen = LPt[Period]) >= Period) {
                if ((rightLen < curFacLen) ||
                    (seq_original[nxtFacBgn] != seq_original[nxtFacBgn - Period]))
                    NewFirstTypeExactRepeat(curFacBgn - (LSt[Period] + Period),
                                            curFacBgn + rightLen);
                Period++;
                continue;
            } else {
                if ((rightLen > 0) ||
                    (LSt[Period] + Period < prvFacLen))
                    NewFirstTypeExactRepeat(curFacBgn - (LSt[Period] + Period),
                                            curFacBgn + rightLen);
            }
        }

        /* reps with a full period in current factor */

        if (LSu[Period] + LPu[Period] >= Period) {
            if ((LPu[Period] + Period < curFacLen) ||
                (seq_original[nxtFacBgn] != seq_original[nxtFacBgn - Period]))
                NewFirstTypeExactRepeat(curFacBgn - LSu[Period],
                                        curFacBgn + LPu[Period] + Period);
        }
        Period++;
    }

    while (Period < curFacLen && Period <= maxPer)    // LIMIT
    {
        /* reps with a full period in current factor only */
        if (LSu[Period] + LPu[Period] >= Period) {
            if ((LPu[Period] + Period < curFacLen) ||
                (seq_original[nxtFacBgn] != seq_original[nxtFacBgn - Period]))
                NewFirstTypeExactRepeat(curFacBgn - LSu[Period],
                                        curFacBgn + LPu[Period] + Period);
        }
        Period++;
    }
    free(LPu);
    free(LSu);
    free(LPt);
    free(LSt);
}


void finalRepeats(int prvFacBgn, int curFacBgn)
/* called from maxreps;
searching for exact reps of the first type
containing the end of the current window,
but not the last border between s-factors (?)
(in other words, all exact reps of the first
type in the end of word, which
have not been found by mainRepeats) */
{
    int *LSt;
    prvFacLen = curFacBgn - prvFacBgn;
    maxLftRoot = (1 + prvFacLen) / 2;


    /* computing longest extension functions */
    LSt = computeLIx(seq_original + curFacBgn - 1, prvFacLen, -1, maxLftRoot);

    Period = lim.min_period;
    while (Period < maxLftRoot && Period <= maxPer)    // LIMIT
    {
        /* reps in previous factor completely */
        if (LSt[Period] >= Period) {
            if (LSt[Period] + Period < prvFacLen)
                NewFirstTypeExactRepeat(curFacBgn - (LSt[Period] + Period),
                                        curFacBgn);
        }
        Period++;
    }
    free(LSt);
}
