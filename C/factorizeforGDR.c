#include <math.h>
#include "defs.h"

/* #include "factorizeforGDR.h" */

#define SECARRSIZ 10000

extern int wordLength, maxNumErr;

void InitPrTree(void);

void MovePrFacNode(void);

void MakeSufTree(void);

void EnterSecTree(void);

void TreatInitSecNode(void);

void MakeInitSecNode(void);

void GoInitFacNode(void);

void MoveSecFacNode(void);

void TreatOldSecNode(void);

void MakeNewSecNode(void);

void GoOldSecNode(void);

void ComputeFinFactor(void);

void ComputeGlobalFactors(void);   // called only in FndReps.c
void MakePowof4(void);

void guessDepth(int length);   // called in superstr.c FndReps.c
int ComputeSizeFact(void);


int pow4[16], sumpow4[16], prPow4, prSumpow4;
char *seq_factor_copy;
int maxNumErr2, actPerBnd;
int *prefixTree, *lvList, prLeaf;
int (**secTree)[5], numSecNod, numSecArr;
int blockSize, finBlockSize, *maxSumFact, sizeFact;
int *copyEnds, *copyBgns, *factorEnds;
int *factorStarts;    // factorStarts[i] is the last (head) symbol
// (in natural enumeration) of i-th factor
int factorNum, numGlobFact, numFullBlock;
int curPrNod, curFacNod, currentPosition, currentSymbol;
int facSecDep, facSecLeaf, *facSecNod;
int newList[5], newSecLeaf[5], nxtSecLeaf, curSecLeaf;
int prevCopyPstn, minSecDep;
int PRIMEDEPTH, FACTYPE, offset_primedepth;

// TODO: Document this out later
void guessDepth(int length) {
    double logSize;
    int depth;
    logSize = log10((double) (length / 1000.0));
    //depth is approximately linear with logSize. depth = 2x + 2;
    depth = (int) (logSize * 2 + 2);
    if (depth > 10) depth = 10;
    else if (depth < 4) depth = 4;
    PRIMEDEPTH = depth;
    offset_primedepth = PRIMEDEPTH + 1;
}

// Factorise to get globally defined repetition per Kolpakov/Kucherov 2013
sfactorization FactorizeforGDR(void) {
    // Estimate the size of the s-factorisation data structure
    // ********************************************************
    MakePowof4(); // TODO: Can this be split so we're not repeatedly computing powers of 4?
    // Declare and initialise an sfactorisation data structure
    sfactorization fStruct;
    fStruct = (sfactorization) malloc(sizeof(struct s_sfactorization));
    sizeFact = ComputeSizeFact();
    int indx;
    // Try and create the factorisation structure (fStruct), checking for out-of-memory conditions.
    // For LZ factorisation:
    if (FACTYPE == LZ_K_MISMATCH) {
        // Store previous run.
        if ((fStruct->previous = copyEnds = (int *) malloc((sizeFact + 1) * sizeof(int))) == NULL) { // TODO: Exception
            printf("Not enough memory for factorization\n");
            exit(29);
        }
        // Create space for new runs.
        if ((fStruct->factor = factorEnds = (int *) malloc((sizeFact + 1) * sizeof(int))) == NULL) { // TODO: Exception
            printf("Not enough memory for factorization\n");
            exit(30);
        }
        *(factorEnds++) = (*(copyEnds++)) = 0;
        prevCopyPstn = (*(factorEnds++)) = (*(copyEnds++)) = 1;
    // For s-factorisation:
    } else {
        if ((fStruct->previous = copyBgns = (int *) malloc((sizeFact + 1) * sizeof(int))) == NULL) { // TODO: Exception
            printf("Not enough memory for factorization\n");
            exit(29);
        }
        if ((fStruct->factor = factorStarts = (int *) malloc((sizeFact + 1) * sizeof(int))) == NULL) { // TODO: Exception
            printf("Not enough memory for factorization\n");
            exit(30);
        }
        *(factorStarts++) = (*(copyBgns++)) = 0;
        *factorStarts = 1;
    }
    // TODO: Factor num?
    factorNum = 2;
    // Create a prefix tree
    if ((prefixTree = (int *) calloc(sumpow4[PRIMEDEPTH + 2], sizeof(int))) == NULL) { // TODO: Exception
        printf("Not enough memory for factorization\n");
        exit(9);
    }
    // Secondary array
    secTree = (int (**)[5]) calloc(wordLength / SECARRSIZ + 1, sizeof(int (*)[5]));
    if ((secTree[numSecArr = 0] = (int (*)[5]) calloc(SECARRSIZ, sizeof(int[5]))) == NULL) { // TODO: Exception
        printf("Can't allocate memory for first secondary array\n");
        exit(11);
    }
    // Create a level list.
    if ((lvList = (int *) calloc(wordLength + 1, sizeof(int))) == NULL) { // TODO: Exception
        printf("Can't allocate memory for factorization\n");
        exit(10);
    }
    facSecDep = 0;
    curFacNod = 0;
    numSecNod = 0;
    InitPrTree();
    while (currentPosition <= wordLength) {
        MakeSufTree();
        currentPosition++;
    }
    ComputeFinFactor();
    free(lvList);
    for (indx = 0; indx <= numSecArr; indx++)
        free(secTree[indx]);
    free(secTree);
    free(prefixTree);
    fStruct->nbfactors = factorNum;
    if (FACTYPE == LZ_K_MISMATCH) {
        copyEnds = fStruct->previous;
        factorEnds = fStruct->factor;
    } else {
        copyBgns = fStruct->previous;
        factorStarts = fStruct->factor;
    }
    return (fStruct);
}


void InitPrTree(void) {
    int sufNode, suffixLength, maxPosition;
    prefixTree[0] = -1;
    prefixTree[curPrNod = (int) seq_factor_copy[1]] = 1;
    maxPosition = (PRIMEDEPTH < wordLength) ? PRIMEDEPTH : wordLength;
    for (currentPosition = 2; currentPosition <= maxPosition; currentPosition++) {
        currentSymbol = (int) seq_factor_copy[currentPosition];
        MovePrFacNode();
        curPrNod = 4 * curPrNod + currentSymbol;
        prefixTree[curPrNod] = currentPosition;
        suffixLength = currentPosition - 1;
        sufNode = (curPrNod - sumpow4[suffixLength]) % pow4[suffixLength] + sumpow4[suffixLength];
        while (prefixTree[sufNode] == 0) {
            prefixTree[sufNode] = currentPosition;
            suffixLength--;
            sufNode = (sufNode - sumpow4[suffixLength]) % pow4[suffixLength] + sumpow4[suffixLength];
        }
    }
}


void MovePrFacNode(void) {
    int curCopyPstn, curFacBgn;
    curFacNod = 4 * curFacNod + currentSymbol;
    if ((curCopyPstn = prefixTree[curFacNod]))
        prevCopyPstn = curCopyPstn;
    else {
        if (FACTYPE == LZ_K_MISMATCH) {
            *(copyEnds++) = prevCopyPstn + 1;
            *(factorEnds++) = prevCopyPstn = currentPosition;
            factorNum++;
            curFacNod = 0;
        } else {
            if ((curFacBgn = currentPosition - 1) > (*factorStarts)) {
                *(copyBgns++) = prevCopyPstn + (*factorStarts) - curFacBgn;
                *(++factorStarts) = curFacBgn;
                factorNum++;
                if ((prevCopyPstn = prefixTree[curFacNod = currentSymbol]) == 0) {
                    *(copyBgns++) = curFacBgn;
                    *(++factorStarts) = currentPosition;
                    factorNum++;
                    curFacNod = 0;
                }
            } else {
                *(copyBgns++) = curFacBgn;
                *(++factorStarts) = currentPosition;
                factorNum++;
                curFacNod = 0;
            }
        }
    }
}


void MakeSufTree(void) {
    int sufNode, sufLen;
    currentSymbol = (int) seq_factor_copy[currentPosition];
    prLeaf = 4 * curPrNod + currentSymbol;
    if (curFacNod < prSumpow4) {
        if (curFacNod < 0)
            MoveSecFacNode();
        else
            MovePrFacNode();
        lvList[currentPosition] = prefixTree[prLeaf];
        prefixTree[prLeaf] = currentPosition;
    } else
        EnterSecTree();

    sufNode = curPrNod = (prLeaf - prSumpow4) % prPow4 + prSumpow4;
    sufLen = PRIMEDEPTH;
    while (prefixTree[sufNode] == 0) {
        prefixTree[sufNode] = currentPosition;
        sufLen--;
        sufNode =
                (sufNode - sumpow4[sufLen]) % pow4[sufLen] + sumpow4[sufLen];
    }
}


void EnterSecTree(void) {
    int indx, actSmb;
    if ((curFacNod = prefixTree[prLeaf]) > 0) {
        for (indx = 1; indx <= 4; indx++)
            newList[indx] = 0;
        do {
            curSecLeaf = curFacNod;
            actSmb = seq_factor_copy[curSecLeaf + 1];
            if (newList[actSmb])
                newSecLeaf[actSmb] = lvList[newSecLeaf[actSmb]] =
                        curSecLeaf;
            else
                newList[actSmb] = newSecLeaf[actSmb] = curSecLeaf;
        } while ((curFacNod = lvList[curFacNod]) > 0);
        if (curFacNod < 0)
            TreatInitSecNode();
        else
            MakeInitSecNode();
    } else {
        if (curFacNod < 0)
            GoInitFacNode();
        else {
            if (FACTYPE == LZ_K_MISMATCH) {
                *(copyEnds++) = prevCopyPstn + 1;
                *(factorEnds++) = prevCopyPstn =
                prefixTree[prLeaf] = currentPosition;
                factorNum++;
                curFacNod = 0;
            } else {
                *(copyBgns++) = prevCopyPstn - PRIMEDEPTH;
                *(++factorStarts) = currentPosition - 1;
                prefixTree[prLeaf] = currentPosition;
                factorNum++;
                if ((prevCopyPstn = prefixTree[curFacNod = currentSymbol]) == 0) {
                    *(copyBgns++) = (*factorStarts);
                    *(++factorStarts) = currentPosition;
                    factorNum++;
                    curFacNod = 0;
                }
            }
        }
    }
}


void TreatInitSecNode(void) {
    int indx;
    facSecNod = (int *)
            (secTree[(-curFacNod) / SECARRSIZ] + (-curFacNod) % SECARRSIZ);
    for (indx = 1; indx <= 4; indx++) {
        if (newList[indx]) {
            lvList[newSecLeaf[indx]] = facSecNod[indx];
            facSecNod[indx] = newList[indx];
        }
    }
    prefixTree[prLeaf] = curFacNod;
    facSecLeaf = currentPosition;
    facSecDep = 1;
}


void MakeInitSecNode(void) {
    int nodIndx, indx;
    if ((nodIndx = (++numSecNod) % SECARRSIZ) == 0) {
        if ((secTree[++numSecArr] = (int (*)[5]) calloc(SECARRSIZ, sizeof(int[5])))
            == NULL) {
            printf("Can't allocate memory for %d-th secondary array\n", numSecArr + 1);
            exit(11);
        }
    }
    facSecNod = (int *) (secTree[numSecArr] + nodIndx);
    *facSecNod = curSecLeaf;
    for (indx = 1; indx <= 4; indx++) {
        if (newList[indx]) {
            facSecNod[indx] = newList[indx];
            lvList[newSecLeaf[indx]] = 0;
        }
    }
    curFacNod = prefixTree[prLeaf] = -numSecNod;
    facSecLeaf = currentPosition;
    facSecDep = 1;
}


void GoInitFacNode(void) {
    facSecNod = (int *)
            (secTree[(-curFacNod) / SECARRSIZ] + (-curFacNod) % SECARRSIZ);
    facSecLeaf = currentPosition;
    facSecDep = 1;
}


void MoveSecFacNode(void) {
    int indx, actSmb;
    if ((nxtSecLeaf = facSecNod[currentSymbol]) > 0) {
        facSecDep++;
        for (indx = 1; indx <= 4; indx++)
            newList[indx] = 0;
        do {
            curSecLeaf = nxtSecLeaf;
            actSmb = seq_factor_copy[curSecLeaf + facSecDep];
            if (newList[actSmb])
                newSecLeaf[actSmb] = lvList[newSecLeaf[actSmb]] =
                        curSecLeaf;
            else
                newList[actSmb] = newSecLeaf[actSmb] = curSecLeaf;
        } while ((nxtSecLeaf = lvList[nxtSecLeaf]) > 0);
        if (nxtSecLeaf < 0)
            TreatOldSecNode();
        else
            MakeNewSecNode();
    } else {
        if (nxtSecLeaf < 0)
            GoOldSecNode();
        else {
            facSecNod[currentSymbol] = facSecLeaf;
            if (FACTYPE == LZ_K_MISMATCH) {
                *(copyEnds++) = (*facSecNod) + facSecDep;
                *(factorEnds++) = prevCopyPstn = currentPosition;
                factorNum++;
                curFacNod = 0;
            } else {
                *(copyBgns++) = (*facSecNod) - offset_primedepth;
                *(++factorStarts) = currentPosition - 1;
                factorNum++;
                if ((prevCopyPstn = prefixTree[curFacNod = currentSymbol]) == 0) {
                    *(copyBgns++) = (*factorStarts);
                    *(++factorStarts) = currentPosition;
                    factorNum++;
                    curFacNod = 0;
                }
            }
        }
    }
}


void TreatOldSecNode(void) {
    int indx, *nxtSecNod;
    nxtSecNod = (int *)
            (secTree[(-nxtSecLeaf) / SECARRSIZ] + (-nxtSecLeaf) % SECARRSIZ);
    for (indx = 1; indx <= 4; indx++) {
        if (newList[indx]) {
            lvList[newSecLeaf[indx]] = nxtSecNod[indx];
            nxtSecNod[indx] = newList[indx];
        }
    }
    facSecNod[currentSymbol] = nxtSecLeaf;
    facSecNod = nxtSecNod;
}


void MakeNewSecNode(void) {
    int nodIndx, indx, *nxtSecNod;
    if ((nodIndx = (++numSecNod) % SECARRSIZ) == 0) {
        if ((secTree[++numSecArr] = (int (*)[5]) calloc(SECARRSIZ, sizeof(int[5])))
            == NULL) {
            printf("Can't allocate memory for %d-th secondary array\n", numSecArr + 1);
            exit(12);
        }
    }
    nxtSecNod = (int *) (secTree[numSecArr] + nodIndx);
    *nxtSecNod = curSecLeaf;
    for (indx = 1; indx <= 4; indx++) {
        if (newList[indx]) {
            nxtSecNod[indx] = newList[indx];
            lvList[newSecLeaf[indx]] = 0;
        }
    }
    facSecNod[currentSymbol] = -numSecNod;
    facSecNod = nxtSecNod;
}


void GoOldSecNode(void) {
    int *nxtSecNod;
    nxtSecNod = (int *)
            (secTree[(-nxtSecLeaf) / SECARRSIZ] + (-nxtSecLeaf) % SECARRSIZ);
    facSecNod = nxtSecNod;
    facSecDep++;
}


void ComputeFinFactor(void) {
    if (FACTYPE == LZ_K_MISMATCH) {
        if (curFacNod < 0)
            *copyEnds = (*facSecNod) + facSecDep;
        else
            *copyEnds = prevCopyPstn + 1;
        *factorEnds = currentPosition;
    } else {
        if (curFacNod) {
            if (curFacNod < 0)
                *copyBgns = (*facSecNod) - offset_primedepth;
            else
                *copyBgns = prevCopyPstn + (*factorStarts) - wordLength;
            *(++factorStarts) = wordLength;
        } else
            factorNum--;
    }
}


void ComputeGlobalFactors(void) {
    int blockSizeBnd, fctIndx, subFctIndx, blockIndx, sumFact;
    actPerBnd = 0;
    blockSize = 2;

    if ((blockSizeBnd = maxNumErr + 1) > 2) {
        while (2 * blockSize <= blockSizeBnd)
            blockSize *= 2;
    }
    numFullBlock = (factorNum - 1) / blockSize;
    numGlobFact = numFullBlock + 1;

    finBlockSize = factorNum - numFullBlock * blockSize;

    maxSumFact = (int *) calloc(numGlobFact + 1, sizeof(int));
    fctIndx = (maxNumErr2 > blockSize) ? maxNumErr2 : blockSize;
    subFctIndx = fctIndx - maxNumErr2;

    while (fctIndx <= factorNum) {
        sumFact = factorEnds[fctIndx] - factorEnds[subFctIndx];
        if (sumFact > actPerBnd)
            actPerBnd = sumFact;
        blockIndx = fctIndx / blockSize;
        while (blockIndx * blockSize > subFctIndx) {
            if (sumFact > maxSumFact[blockIndx])
                maxSumFact[blockIndx] = sumFact;
            blockIndx--;
        }
        fctIndx++;
        subFctIndx++;
    }
    maxSumFact[numGlobFact] =
            factorEnds[factorNum] - factorEnds[factorNum - maxNumErr2];
}

// Method computes powers of 4 and sum_{i=0-n} 4^n, then assigns prPow4 and prSumpow4 using PRIMEDEPTH.
void MakePowof4(void) {
    int indx;
    sumpow4[0] = 0;
    pow4[0] = 1;
    for (indx = 1; indx < 16; indx++) {
        pow4[indx] = 4 * pow4[indx - 1];
        sumpow4[indx] = sumpow4[indx - 1] + pow4[indx - 1];
    }
    prPow4 = pow4[PRIMEDEPTH];
    prSumpow4 = sumpow4[PRIMEDEPTH];
}


int ComputeSizeFact(void) {
    int depth, sumLen = 0;
    if (FACTYPE == LZ_K_MISMATCH) {
        for (depth = 1; depth < 16; depth++) {
            if ((sumLen + depth * pow4[depth]) >= wordLength)
                return (sumpow4[depth] + (wordLength - sumLen) / depth);
            else
                sumLen += depth * pow4[depth];
        }
    } else {
        for (depth = 0; depth < 15; depth++) {
            if ((sumLen + depth * pow4[depth + 1]) >= wordLength)
                return (sumpow4[depth + 1] + (wordLength - sumLen) / depth);
            else
                sumLen += depth * pow4[depth + 1];
        }
    }
    /* this situation cannot happen unless the sequence length is about 20 10^9 */
    fprintf(stderr, "the input sequence is too long\n"); // TODO: Exception conversion
    exit(100);
} 
