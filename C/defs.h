/* defs.h : Macros */
/*******************/

/* All global definitions are here */

#include <stdio.h>
#include <stdlib.h>
#include "queue.h"

#ifndef __REPS_DEFS
#define __REPS_DEFS

#define ACGT_FILE_FASTA   1
/* #define ACGT_SEQ_PLAIN    2 */
/* #define ACGT_FILE_PLAIN   3 */
#define ASCII_FILE_SSPACE 5 
/* #define ASCII_SEQ         6 */
/* #define ASCII_FILE        7 */

/* Constants */
/* #define MAXDISPLAY 2000 */

/* Factorization options */

#define WITH_OVERLAP 0
#define WITHOUT_OVERLAP 1

/* Used for parameter FACTYPE in roman/factorizeforGDR.c */
/* Lempel-Ziv factorization for k-mismatch */
#define LZ_K_MISMATCH 0

/* s-factorization for exact case */
#define S_FACT_EXACT_MATCH 1


/* Globals */

//extern int verbose ;         /* This externs are set in main.c */
//extern int verbose_long ;

/* extern int allow_non_acgt ;  // This one is in dawg.c  */

#define YES  1
#define NO   0

#define VIRTUAL_STEP ((LastWindow) ? 2*step : step)

/* Types */

struct s_limits
{
  int err_number ;
  int min_period ;
  int max_period ;        
  float min_exponent ;
  int max_exponent ;       /* Not implemented */
  int min_size ;
  int max_size ;           
} ;

typedef struct s_limits limits ;

struct s_sfactorization
{
  int nbfactors ;
  int *factor ;
  int *previous ;
} ;

typedef struct s_sfactorization *sfactorization ;

struct s_repetition
{
  char type ;
  int period ;
  int initpos ;
  int endpos ;
} ;

typedef struct s_repetition repetition ;

struct s_listreps
{
  repetition rep ;
  struct s_listreps *next ;
  
} ;

typedef struct s_listreps *listreps ;


typedef struct s_detected_repeat detected_repeat;

struct s_detected_repeat {
    int period;
    int size;
    float errorRate;
    int start;
    int end;
    char *motif;
};


struct s_repeat_queue_node {
    detected_repeat rpt;
    TAILQ_ENTRY(s_repeat_queue_node) nodes;
};

typedef struct s_repeat_queue_node repeat_queue_node;

// This typedef creates a head_t that makes it easy for us to pass pointers to
// head_t without the compiler complaining.
typedef TAILQ_HEAD(head_s, s_repeat_queue_node) head_t;

/* for accounting repetitions crossing the bounds
 * between steps
 */
struct s_transreps
{
  int period ;
  struct s_transreps *next ;
  
} ;

typedef struct s_transreps *transreps ;


/* Now the functions */
/* ----------------- */

/* let2num */

#define LET_TO_NUM(c) ((unsigned char) (c))

#endif /* __REPS_DEFS */


#define PROCESS_ERROR(output_file,output_file_name,error_code)   \
{                                                            \
  if (xmloutput == YES)                                      \
    {                                                        \
      fprintf(output_file,"<errorcode>%d</errorcode>\n",error_code);  \
      fprintf(output_file,"</mreps>\n");                      \
      printf("error %d is registered in xml file %s\n",error_code,output_file_name);  \
    }                                                        \
}

#define LIMIT_N_PROPORTION 0.05
