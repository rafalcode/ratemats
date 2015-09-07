/*  Returning to this after a mont or two, it was incredibly hard to work out
 *  Let's see: a certain site undergoes changes of base. Each change can be seen as a jump, and we don't know how many
 *  jumps there will per site. Normally the number of positions coming from the jumps would be one more than 
 *  the number of jumps, but in this cas,e the final position is actually known, because cute off the sites development
 *  after a fixed time. This means that the number of positions and the number of changes (i.e. bases adopted) is the same */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

#define NUMSITES 100 /*  number of sites running*/
#define NSTATES 4 /*  number of states our sites can be in: 4 for DN, 20 for AA, etc. */
#define LENC 0.1 /* length of course */
#define BUF 8

#define SETSEED 5 /*  if -DUNPREDRA is not fed in */

typedef enum /*  simple enum to hold our four bases */
{
    A,
    C,
    G,
    T
} base;

typedef struct /* our definition of a site is a struct */
{
    int currp; /* curr num of jumps: also uses as index to posarr[] */
    int jbf; /* buffer allow for jumps */
    float mxdisp; /* current maximum displacement */
    float *posarr; /* the array of position values: to be a cumulative */
    base *brec; /* the array of past and current bases: one less than total of posarr */
    base latestb; /* the next base */
    char starsymb; /* the starting symbol (character) */
    char endsymb; /* the end symb */
} sitedef;/* our definition of a site is a struct */ 

int *memind(int n) /* return an n x n-1 matrix of indices which exclude the diags */
{
    /* build our excluding-diag-index array first */
    int i, j, k, m, l=0;
    int *excind = calloc(n*(n-1), sizeof(int));

    for(i=0;i<n;++i)
        for(j=0;j<n;++j) { 
            k=n*i+j;
            m=k%n;
            if(k!=i*n+i)
                excind[l++] = m; 
        } /* excind is built */
    return excind;
}

base getnextrbase(float *nsf, int n, int *excind, base sta) /* get our next random base */
{
    int i;
    base retbase;

    float r= (float)rand()/RAND_MAX * nsf[n-1]; /* multiplied by final entry in nsf .... therefore stretched, so to speak */

    for(i=0;i<n-1;++i)
        if(r<nsf[i+1]) { /* shouldn't this be if(r<nsf[sta*(n-1)+i+1])? nipe just one row */
            retbase = (base) excind[sta*(n-1) + i]; 
            break;
        }
    return retbase;
}

float *mat2cum(float *arr, int n, int *excind)
{
    int i, j; /* note cumulative array will calculate second values onward */

    float *cumarr;
    cumarr=calloc(n*n, sizeof(float)); 
    for(i=0;i<n;++i)
        for(j=1;j<n;++j)
            cumarr[n*i+j] = cumarr[n*i+j-1] + (arr[i*n+excind[i*(n-1) +j-1]] / -arr[i*n + i]);

    return cumarr;
}

void sitesubproc(sitedef* sites, char symb) /* the function for the site substitution process */
{
    int i;
#ifdef UNPREDRA
    struct timeval tnow;
    gettimeofday(&tnow, NULL);
    /*  use the five last digits of the current microsecond reading to generate a seed */
    unsigned int tseed=(unsigned int)((tnow.tv_usec/100000.0 - tnow.tv_usec/100000)*RAND_MAX);
    srand(tseed);
#else
    srand(SETSEED);
#endif

    /* what's the starting symbol? */
    for(i=0;i<NUMSITES;++i) {
        switch (symb) {
            case 'A':
                sites[i].latestb=A; /*  */
                sites[i].starsymb='A'; /*  */
                break;
            case 'C':
                sites[i].latestb=C; /*  */
                sites[i].starsymb='C'; /*  */
                break;
            case 'G':
                sites[i].latestb=G; /*  */
                sites[i].starsymb='G'; /*  */
                break;
            case 'T':
                sites[i].latestb=T; /*  */
                sites[i].starsymb='T'; /*  */
                break;
        }
    }

    float ar[16] = {-0.886, 0.190, 0.633, 0.063, 0.253, -0.696, 0.127, 0.316, 1.266, 0.190, -1.519, 0.063, 0.253, 0.949, 0.127, -1.329}; /*  this is a given: substitution matrix */

    int *excind=memind(4); /* a matrix of indices will be 4 rows by 3 columns */
    float *nsf = mat2cum(ar,4, excind); /* an array to hold the off diagonal entries, divided by the diagonal entry, all made negative */

    float ura; /*  variable to hold one uniform random variable 0-1 */
    base currba;
    for(i=0;i<NUMSITES;++i) {
        sites[i].latestb = sites[i].brec[0] = G;
        while(1) { /* infinite loop to be broken out of when maxlen reaches a certain condition */
            ura= (float)rand()/RAND_MAX;
            currba = sites[i].latestb;
            //            printf("currba: %u\n",currba);
            /*  armed with our ura, we can generate the Exp() */
            sites[i].posarr[sites[i].currp + 1] = sites[i].posarr[sites[i].currp] + (-1.0/ar[4*currba+currba])*log(1+ura); 
            // sites[i].posarr[sites[i].currp + 1] = sites[i].posarr[sites[i].currp] + (-1.0/-ar[4*currba+currba])*log1p(-ura); 
            sites[i].brec[sites[i].currp + 1] = getnextrbase(nsf + 4*currba, 4, excind, currba);  /* note a single row of nsf is "sent up"

            sites[i].latestb = sites[i].brec[sites[i].currp + 1];
            sites[i].mxdisp = sites[i].posarr[sites[i].currp + 1];

            sites[i].currp++;
            /* check posarr buf sz */
            if(sites[i].currp==sites[i].jbf-1) {
                sites[i].jbf += BUF;
                sites[i].posarr=realloc(sites[i].posarr, sites[i].jbf * sizeof(float));
                memset(sites[i].posarr+sites[i].jbf-BUF, 0, BUF*sizeof(float));
                sites[i].brec=realloc(sites[i].brec, sites[i].jbf * sizeof(base));
                memset(sites[i].brec+sites[i].jbf-BUF, 0, BUF*sizeof(base));
            }
            /*  breaking out when condition met */
            if(sites[i].mxdisp >= LENC)
                break; /*  this site has now crossed finishing line, go to next site*/
        }
        /*  rationalise posarr array size here */
        sites[i].posarr=realloc(sites[i].posarr, sites[i].currp * sizeof(float));
        sites[i].brec=realloc(sites[i].brec, sites[i].currp * sizeof(base));
        /*  FInd out at what symbol site finished at */
        switch(sites[i].brec[sites[i].currp-1]) {
            case 0:
                sites[i].endsymb='A'; break;
            case 1:
                sites[i].endsymb='C'; break;
            case 2:
                sites[i].endsymb='G'; break;
            case 3:
                sites[i].endsymb='T'; break;
        }
    }
    free(nsf);
    free(excind);
}

int main(int argc, char *argv[])
{
    int i;
    sitedef* sitearr; /*  declare */
    /* now, initialise */
    sitearr=malloc(sizeof(sitedef) * NUMSITES);
    for(i=0;i<NUMSITES;++i) {
        sitearr[i].mxdisp=0.0;
        sitearr[i].currp=0;
        sitearr[i].jbf=BUF;
        sitearr[i].posarr=calloc(sizeof(float), sitearr[i].jbf);
        sitearr[i].brec=calloc(sizeof(base), sitearr[i].jbf);
        /* latestb, starsymb and endsymb will all be set by the subroutine */
    }

    /*  OK, the race takes place here: */
    sitesubproc(sitearr, 'G');

    int *nucrec=calloc(sizeof(int), NSTATES);
    for(i=0;i<NUMSITES;++i)
        switch(sitearr[i].endsymb) {
            case 'A':
                nucrec[0]++; break;
            case 'C':
                nucrec[1]++; break;
            case 'G':
                nucrec[2]++; break;
            case 'T':
                nucrec[3]++; break;
        }

    printf("Final Distribution:\n"); 
    for(i=0;i<NSTATES;++i) 
        printf("%6i ",nucrec[i]); 
    printf("\n"); 
    for(i=0;i<NSTATES;++i) 
        printf("%3.4f ",(float)nucrec[i]/NUMSITES); 
    printf("\n"); 

    for(i=0;i<NUMSITES;++i) {
        free(sitearr[i].posarr);
        free(sitearr[i].brec);
    }
    free(sitearr);
    free(nucrec);
    return 0;
}
