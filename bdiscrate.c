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

#define BUF 8
#define GBUF 8
#define WBUF 8

#define HCNSTATES 2

typedef unsigned char boole;

typedef struct /* r_t, rate type: convenience struct for the two rates and the mxextent */
{
    float s1r, s2r;
    int mxt;
} r_t;

typedef enum /*  simple enum to hold our two bases */
{
    A,
    B
} base;

typedef struct /* our definition of a site is a struct */
{
    int currp; /* curr num of jumps: also uses as index to posarr[], i.e. gives the size of posarr */
    int jbf; /* buffer allow for jumps */
    float mxdisp; /* current maximum displacement */
    float *posarr; /* the array of position values: to be a cumulative */
    base *brec; /* the array of past and current bases: one less than total of posarr */
    base latestb; /* the next base */
    char starsymb; /* the starting symbol (character) */
    char endsymb; /* the end symb */
} sitedef;/* our definition of a site is a struct */ 

void usage(char *progname)
{
    printf("Usage: \"%s\". This is a two state discrete exponential time simulator.\n", progname);
    printf("Please supply 4 or 5 arguments:\n");
    printf("\t1) Substitution rate of first state, either as fraction or float. Eg. 1/20 or .05 would be a rate of 1 in 20.\n");
    printf("\t2) Substitution rate of second state, same format as above.\n");
    printf("\t3) Number of sites, i.e. number of separate replications.\n");
    printf("\t4) Number of minimum extent units. I.e duration of each replication, it's calculated by");
    printf("\t\tThe smallest substrate by average will occur once per 1 extent unit,which is an int, not a float.\n");
    printf("\t5) specified random number seed (optional, if not supplied, a random seed will be generated).\n");
}

r_t *fill_rt(char *s1, char *s2);
{
   r_t *rts=malloc(sizeof(r_t));
   int n, d1=0, d2=0;
   char *tstr;
   char ttstr[32]={0};
   float min = 9999999.0;

    /* Note, this function turned out to be tricky because we allow the rate to be given as fraction or float
     * and even one rate may be a float, and the other a fraction */

   /* first string */
   tstr=strchr(s1, '/');
   if(!tstr) {
       rts->s1r=atof(s1);
       if(rts->s1r<min)
            min=rts->s1r;
   } else {
       strncpy(ttstr, s1, (tstr-s1)*sizeof(char));
       n=atoi(ttstr);
       strcpy(ttstr, tstr+1);
       d1=atoi(ttstr);
       rts->s1r=(float)n/d1;
   }

   /* second string */
   tstr=strchr(s2, '/');
   if(!tstr) {
       rts->s2r=atof(s2);
       if(rts->s2r<min)
           min=rts->s2r;
   } else {
       strncpy(ttstr, s2, (tstr-s2)*sizeof(char));
       n=atoi(ttstr);
       strcpy(ttstr, tstr+1);
       d2=atoi(ttstr);
       rts->s2r=(float)n/d2;
   }

   int mxt2=(int)(.5+1./min);
   if( d1 & d2) 
       rts->mxt=(d1>d2)? d1:d2;
   else if( (d1 & (d2==0) ) | (d2 & (d1==0) ) )
       if(rts->mxt < mxt2)
           rts->mxt = mxt2;
   else
       rts->mxt=mxt2;

    return rts;
}

sitedef *crea_sd(int numsites)
{
    int i;
    sitedef *sitearr=malloc(sizeof(sitedef) * numsites);
    for(i=0;i<numsites;++i) {
        sitearr[i].mxdisp=0.0;
        sitearr[i].currp=0;
        sitearr[i].jbf=BUF;
        sitearr[i].posarr=calloc(sizeof(float), sitearr[i].jbf);
        sitearr[i].brec=calloc(sizeof(base), sitearr[i].jbf);
    }
    return sitearr;
}

void summarysites(sitedef *sitearr, int numsites, int nstates, char *idstrng)
{
    int i;
    int *nucrec=calloc(sizeof(int), nstates);
    for(i=0;i<numsites;++i)
        switch(sitearr[i].endsymb) {
            case 'A':
                nucrec[0]++; break;
            case 'B':
                nucrec[1]++; break;
        }

    float avgnj=.0;
    for(i=0;i<numsites;++i) {
        avgnj+=sitearr[i].currp-1; /* why -1? the first currp is not a change, it's the initial state */
#ifdef DBG
        int j;
        printf("site:%3d pos) ", i); 
        for(j=0;j<sitearr[i].currp;++j) 
            printf("%.4f ", sitearr[i].posarr[j]); 
        printf("\n"); 
#endif
    }
    avgnj /=numsites;
    printf("%s.. avg #jumps %.2f:\n", idstrng, avgnj); 
    for(i=0;i<nstates;++i) 
        printf("%6i ",nucrec[i]); 
    printf("\n"); 
    for(i=0;i<nstates;++i) 
        printf("%3.4f ",(float)nucrec[i]/numsites); 
    printf("\n"); 
    free(nucrec);
}

void sitesubproc(sitedef* sites, r_t *rts, int nstates, int numsites, char symb, float lenc, int rsee)
{
    int i;
    /* what's the starting symbol? All sites given the same one */
    for(i=0;i<numsites;++i) {
        switch (symb) {
            case 'A':
                sites[i].latestb=A; /*  */
                sites[i].starsymb='A'; /*  */
                break;
            case 'B':
                sites[i].latestb=B; /*  */
                sites[i].starsymb='B'; /*  */
                break;
        }
    }

    float ura; /*  variable to hold one uniform random variable 0-1 */
    srand(rsee);
    base currba;
    float srate;
    for(i=0;i<numsites;++i) {
        sites[i].latestb = sites[i].brec[0] = A;
        while(1) { /* infinite loop to be broken out of when maxlen reaches a certain condition */
            ura= (float)rand()/RAND_MAX;
            currba = sites[i].latestb;
            if(currba==A) {
            sites[i].brec[sites[i].currp + 1] = B;
            srate=rts->s1r;
            } else {
            sites[i].brec[sites[i].currp + 1] = A;
            srate=rts->s2r;
            }
            sites[i].posarr[sites[i].currp + 1] = sites[i].posarr[sites[i].currp] + (-1.0/srate)*log1p(-ura); 

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
            if(sites[i].mxdisp >= lenc)
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
                sites[i].endsymb='B'; break;
        }
    }
}

int main(int argc, char *argv[])
{
    int i, rsee;
    /* argument accounting */
    if((argc!=6) & (argc!=5)) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    } else if (argc==5) { /* no seed given, use a random one */
        struct timeval tnow;
        gettimeofday(&tnow, NULL);
        rsee=(int)((tnow.tv_usec/100000.0 - tnow.tv_usec/100000)*RAND_MAX);
    } else
        rsee=atoi(argv[5]);
    printf("rsee=%d\n", rsee); 

    nstates=HCNSTATES;
    r_t *rts=fill_rt(argv[1], argv[2]);
#ifdef
    printf("rts: %.4f %.4f %d\n", rts->s1r, rts->s2r, rts->mxt); 
#endif
    int numsites=atoi(argv[3]);
    int lenc=atoi(argv[4])*rts->mxt;

#ifdef DBG
    int j;
    for(i=0;i<nstates;++i) {
        for(j=0;j<nstates;++j) 
            printf("%f ", mat[i*nstates+j]);
        printf("\n"); 
    }
#endif

    sitedef *sitearr=crea_sd(numsites);
    sitesubproc(sitearr, rts, nstates, numsites, 'A', lenc, rsee);
    summarysites(sitearr, numsites, nstates, "Final dist: -1/-rate log1p(-ura)");

    for(i=0;i<numsites;++i) {
        free(sitearr[i].posarr);
        free(sitearr[i].brec);
    }
    free(sitearr);
    free(rts);
    return 0;
}
