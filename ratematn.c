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

typedef unsigned char boole;
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

void usage(void)
{
    printf("Usage. Pls supply 3 or 4 arguments:\n");
    printf("\t1) name of text file with matrix.\n");
    printf("\t2) number of sites.\n");
    printf("\t3) branch length limit, i.e. duration of simulation: 0.1 is small, 10 big, 1000 very big.\n");
    printf("\t4) (optional, if not supplied, a random seed will be generated) specified random number seed.\n");
}

typedef struct /* wseq_t */
{
    size_t *wln;
    size_t wsbuf;
    size_t quan;
    size_t lbuf; /* a buffer for the number of lines */
    size_t numl; /* number of lines, i.e. rows */
    size_t *wpla; /* words per line array: the number of words on each line */
} wseq_t;

wseq_t *create_wseq_t(size_t initsz)
{
    wseq_t *words=malloc(sizeof(wseq_t));
    words->wsbuf = initsz;
    words->quan = initsz;
    words->wln=calloc(words->wsbuf, sizeof(size_t));
    words->lbuf=WBUF;
    words->numl=0;
    words->wpla=calloc(words->lbuf, sizeof(size_t));
    return words;
}

void free_wseq(wseq_t *wa)
{
    free(wa->wln);
    free(wa->wpla);
    free(wa);
}

float *processinpf(char *fname, int *m, int *n)
{
    /* In order to make no assumptions, the file is treated as lines containing the same amount of words each,
     * except for lines starting with #, which are ignored (i.e. comments). These words are checked to make sure they contain only floating number-type
     * characters [0123456789+-.] only, one string variable is icontinually written over and copied into a growing floating point array each time */

    /* declarations */
    FILE *fp=fopen(fname,"r");
    int i;
    size_t couc /*count chars per line */, couw=0 /* count words */, oldcouw = 0;
    int c;
    boole inword=0;
    wseq_t *wa=create_wseq_t(GBUF);
    size_t bwbuf=WBUF;
    char *bufword=calloc(bwbuf, sizeof(char)); /* this is the string we'll keep overwriting. */

    float *mat=malloc(GBUF*sizeof(float));

    while( (c=fgetc(fp)) != EOF) {
        /*  take care of  */
        if( (c== '\n') | (c == ' ') | (c == '\t') | (c=='#')) {
            if( inword==1) { /* we end a word */
                wa->wln[couw]=couc;
                bufword[couc++]='\0';
                bufword = realloc(bufword, couc*sizeof(char)); /* normalize */
                mat[couw]=atof(bufword);
                couc=0;
                couw++;
            }
            if(c=='#') {
                while( (c=fgetc(fp)) != '\n') ;
                continue;
            } else if(c=='\n') {
                if(wa->numl == wa->lbuf-1) {
                    wa->lbuf += WBUF;
                    wa->wpla=realloc(wa->wpla, wa->lbuf*sizeof(size_t));
                    memset(wa->wpla+(wa->lbuf-WBUF), 0, WBUF*sizeof(size_t));
                }
                wa->wpla[wa->numl] = couw-oldcouw;
                oldcouw=couw;
                wa->numl++;
            }
            inword=0;
        } else if( (inword==0) && ((c == 0x2B) | (c == 0x2D) | (c == 0x2E) | ((c >= 0x30) && (c <= 0x39))) ) { /* deal with first character of new word, + and - also allowed */
            if(couw == wa->wsbuf-1) {
                wa->wsbuf += GBUF;
                wa->wln=realloc(wa->wln, wa->wsbuf*sizeof(size_t));
                mat=realloc(mat, wa->wsbuf*sizeof(float));
                for(i=wa->wsbuf-GBUF;i<wa->wsbuf;++i)
                    wa->wln[i]=0;
            }
            couc=0;
            bwbuf=WBUF;
            bufword=realloc(bufword, bwbuf*sizeof(char)); /* don't bother with memset, it's not necessary */
            bufword[couc++]=c; /* no need to check here, it's the first character */
            inword=1;
        } else if( (c == 0x2E) | ((c >= 0x30) && (c <= 0x39)) ) {
            if(couc == bwbuf-1) { /* the -1 so that we can always add and extra (say 0) when we want */
                bwbuf += WBUF;
                bufword = realloc(bufword, bwbuf*sizeof(char));
            }
            bufword[couc++]=c;
        } else {
            printf("Error. Non-float character detected. This program is only for reading floats\n"); 
            free_wseq(wa);
            exit(EXIT_FAILURE);
        }

    } /* end of big for statement */
    fclose(fp);
    free(bufword);

    /* normalization stage */
    wa->quan=couw;
    wa->wln = realloc(wa->wln, wa->quan*sizeof(size_t)); /* normalize */
    mat = realloc(mat, wa->quan*sizeof(float)); /* normalize */
    wa->wpla= realloc(wa->wpla, wa->numl*sizeof(size_t));

    *m= wa->numl;
    int k=wa->wpla[0];
    for(i=1;i<wa->numl;++i)
        if(k != wa->wpla[i])
            printf("Warning: Numcols is not uniform at %i words per line on all lines. This file has one with %zu.\n", k, wa->wpla[i]); 
    *n= k; 
    free_wseq(wa);

    return mat;
}

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

void sitesubproc(sitedef* sites, float *ar, int nstates, int numsites, char symb, float lenc, int rsee)
{
    int i;
    /* what's the starting symbol? */
    for(i=0;i<numsites;++i) {
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

    int *excind=memind(nstates); /* a matrix of indices will be 4 rows by 3 columns */
    float *nsf = mat2cum(ar, nstates, excind); /* an array to hold the off diagonal entries, divided by the diagonal entry, all made negative */

    float ura; /*  variable to hold one uniform random variable 0-1 */
    base currba;
    for(i=0;i<numsites;++i) {
        sites[i].latestb = sites[i].brec[0] = G;
        while(1) { /* infinite loop to be broken out of when maxlen reaches a certain condition */
            ura= (float)rand()/RAND_MAX;
            currba = sites[i].latestb;
            //            printf("currba: %u\n",currba);
            /*  armed with our ura, we can generate the Exp() */
            sites[i].posarr[sites[i].currp + 1] = sites[i].posarr[sites[i].currp] + (-1.0/ar[4*currba+currba])*log(1+ura); 
            // sites[i].posarr[sites[i].currp + 1] = sites[i].posarr[sites[i].currp] + (-1.0/-ar[4*currba+currba])*log1p(-ura); 
            sites[i].brec[sites[i].currp + 1] = getnextrbase(nsf + 4*currba, 4, excind, currba);  /* note: just a single row of nsf is "sent up" */

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
    int i, nstates, ncols, rsee;
    float lenc;
    /* argument accounting */
    if((argc!=5) & (argc!=4)) {
        usage();
        exit(EXIT_FAILURE);
    } else if (argc==4) { /* no seed given, use a random one */
        struct timeval tnow;
        gettimeofday(&tnow, NULL);
        rsee=(int)((tnow.tv_usec/100000.0 - tnow.tv_usec/100000)*RAND_MAX);
    } else
        rsee=atoi(argv[4]);
    lenc=atof(argv[3]);

    float *mat=processinpf(argv[1], &nstates, &ncols);
    if(nstates!=ncols) {
        printf("Error: only square matrices accepted.\n"); 
        exit(EXIT_FAILURE);
    }
#ifdef DBG
    for(i=0;i<nstates;++i) {
        for(j=0;j<nstates;++j) 
            printf("%f ", mat[i*nstates+j]);
        printf("\n"); 
    }
#endif
    int numsites=atoi(argv[2]);

    sitedef* sitearr; /*  declare */
    /* now, initialise */
    sitearr=malloc(sizeof(sitedef) * numsites);
    for(i=0;i<numsites;++i) {
        sitearr[i].mxdisp=0.0;
        sitearr[i].currp=0;
        sitearr[i].jbf=BUF;
        sitearr[i].posarr=calloc(sizeof(float), sitearr[i].jbf);
        sitearr[i].brec=calloc(sizeof(base), sitearr[i].jbf);
    }

    /*  OK, the race takes place here: */
    sitesubproc(sitearr, mat, nstates, numsites, 'G', lenc, rsee);

    int *nucrec=calloc(sizeof(int), nstates);
    for(i=0;i<numsites;++i)
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
    for(i=0;i<nstates;++i) 
        printf("%6i ",nucrec[i]); 
    printf("\n"); 
    for(i=0;i<nstates;++i) 
        printf("%3.4f ",(float)nucrec[i]/numsites); 
    printf("\n"); 

    for(i=0;i<numsites;++i) {
        free(sitearr[i].posarr);
        free(sitearr[i].brec);
    }
    free(sitearr);
    free(nucrec);
    free(mat);
    return 0;
}
