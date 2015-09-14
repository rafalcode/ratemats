/* helper program to print out exponential waiting time variates, using Knuth's method */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

float *nxta(float rateparam, float *rfa, int nsteps)
{
    int i;
    float *erfa=malloc(nsteps*sizeof(float));
    for(i=0;i<nsteps;++i) 
        erfa[i] = -log1p(-rfa[i]) / rateparam;
    return erfa;
}

int main(int argc, char *argv[])
{
    /* argument accounting: remember argc, the number of arguments, _includes_ the executable */
    if(argc!=4) {
        printf("Error. Pls supply 2 arguments: 1) int random seed 2) rate value (float or pure fraction) 3) Number of steps (duration of simulation).\n");
        exit(EXIT_FAILURE);
    }

    float rateparam;
    srand(atoi(argv[1]));
    char ttstr[32]={0};
    char *tstr=strchr(argv[2], '/');
    int n, d;
    if(!tstr)
        rateparam=atof(argv[2]);
    else {
        strncpy(ttstr, argv[2], (tstr-argv[2])*sizeof(char));
        ttstr[tstr-argv[2]]='\0';
        n=atoi(ttstr);
        strcpy(ttstr, tstr+1);
        d=atoi(ttstr);
        rateparam=(float)n/d;
    }
    int nsteps=atoi(argv[3]);
    float summit;
    int dsummit;
    int i;

    float *rfa=malloc(nsteps*sizeof(float));
    for(i=0;i<nsteps;++i) 
        rfa[i]=(float) random() / (RAND_MAX + 1.);

    float *erfa=nxta(rateparam, rfa, nsteps);

    printf("Random numbers drawn:\n"); 
    for(i=0;i<nsteps;++i) /* random numbers */
        printf("%.3f ", rfa[i]);
    printf("\n"); 
    summit=.0;
    printf("Conversion to expectation times:\n"); 
    for(i=0;i<nsteps;++i) /* waiting times */
        printf("%.3f ", erfa[i]);
    printf("\n"); 
    printf("Cumulative expectation times:\n"); 
    for(i=0;i<nsteps;++i) { /* cumulative */
        summit += erfa[i];
        printf("%.3f ", summit);
    }
    printf("\nAvg floating point waiting time: %f\n", summit/(float)nsteps);
    printf("Discrete version\n"); 
    for(i=0;i<nsteps;++i)
        printf("%3d ", (int)(.5+erfa[i]));
        // printf("%3d ", (int)ceilf(erfa[i]));
    printf("\n"); 
    for(i=0;i<nsteps;++i) {
        dsummit+=(int)(.5+erfa[i]);
        // dsummit+=(int)ceilf(erfa[i]);
        printf("%3d ", dsummit);
    }
    printf("\nAvg discrete waiting time: %f\n", dsummit/(float)nsteps);

    free(rfa);
    free(erfa);
    return 0;
}
