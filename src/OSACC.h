#ifndef __MAIN_H__
#define __MAIN_H__

#include "random.h"
#include "perm.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include<math.h>
#define BUFFER 50000
#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;

#define NSTACK 50

// global variables
Perm myperm;
Random r;

void permutation1(float allele1[],float allele2[],float pallele1[],float pallele2[],float status[],float pstatus[],int numind);
void permutation2(float allele1[],float allele2[],float pallele1[],float pallele2[],float status[],float pstatus[],int numind);
void permutation3(float allele1[],float allele2[],float pallele1[],float pallele2[],float status[],float pstatus[],int numind);
void permutation4(float allele1[],float allele2[],float pallele1[],float pallele2[],float status[],float pstatus[],int numind);

void make_ped_file(char *pedfile,int index,int* num_row,int num_mrk);
int read_dat_data(char datafile[],char *marker_names[]);




#ifdef __cplusplus
extern "C"
#endif
{
void indexx(unsigned long n, float arr[], unsigned long indx[]);
void sort3(unsigned long n, float ra[], float rb[], float rc[]);
void sort3wholedata(unsigned long n, float ra[], float rb[], float rc[],float rd[]);
void crank(unsigned long n,float v[],float w[],float u[],int *k);
void crankwholedata(unsigned long n,float v[],float w[],float s[],float u[],int *k);
void cntab1(int **nn, int ni, int nj, float *chisq, float *df, float *prob,float *cramrv, float *ccc);
void fexact(int *nrow, int *ncol, double *table, int *ldtabl,double *expect, double *percnt, double *emin, double *prt,
       double *pre, /* new in C : */ int *workspace);

}
#ifdef __cplusplus

#endif
#endif
