#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define nbin 50
#define PI 3.14159265359

int main(void)
{
  int nprocs, a, b, i, j, nstep, nframe, is, ii, il; 
  int *nary, nglob, currCount;
  int *tary, histgram[nbin];
  double *rary, *vary;
  double rho, fac, gl[3], RC, DRC, rr[3], dist, gr, cn; 
  FILE *fp; 
  char fname[20];

  // reset gr array
  for(a=0; a<nbin; a++) histgram[a]=0; 

  printf("Input: nprocs, first, interval, last\n"); 
  scanf("%d %d %d %d", &nprocs, &is, &ii, &il);
  printf("nprocs = %d, first = %d, interval = %d, last = %d\n",nprocs, is,ii,il); 

  for(nframe=0, nstep=is; nstep<=il; nstep+=ii, nframe++) {
     sprintf(fname, "data/pmd%06d",nstep);

     fp = fopen(fname,"ro");
     fscanf(fp,"%d %d %lf %lf %lf\n",&nglob,&currCount,&gl[0],&gl[1],&gl[2]);
 
     // initialization, compute cutoff length, allocate memory
     if(nstep==is) {

       RC  = gl[0]; 
       RC  = gl[1] > RC ? gl[1] : RC; 
       RC  = gl[2] > RC ? gl[2] : RC; 
       RC  = 0.5*RC; 
       DRC = RC/nbin; 
       rho = nglob/(gl[0]*gl[1]*gl[2]);

       printf("nglob = %d currCount = %d gl[0:2] = %lf %lf %lf\n",nglob,currCount,gl[0],gl[1],gl[2]);
       printf("RC = %f, DRC = %f, nbin = %d, rho = %f\n", RC, DRC, nbin, rho);

       nary=malloc(sizeof(long int)*nprocs);
       tary=malloc(sizeof(int)*nglob);
       rary=malloc(sizeof(double)*3*nglob);
       vary=malloc(sizeof(double)*3*nglob);
     }

     //printf("--- %s ---\n",fname);

     for(a=0; a<nprocs; a++) fscanf(fp,"%d",&nary[a]);
     for(a=0; a<nglob; a++) {
        fscanf(fp,"%d %lf %lf %lf %lf %lf %lf\n",
        &tary[a], &rary[3*a],&rary[3*a+1],&rary[3*a+2],&vary[3*a],&vary[3*a+1],&vary[3*a+2]);
     }

     // compute gr
     for(i=0; i<nglob; i++){
        for(j=0; j<nglob; j++) {
           for(dist=0.0, a=0; a<3; a++) {
              // compute interatomic distance
              rr[a]=rary[3*j+a]-rary[3*i+a]; 

              // apply periodic boundary condition
              if(rr[a]>=0.5*gl[a]) rr[a]-=gl[a];
              if(rr[a]<-0.5*gl[a]) rr[a]+=gl[a];

              dist+=rr[a]*rr[a];
            }
            // get interatomic distance for i-j pair
            dist=sqrt(dist);
            // if dist is less than the cutoff, increment histgram array
            if(dist<RC) {
               b = dist/DRC; 
               histgram[b]++;
            }
        }
     }
     
     fclose(fp);
  }

  // finalize, print gr, deallocate arrays 
  for(cn=0.0,a=1; a<nbin; a++){
     dist=a*DRC; 
     fac=4.0*PI*dist*dist*DRC*rho;
     gr=histgram[a]/fac/nframe/nglob;
     cn+=histgram[a]/nframe/nglob;
     printf("%f %f %f\n", dist, gr, cn); 
  }

  free(nary); 
  free(tary); 
  free(rary); 
  free(vary); 
  exit(0); 
}
