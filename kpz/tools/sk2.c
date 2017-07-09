#include <stdio.h>
#include <math.h>

#define tmax 100001
#define Hmax 200

int tin,fnum,count;
int i,j,k,ca,navar,np,ndum;
float cin,dum;
double W2[tmax+1],dH;
double w2a,w2,k1,k2,k3,k4;
float p,z,time;
double S,Q,norm;
double natl;
FILE *fopen(),*fp1, *fp0;
char fst[80];

main (int argc, const char* argv[]) {

/*
  printf(" how many files ? ");
  scanf("%d",&fnum);
*/
  fnum = 1;
  if(argc < 2)
  {
      printf(" avar time ? \n");
      scanf("%d",&navar);
  }
  else
      navar = atoi(argv[1]);
//  printf(" how many points ? ");
//  scanf("%d",&np);
  np = 100000;

//  fp0 = fopen("files","r");
  k2 = 0; k2=0; k3=0; dH=0;
  norm = np*fnum;

  for (j = 1; j <= fnum; j++) {

//    fscanf(fp0,"%s",&fst);
//    printf("%s\n",fst);
    if(argc < 3)
        fp1 = fopen("kpzSteady-12-w2","r");
    else
        fp1 = fopen(argv[2],"r");

    w2 = 0; w2a=0;
    for (i = 0; i <= tmax; i++) W2[i] = 0;
    count = 1;

    for (i = 1; i <= np; i++) {
      if(fscanf(fp1,"%d %f %f %f %d\n",&tin,&cin,&dum,&dum,&ndum) > 0 ) {
         if(tin > navar) {
           count ++;
           w2 = (long double)cin - 0.25;    
           W2[count] = w2;
           w2a += w2;
//           printf("%d %f\n",count,W1[count]);
         }
       }
      }
      fclose(fp1);
    }
//    fclose(fp0);

    w2a /= (double)count;
    printf("avar : %f\n",w2);
  
   for (i = 1; i <= count; i++) {

     dH = (W2[i] - w2);

     k1 += dH / (double)count / (double)fnum;
     k3 += (dH * dH *dH) / (double)count / (double)fnum;
     k2 += (dH * dH) / (double)count / (double)fnum;
     k4 += (dH * dH * dH *dH) / (double)count /(double)fnum;
 
    }
 
    S = k3 / pow(k2,1.5);
    Q = k4 / (k2 * k2) - 3.0;

//      S = k3 / (k1 * k1 * k1);
//      Q = k4 / (k1 * k1 * k1 * k1);


  printf("Skewness : %f, Kurtosis: %f\n",S,Q);
  
  
}
