#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#define N 256
double* linspace(double a, double b, int n)
{
    double c; int i;
    double *u = malloc(n*sizeof(double));
    c = (b - a)/(n - 1);
    for(i = 0; i < n - 1; ++i)
        u[i] = a + i*c;
    u[n-1] = b;
    return u;
}
double func(double x, double y){
    return exp(x+y)*((x*x + 3*x)*(y*y - y) + (y*y + 3*y)*(x*x-x));
}
double sol(double x, double y){
    return exp(x+y)*(x*x - x)*(y*y - y);
}
double finderr(int m, int n, double a[m][n], double b[m][n]){
    int i,j;
    double maxi = -1e13;
    for(i=1; i<m-1; i++){
        for(j=1; j<n-1; j++)
            maxi = (fabs(a[i][j]-b[i][j])>maxi)?fabs(a[i][j]-b[i][j]):maxi;
    }
    return maxi;
}
int main(int argc, char *argv[])
{
    double threshold; int jet0=0, jet1=0, iterCount;
    if (!(argc==3)){
        printf("Error Wrong Argument Count. Assuming threshold 1e-4 and iteration limit 4000.\n");
        threshold = 1e-4; iterCount = 4000;
    }
    else{
      jet0 = sscanf(argv[1], "%i" , &iterCount);
      jet1 = sscanf(argv[2], "%lg" , &threshold);
    }
    if ((jet0==-1) || (jet1==-1)){
        printf("Error Wrong Argument Type\n"); return -1;
    }
    int i,j,k = 0;
    double *x, *y;
    double f[N+1][N+1];
    double u[N+1][N+1];
    double v[N+1][N+1];
    double s[N+1][N+1];
    double b[N+1][N+1];
    x = linspace(0.0,1.0,N+1);
    y = linspace(0.0,1.0,N+1);
    for(i=0; i<N+1; i++){
        for(j=0; j<N+1; j++){
            f[i][j] = func(x[i],y[j]);
            u[i][j] = 0;
            v[i][j] = 0;
            s[i][j] = sol(x[i],y[j]);
            b[i][j] = 0;
        }
    }

    double error, maximum;
    while (1){
      
      if((k%10)==0){
          if (k!=0){
          error = finderr(N+1,N+1,b,u);
          for(i=0; i<N+1; i++){
            for(j=0; j<N+1; j++)
              b[i][j] = u[i][j];
          }
          maximum = finderr(N+1,N+1,u,s);
          //printf("%i\t%e\t%e\n",k,maximum,error);
          if ((k >= iterCount) && (error<threshold)){
            printf("d = %e\tc = %e\n\n",error,maximum);
            break;
          }
          }
      }
      k++;  
        for(i=1; i<N; i++){
            for(j=1; j<N; j++){
                v[i][j] = 0.25*(u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] - 1.0/(N*N)*f[i][j]);
            }
        }
        //double maximum = finderr(N+1,N+1,u,v);
        //double error = finderr(N+1,N+1,u,s);
        //printf("%i\t%e\t%e\n",k,maximum,error);
        memcpy(u, v, (N+1)*(N+1)*sizeof(double));
    }
    
    
    //for(i = 0; i < N+1;i++) printf("%i\t%.20f\n",i,f[i][200]);
    free(x); free(y);
    return 0;
}