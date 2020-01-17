/* ===============================================
 * File Name: subroutines.c
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-05-19 10:18:37
 * =============================================== 
 */

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

// def some struct
struct DATA {
    int n;
    gsl_vector * x;
    gsl_vector * y;
    gsl_matrix * cov;
    double sigf;
    double lenf;
};

// global data
struct DATA obs;

// initial data
int initial_gp(int N, char * filename, double sigf, double lenf)
{
    int i, j;
    double x, y, sig;
    FILE * fp = NULL;

    obs.n = N;
    obs.sigf = sigf;
    obs.lenf = lenf;
    obs.x = gsl_vector_alloc(N);
    obs.y = gsl_vector_alloc(N);
    obs.cov = gsl_matrix_alloc(N,N);

    // read in data
    fp = fopen(filename, "r");

    for (i = 0; i< N; i++){
        fscanf(fp, "%lf %lf %lf\n", &x, &y, &sig);
        gsl_vector_set(obs.x, i, x);
        gsl_vector_set(obs.y, i, y);
        for (j = 0; j<N; j++){
            gsl_matrix_set(obs.cov, i, j, 0.0 );
        }
        gsl_matrix_set(obs.cov, i, i, sig*sig);
    }

    fclose(fp);

    return 0;
}

int test_struct_print()
{
    int i;
    double x,y,cov;

    {
        printf("sigf is %5.3f\n", obs.sigf);
        printf("lenf is %5.3f\n", obs.lenf);
        printf("numb is %d\n", obs.n);
        for (i=0; i<obs.n; i++){
            x = gsl_vector_get(obs.x, i);
            y = gsl_vector_get(obs.y, i);
            cov = gsl_matrix_get(obs.cov, i, i);
            printf("x = %7.3f; y = %7.3f; sig= %7.3f\n", x,y, sqrt(cov) );
        }
    }
    
    return 0;
}


// kernel
double kernel(double x, double y) 
{
    return pow(obs.sigf, 2) * exp(-pow(x-y, 2)/2/pow(obs.lenf, 2));
}

//double rec_mu_x(double x, )
//

// free obs
int free_obs()
{
    gsl_vector_free(obs.x);
    gsl_vector_free(obs.y);
    gsl_matrix_free(obs.cov);
    return 0;
}

//===========================
int main(void)
{
    int status;
    int N = 31;
    char * filename = "/home/ekli/myworks/cosmodata/OHD_31_CC.txt";
    double sigf=133.80863893, lenf=1.93610328;

    printf("Start");

    test_struct_print();

    if((status = initial_gp(N,filename,sigf,lenf)) != 0)
        printf("Can not initial gp!\n");

    test_struct_print();

    return 0;
}
