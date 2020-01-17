/* ===============================================
 * File Name: gaussian_process.c
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-05-16 10:13:41
 * =============================================== 
 */

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_integration.h>

// def some struct
struct DATA {
    int n;
    gsl_vector * x;
    gsl_vector * y;
    gsl_vector * sig;
    gsl_matrix * IMcov;
    double sigf;
    double lenf;
};

// global data
struct DATA obs;

/**************************************************/
/* Using Gauss method to integrate function */

double qgaus(double (*func)(double), double a, double b)
{
	int j;
	double xr,xm,dx,s;
	static double x[]={0.0,0.1488743389,0.4333953941,
		0.6794095682,0.8650633666,0.9739065285};
	static double w[]={0.0,0.2955242247,0.2692667193,
		0.2190863625,0.1494513491,0.0666713443};

	xm=0.5*(b+a);
	xr=0.5*(b-a);
	s=0;
	for (j=1;j<=5;j++) {
		dx=xr*x[j];
		s += w[j]*((*func)(xm+dx)+(*func)(xm-dx));
	}
	return s *= xr;
}

double qgaus_arg(double (*func)(double, double), double a, double b, double arg)
{
	int j;
	double xr,xm,dx,s;
	static double x[]={0.0,0.1488743389,0.4333953941,
		0.6794095682,0.8650633666,0.9739065285};
	static double w[]={0.0,0.2955242247,0.2692667193,
		0.2190863625,0.1494513491,0.0666713443};

	xm=0.5*(b+a);
	xr=0.5*(b-a);
	s=0;
	for (j=1;j<=5;j++) {
		dx=xr*x[j];
		s += w[j]*((*func)(xm+dx, arg)+(*func)(xm-dx, arg));
	}
	return s *= xr;
}

double qgaus_args(double (*func)(double, double, double), double a, double b, double arg1, double arg2)
{
	int j;
	double xr,xm,dx,s;
	static double x[]={0.0,0.1488743389,0.4333953941,
		0.6794095682,0.8650633666,0.9739065285};
	static double w[]={0.0,0.2955242247,0.2692667193,
		0.2190863625,0.1494513491,0.0666713443};

	xm=0.5*(b+a);
	xr=0.5*(b-a);
	s=0;
	for (j=1;j<=5;j++) {
		dx=xr*x[j];
		s += w[j]*((*func)(xm+dx, arg1, arg2)+(*func)(xm-dx, arg1, arg2));
	}
	return s *= xr;
}

double qgaus2d_fix(double (*func)(double, double), double a, double b, double c, double d)
{
	int i,j;
	double xr,xm,dx,s;
    double yr,ym,dy;
	static double x[]={0.0,0.1488743389,0.4333953941,
		0.6794095682,0.8650633666,0.9739065285};
	static double w[]={0.0,0.2955242247,0.2692667193,
		0.2190863625,0.1494513491,0.0666713443};

	xm=0.5*(b+a);
	xr=0.5*(b-a);
    ym=0.5*(d+c);
    yr=0.5*(d-c);
	s=0;
	for (j=1;j<=5;j++) {
		dx=xr*x[j];
        for(i=1;i<=5;i++){
            dy = yr*x[i];
            s+= w[j]* w[i] * ((*func)(xm+dx, ym+dy) +(*func)(xm-dx, ym-dy));
        }
	}
    s = s*xr*yr;
	return s;
}

/**************************************************/

// initial data
int initial_gp(int N, char * filename, double sigf, double lenf)
{
    int i;
    double x, y, sig;
    FILE * fp = NULL;

    obs.n = N;
    obs.sigf = sigf;
    obs.lenf = lenf;
    obs.x = gsl_vector_alloc(N);
    obs.y = gsl_vector_alloc(N);
    obs.sig = gsl_vector_alloc(N);
    obs.IMcov = gsl_matrix_alloc(N,N);
    obs.IMcov = gsl_matrix_alloc(N,N);

    // read in data
    fp = fopen(filename, "r");

    for (i = 0; i< N; i++){
        fscanf(fp, "%lf %lf %lf\n", &x, &y, &sig);
        gsl_vector_set(obs.x, i, x);
        gsl_vector_set(obs.y, i, y);
        gsl_vector_set(obs.sig, i, sig);
    }
    fclose(fp);

    return 0;
}

// kernel
double kernel(double x, double y) 
{
    return pow(obs.sigf, 2) * exp(-pow(x-y, 2)/2/pow(obs.lenf, 2));
}

// set_up covariance
int setup_gp()
{
    int i, j;
    double x,y,sig;

    for (i = 0; i< obs.n; i++){
        x = gsl_vector_get(obs.x, i);
        sig = gsl_vector_get(obs.sig, i);
        for (j = 0; j< obs.n; j++){
            y = gsl_vector_get(obs.x, j);
            if(i == j)
                gsl_matrix_set(obs.IMcov, i, i, sig*sig +kernel(x,y) );
            else
                gsl_matrix_set(obs.IMcov, i, j, kernel(x,y) );
        }
    }

    gsl_linalg_cholesky_decomp(obs.IMcov);
    gsl_linalg_cholesky_invert(obs.IMcov);

    return 0;
}

// print obs dat
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
            cov = gsl_matrix_get(obs.IMcov, i, i);
            printf("x = %7.3f; y = %7.3f; sig= %12.8f\n", x,y, sqrt(cov) );
        }
    }
    
    return 0;
}

/*=====================================*/
// rec mu
double rec_mu(double x)
{
    double mu, xi, yj, M;
    int i,j;

    mu = 0.0;

    for (i=0; i<obs.n; i++){
        xi = gsl_vector_get(obs.x, i);
        for (j=0; j<obs.n; j++){
            yj = gsl_vector_get(obs.y, j);
            M = gsl_matrix_get(obs.IMcov, i, j);
            mu += kernel(x, xi) * M * yj;
        }
    }
    return mu;
}

// rec covariance
double rec_covariance(double x, double y)
{
    double cov, xi, xj, M;
    int i,j;

    cov = 0.0;

    for (i=0; i<obs.n; i++){
        xi = gsl_vector_get(obs.x, i);
        for (j=0; j<obs.n; j++){
            xj = gsl_vector_get(obs.x, j);
            M = gsl_matrix_get(obs.IMcov, i, j);
            cov += kernel(x, xi) * M * kernel(xj, y);
        }
    }
    cov = kernel(x,y) - cov;
    return cov;
}

// rec covariance_over_mumu
double rec_covariance_over_mumu(double x, double y)
{
    return rec_covariance(x,y)/rec_mu(x)/rec_mu(y);
}

/*=====================================*/

// func: 1/mu
double func_one_over_mu(double x)
{
    return 1.0/rec_mu(x);
}

// integrate func: 1/mu
double integrate_one_over_mu(double x)
{
    return qgaus(&func_one_over_mu, 0.0, x);
}

/* covariance of integrate_one_over_mu */
// func: cov(x,y)/x^2/y^2

double func_cov_over_mumu(double x, double y)
{
    return rec_covariance_over_mumu(x,y)*func_one_over_mu(x)*func_one_over_mu(y);
}

double cov_int_one_over_mu(double x, double y)
{
    return qgaus2d_fix(&func_cov_over_mumu, 0.0, x, 0.0, y);
}

double cov_mu_int_one_over_mu(double x, double y)
{
    return qgaus_arg(&func_cov_over_mumu, 0.0, x, y);
}

/*=====================================*/

// func: mu* int 1/mu dx

double rec_mu_int_one_over_mu(double x)
{
    return rec_mu(x) * integrate_one_over_mu(x);
}

double rec_cov_mu_int_one_over_mu(double x, double y)
{
    double fi, fj, mui, muj;
    double cov00, cov01, cov10, cov11;

    mui = rec_mu(x);
    muj = rec_mu(y);

    fi = rec_mu_int_one_over_mu(x);
    fj = rec_mu_int_one_over_mu(y);

    cov00 = rec_covariance_over_mumu(x,y);
    cov01 = cov_mu_int_one_over_mu(y,x);
    cov10 = cov_mu_int_one_over_mu(x,y);
    cov11 = cov_int_one_over_mu(x,y);

    return fi*fj*cov00 - fi*mui*muj*cov01 - fj*mui*muj*cov10 +mui*muj*cov11;

}

/*=====================================*/

// free obs
int free_obs()
{
    gsl_vector_free(obs.x);
    gsl_vector_free(obs.y);
    gsl_vector_free(obs.sig);
    gsl_matrix_free(obs.IMcov);
    return 0;
}

/*=====================================*/

double lcdm_hz(double z, double omm, double h0)
{
    return h0*sqrt(omm*pow(1+z, 3) + 1-omm);
}

double lcdm_one_over_hz(double z, double omm, double h0)
{
    return 1.0/lcdm_hz(z, omm, h0);
}

double lcdm_dz(double z, double omm, double h0)
{
    return qgaus_args(&lcdm_one_over_hz, 0.0, z, omm, h0);
}

//===========================
int main(void)
{
    int status;
    int N = 31, i;
    double x;
    char * filename = "/home/ekli/myworks/cosmodata/OHD_31_CC.txt";
    double sigf=133.80863893, lenf=1.93610328;

    if((status = initial_gp(N,filename,sigf,lenf)) != 0)
        printf("Can not initial gp!\n");

    if((status = setup_gp()) != 0)
        printf("Can not initial gp!\n");

    // test_struct_print();
    //
    for (i=0; i< 41; i++){
        x = 2.0*i/40.0;
        if(i == 0)
            printf("%5s %8s %9s %8s %9s %9s %9s %9s %9s %9s \n", 
                    "z", "H(z)", "sig", "lc_H", "D", "sig", "lc_D", "HD", "sig", "lc_HD");

        //printf("%5.2f %8.3f %9.6f ", x, rec_mu(x), sqrt( rec_covariance(x,x) ) );
        //printf("%8.3f ", lcdm_hz(x, 0.3, rec_mu(0.0)));
        
        //printf("intg = %9.6f sig = %9.6f ", integrate_one_over_mu(x),
        //        sqrt( cov_int_one_over_mu(x,x) ) );
        //printf("%9.6f ", lcdm_dz(x, 0.3, rec_mu(0.0)) );

        printf("%9.6f %9.6f ", rec_mu_int_one_over_mu(x),
                sqrt( rec_cov_mu_int_one_over_mu(x,x) ) );
        printf("%9.6f \n", lcdm_dz(x, 0.3, rec_mu(0.0))*lcdm_hz(x, 0.3, rec_mu(0.0)) );
    }

    free_obs();

    return 0;
}
