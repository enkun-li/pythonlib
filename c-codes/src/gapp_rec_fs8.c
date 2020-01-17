/* ===============================================
 * File Name: gapp_rec_fs8.c
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-06-10 22:02:09
 * =============================================== 
 */

#include "common.h"
#include "utils.h"
#include "gapp.h"
#include "gapp_rec.h"

/* typedef fs8 data struct */
struct fs8_data {
    int n;
    gsl_vector *z;
    gsl_vector *f;
    gsl_matrix *cov;
};

struct fs8_data fs8obs;

/* @@@@@@@@@@@@@@@@@@@ */
/* some integrate func */
/* @@@@@@@@@@@@@@@@@@@ */

double rec_fs8_Ii_Hz(double z, double gam)
{
    double pars[] = {3*gam-1, 2*gam};
    return integrate_fx_over_fmu(0.0, z, &pars[0]);    
}

double rec_fs8_Mij_Hz(double x, double y, double gam)
{
    double pars[] = {y, 3*gam-1, 2*gam+1, 1.0};
    return integrate_fxcovxj_over_fmux(0.0, x, &pars[0]);    
}

double rec_fs8_Nij_Hz(double x, double y, double gam)
{
    double pars[] = {3*gam-1, 2*gam+1};
    return integrate_fxcovxy_over_fmuxy(0.0,x,0.0,y,&pars[0]);
}

/* @@@@@@@@@@@@@@@@@@ */
/* rec fs8 using H(z) */
/* @@@@@@@@@@@@@@@@@@ */
double rec_fs8_with_Hz(double z, double sig8, double ommh2, double gam)
{
    double pars[] = {3*gam, 2*gam};
    double ommH2g = pow(ommh2*10000, gam);

    return sig8 * ommH2g * func_one_over_mu(z, &pars[0])
        * exp( -ommH2g * rec_fs8_Ii_Hz(z, gam));
} 

double rec_cov_fs8_with_Hz(double zi, double zj, double sig8, double ommh2, double gam)
{
    double ommH2g = pow(ommh2*10000, gam);
    double fi, fj;
    double cov00, covx0, cov0y, covxy;
    double pars1[] = {0.0, 1.0};

    fi = rec_fs8_with_Hz(zi,sig8,ommh2,gam);
    fj = rec_fs8_with_Hz(zj,sig8,ommh2,gam);

    cov00 = func_covxy_over_mux_muy(zi,zj,&pars1[0]);
    covx0 = rec_fs8_Mij_Hz(zi,zj,gam);
    cov0y = rec_fs8_Mij_Hz(zj,zi,gam);
    covxy = rec_fs8_Nij_Hz(zi,zj,gam);

    return fi*fj * pow(2*gam, 2.0) * (cov00 
            - ommH2g*covx0 - ommH2g*cov0y 
            +ommH2g*ommH2g * covxy);
} 

/* @@@@@@@@@@@@@@@@@@@ */
/* some integrate func */
/* @@@@@@@@@@@@@@@@@@@ */

double rec_fs8_Ii_Ez(double z, double gam)
{
    double pars[] = {3*gam-1, 2*gam};
    return integrate_fxmu0_over_fmu(0.0, z, &pars[0]);    
}

double rec_fs8_Mij_Ez(double x, double y, double gam)
{
    double pars[] = {y, 3*gam-1, 2*gam+1, 1.0};
    return integrate_fxmu0covxj_over_fmux(0.0, x, &pars[0]);    
}

double rec_fs8_Nij_Ez(double x, double y, double gam)
{
    double pars[] = {3*gam-1, 2*gam+1};
    return integrate_fxmu0covxy_over_fmuxy(0.0,x,0.0,y,&pars[0]);
}

/* @@@@@@@@@@@@@@@@@@ */
/* rec fs8 using E(z) */
/* @@@@@@@@@@@@@@@@@@ */
double rec_fs8_with_Ez(double z, double sig8, double omm, double gam)
{
    double pars[] = {3*gam, 2*gam};
    //double pars2[] = {3*gam-1, 2*gam};
    double ommg = pow(omm, gam);

    return sig8 * ommg * func_mu0_over_mu(z, &pars[0])
        * exp(-ommg * rec_fs8_Ii_Ez(z, gam));
        //* exp(-ommg * integrate_fxmu0_over_fmu(0.0,z,&pars2[0]) );
}

double rec_cov_fs8_with_Ez(double zi, double zj, double sig8, double omm, double gam)
{
    double ommg = pow(omm, gam);
    double fi, fj;
    double cov00, covx0, cov0y, covxy;
    double pars1[] = {0.0, 1.0};
    //double parsx0[] = {zj,3*gam-1, 2*gam+1, 1.0};
    //double pars0y[] = {zi,3*gam-1, 2*gam+1, 1.0};
    //double pars2[] = {3*gam-1, 2*gam+1};

    fi = rec_fs8_with_Ez(zi,sig8,omm,gam);
    fj = rec_fs8_with_Ez(zj,sig8,omm,gam);

    cov00 = func_mu0_mu0_covxy_over_mux_muy(zi,zj,&pars1[0]);
    //covx0 = integrate_fxmu0covxj_over_fmux(0.0, zi, &parsx0[0]);
    //cov0y = integrate_fxmu0covxj_over_fmux(0.0, zj, &pars0y[0]);
    //covxy = integrate_fxmu0covxy_over_fmuxy(0.0, zi, 0.0, zj, &pars2[0]);
    covx0 = rec_fs8_Mij_Ez(zi,zj,gam);
    cov0y = rec_fs8_Mij_Ez(zj,zi,gam);
    covxy = rec_fs8_Nij_Ez(zi,zj,gam);

    return fi*fj * pow(2*gam, 2.0) * (cov00 
            - ommg*covx0 - ommg*cov0y 
            +ommg*ommg * covxy);
}

/* @@@@@@@@@@@@@@@@@@@@@@@@@ */
/* loglikelihood of fsigma_8 */
/* @@@@@@@@@@@@@@@@@@@@@@@@@ */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~ */

int initial_fs8_obs(int N, char *file_fz, char *file_cov)
{
    int i,j;
    double x, y, sig;
    FILE *fpf = NULL;
    FILE *fpc = NULL;

    // initial array and matrix
    fs8obs.n = N;
    fs8obs.z = gsl_vector_alloc(N);
    fs8obs.f = gsl_vector_alloc(N);
    fs8obs.cov = gsl_matrix_alloc(N,N);

    // open file
    fpf = fopen(file_fz, "r");
    fpc = fopen(file_cov, "r");

    // read in data
    for (i=0; i<N; i++) {
        fscanf(fpf, "%lf %lf\n", &x, &y);
        gsl_vector_set(fs8obs.z, i, x);
        gsl_vector_set(fs8obs.f, i, y);
        for (j=0; j<N; j++) {
            fscanf(fpc, "%lf\n", &sig);
            gsl_matrix_set(fs8obs.cov, i, j, sig);
        }
    }

    // close file
    fclose(fpf);
    fclose(fpc);

    return _SUCCESS_;
}

double fs8_loglikelihood_Hz(double sig8, double ommh2, double gam)
{
    double ommH2g = pow(ommh2*10000, gam);
    double chisq, z, zj, ff, ffs, fi, fj, c00, c01, c10, c11;
    gsl_vector *diffs8 = gsl_vector_alloc(fs8obs.n);
    gsl_matrix *cov = gsl_matrix_alloc(fs8obs.n,fs8obs.n);
    gsl_vector *A = gsl_vector_alloc(fs8obs.n);
    gsl_vector *B = gsl_vector_alloc(fs8obs.n);
    gsl_matrix *K = gsl_matrix_alloc(fs8obs.n,fs8obs.n);
    gsl_matrix *M = gsl_matrix_alloc(fs8obs.n,fs8obs.n);
    gsl_matrix *N = gsl_matrix_alloc(fs8obs.n,fs8obs.n);
    int i, j;
    double pars1[] = {3*gam, 2*gam};
    double pars2[] = {3*gam-1, 2*gam};
    double pars3[] = {0.0, 1.0};
    double pars4[] = {0.0, 3*gam-1, 2*gam+1, 1.0};
    double pars5[] = {3*gam-1, 2*gam+1};

    for (i=0; i<fs8obs.n; i++) {
        z = gsl_vector_get(fs8obs.z, i);
        gsl_vector_set(A, i, func_one_over_mu(z, &pars1[0]));
        gsl_vector_set(B, i, integrate_fx_over_fmu(0.0, z, &pars2[0]));
        for (j=0; j<=i; j++) {
            zj = gsl_vector_get(fs8obs.z, j);
            
            ff = func_covxy_over_mux_muy(z,zj,&pars3[0]);
            gsl_matrix_set(K, i, j, ff);
            gsl_matrix_set(K, j, i, ff);
            
            ff = integrate_fxcovxy_over_fmuxy(0.0, z, 0.0, zj, &pars5[0]);
            gsl_matrix_set(N, i, j, ff);
            gsl_matrix_set(N, j, i, ff);
        }
        
        for (j=0; j<fs8obs.n; j++) {
            zj = gsl_vector_get(fs8obs.z, j);
            
            pars4[0] = zj;
            ff = integrate_fxcovxj_over_fmux(0.0, z, &pars4[0]);
            gsl_matrix_set(M, i, j, ff);
        }
    }

    // rec fs8 mean
    for (i=0; i<fs8obs.n; i++) {
        ff = gsl_vector_get(A,i);
        ffs = gsl_vector_get(B,i);
        fi = sig8*ommH2g*ff*exp(-ommH2g * ffs);
        fj = gsl_vector_get(fs8obs.f, i);
        gsl_vector_set(diffs8, i, fi-fj);
    }

    // rec fs8 cov
    for (i=0; i<fs8obs.n; i++) {
        fi = gsl_vector_get(diffs8, i);
        for (j=0; j<=i; j++) {
            fj = gsl_vector_get(diffs8, j);
            c00 = gsl_matrix_get(K,i,j);
            c01 = gsl_matrix_get(M,i,j);
            c10 = gsl_matrix_get(M,j,i);
            c11 = gsl_matrix_get(N,i,j);
            ff = fi*fj*pow(2*gam, 2)*( c00 - ommH2g*c01 
                    - ommH2g*c10 + ommH2g*ommH2g * c11);
            ffs = gsl_matrix_get(fs8obs.cov, i, j);
            gsl_matrix_set(cov, i, j, ff+ffs);
            gsl_matrix_set(cov, j, i, ff+ffs);
        }
    }

    // inverse cov
    gsl_linalg_cholesky_decomp(cov);
    gsl_linalg_cholesky_decomp(cov);

    // chisq
    chisq = 0.0;
    for (i=0; i<fs8obs.n; i++) {
        fi = gsl_vector_get(diffs8, i);
        for (j=0; j<fs8obs.n; j++) {
            fj = gsl_vector_get(diffs8, j);
            c00 = gsl_matrix_get(cov, i, j);
            chisq += fi*c00*fj;
        }
    }

    // free the gsl vector and matrix
    gsl_vector_free(diffs8);
    gsl_vector_free(A);
    gsl_vector_free(B);
    gsl_matrix_free(cov);
    gsl_matrix_free(K);
    gsl_matrix_free(M);
    gsl_matrix_free(N);
    return -0.5*chisq;
}

double fs8_loglikelihood_Hz_t(double sig8, double ommh2, double gam)
{
    int i, j;
    double chisq, zi, zj, fi, fj, mij;
    gsl_vector *diffs8 = gsl_vector_alloc(fs8obs.n);
    gsl_matrix *Mcov = gsl_matrix_alloc(fs8obs.n, fs8obs.n);

    // initial dif fs8 & cov
    for (i=0; i<fs8obs.n; i++) {
        zi = gsl_vector_get(fs8obs.z, i);
        fi = rec_fs8_with_Hz(zi,sig8,ommh2,gam) 
            - gsl_vector_get(fs8obs.f, i);
        gsl_vector_set(diffs8, i, fi);
        for (j=0; j<=i; j++) {
            zj = gsl_vector_get(fs8obs.z, j);
            mij = rec_cov_fs8_with_Hz(zi,zj,sig8,ommh2,gam) 
                + gsl_matrix_get(fs8obs.cov, i,j);
            gsl_matrix_set(Mcov, i, j, mij);
            gsl_matrix_set(Mcov, j, i, mij);
        }
    }

    // inverse Mcov
    gsl_linalg_cholesky_decomp(Mcov);
    gsl_linalg_cholesky_decomp(Mcov);
    
    // chisq
    chisq = 0.0;
    for (i=0; i<fs8obs.n; i++) {
        fi = gsl_vector_get(diffs8, i);
        for (j=0; j<fs8obs.n; j++) {
            fj = gsl_vector_get(diffs8, j);
            mij = gsl_matrix_get(Mcov, i, j);
            chisq += fi*mij*fj;
        }
    }

    // free the gsl vector and matrix
    gsl_vector_free(diffs8);
    gsl_matrix_free(Mcov);
    return -0.5*chisq;
}

double fs8_loglikelihood_Hz_not(double sig8, double ommh2, double gam)
{
    int i, j;
    double chisq, zi, fi, fj, mij;
    gsl_vector *diffs8 = gsl_vector_alloc(fs8obs.n);
    gsl_matrix *Mcov = gsl_matrix_alloc(fs8obs.n, fs8obs.n);

    // initial dif fs8 & cov
    for (i=0; i<fs8obs.n; i++) {
        zi = gsl_vector_get(fs8obs.z, i);
        fi = rec_fs8_with_Hz(zi,sig8,ommh2,gam) 
            - gsl_vector_get(fs8obs.f, i);
        gsl_vector_set(diffs8, i, fi);
        for (j=0; j<=i; j++) {
            mij = gsl_matrix_get(fs8obs.cov, i,j);
            gsl_matrix_set(Mcov, i, j, mij);
            gsl_matrix_set(Mcov, j, i, mij);
        }
    }

    // inverse Mcov
    // gsl_linalg_cholesky_decomp(Mcov);
    // gsl_linalg_cholesky_decomp(Mcov);
    gsl_matrix_inv(Mcov);    
    
    // chisq
    chisq = 0.0;
    for (i=0; i<fs8obs.n; i++) {
        fi = gsl_vector_get(diffs8, i);
        for (j=0; j<fs8obs.n; j++) {
            fj = gsl_vector_get(diffs8, j);
            mij = gsl_matrix_get(Mcov, i, j);
            chisq += fi*mij*fj;
        }
    }

    // free the gsl vector and matrix
    gsl_vector_free(diffs8);
    gsl_matrix_free(Mcov);
    return -0.5*chisq;
}

double fs8_loglikelihood_Ez(double sig8, double omm, double gam)
{
    double ommg = pow(omm, gam);
    double chisq, z, zj, ff, ffs, fi, fj, c00, c01, c10, c11;
    gsl_vector *diffs8 = gsl_vector_alloc(fs8obs.n);
    gsl_matrix *cov = gsl_matrix_alloc(fs8obs.n,fs8obs.n);
    gsl_vector *A = gsl_vector_alloc(fs8obs.n);
    gsl_vector *B = gsl_vector_alloc(fs8obs.n);
    gsl_matrix *K = gsl_matrix_alloc(fs8obs.n,fs8obs.n);
    gsl_matrix *M = gsl_matrix_alloc(fs8obs.n,fs8obs.n);
    gsl_matrix *N = gsl_matrix_alloc(fs8obs.n,fs8obs.n);
    int i, j;
    double pars1[] = {3*gam, 2*gam};
    double pars2[] = {3*gam-1, 2*gam};
    double pars3[] = {0.0, 1.0};
    double pars4[] = {0.0, 3*gam-1, 2*gam+1, 1.0};
    double pars5[] = {3*gam-1, 2*gam+1};

    for (i=0; i<fs8obs.n; i++) {
        z = gsl_vector_get(fs8obs.z, i);
        gsl_vector_set(A, i, func_mu0_over_mu(z, &pars1[0]));
        gsl_vector_set(B, i, integrate_fxmu0_over_fmu(0.0, z, &pars2[0]));
        for (j=0; j<=i; j++) {
            zj = gsl_vector_get(fs8obs.z, j);
            gsl_matrix_set(K, i, j, func_mu0_mu0_covxy_over_mux_muy(z,zj,&pars3[0]));
            ff = gsl_matrix_get(K, i, j);
            gsl_matrix_set(K, j, i, ff);
            
            pars4[0] = zj;
            gsl_matrix_set(M, i, j, integrate_fxmu0covxj_over_fmux(0.0, z, &pars4[0]));
            ff = gsl_matrix_get(M, i, j);
            gsl_matrix_set(M, j, i, ff);
            
            gsl_matrix_set(N, i, j, integrate_fxmu0covxy_over_fmuxy(0.0, z, 0.0, zj, &pars5[0]));
            ff = gsl_matrix_get(N, i, j);
            gsl_matrix_set(N, j, i, ff);
        }
    }

    // rec fs8 mean
    for (i=0; i<fs8obs.n; i++) {
        ff = gsl_vector_get(A,i);
        ffs = gsl_vector_get(B,i);
        fi = sig8*ommg*ff*exp(-ommg * ffs);
        fj = gsl_vector_get(fs8obs.f, i);
        gsl_vector_set(diffs8, i, fi-fj);
    }

    // rec fs8 cov
    for (i=0; i<fs8obs.n; i++) {
        fi = gsl_vector_get(diffs8, i);
        for (j=0; j<fs8obs.n; j++) {
            fj = gsl_vector_get(diffs8, j);
            c00 = gsl_matrix_get(K,i,j);
            c01 = gsl_matrix_get(M,i,j);
            c10 = gsl_matrix_get(M,j,i);
            c11 = gsl_matrix_get(N,i,j);
            ff = fi*fj*pow(2*gam, 2)*( c00 - ommg*c01 
                    - ommg*c10 + ommg*ommg * c11);
            ffs = gsl_matrix_get(fs8obs.cov, i, j);
            gsl_matrix_set(cov, i, j, ff+ffs);
        }
    }

    // inverse cov
    gsl_linalg_cholesky_decomp(cov);
    gsl_linalg_cholesky_decomp(cov);

    // chisq
    chisq = 0.0;
    for (i=0; i<fs8obs.n; i++) {
        fi = gsl_vector_get(diffs8, i);
        for (j=0; j<fs8obs.n; j++) {
            fj = gsl_vector_get(diffs8, j);
            c00 = gsl_matrix_get(cov, i, j);
            chisq += fi*c00*fj;
        }
    }
    
    // free the gsl vector and matrix
    gsl_vector_free(diffs8);
    gsl_vector_free(A);
    gsl_vector_free(B);
    gsl_matrix_free(cov);
    gsl_matrix_free(K);
    gsl_matrix_free(M);
    gsl_matrix_free(N);

    return -0.5*chisq;
}

double fs8_loglikelihood_Ez_not(double sig8, double omm, double gam)
{
    int i, j;
    double chisq, zi, fi, fj, mij;
    gsl_vector *diffs8 = gsl_vector_alloc(fs8obs.n);
    gsl_matrix *Mcov = gsl_matrix_alloc(fs8obs.n, fs8obs.n);

    // initial dif fs8 & cov
    for (i=0; i<fs8obs.n; i++) {
        zi = gsl_vector_get(fs8obs.z, i);
        fi = rec_fs8_with_Ez(zi,sig8,omm,gam) 
            - gsl_vector_get(fs8obs.f, i);
        gsl_vector_set(diffs8, i, fi);
        for (j=0; j<=i; j++) {
            mij = gsl_matrix_get(fs8obs.cov, i,j);
            gsl_matrix_set(Mcov, i, j, mij);
            gsl_matrix_set(Mcov, j, i, mij);
        }
    }

    // inverse Mcov
    // gsl_linalg_cholesky_decomp(Mcov);
    // gsl_linalg_cholesky_decomp(Mcov);
    gsl_matrix_inv(Mcov);
    
    // chisq
    chisq = 0.0;
    for (i=0; i<fs8obs.n; i++) {
        fi = gsl_vector_get(diffs8, i);
        for (j=0; j<fs8obs.n; j++) {
            fj = gsl_vector_get(diffs8, j);
            mij = gsl_matrix_get(Mcov, i, j);
            chisq += fi*mij*fj;
        }
    }

    // free the gsl vector and matrix
    gsl_vector_free(diffs8);
    gsl_matrix_free(Mcov);
    return -0.5*chisq;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
// free fs8
int free_fs8obs()
{
    gsl_vector_free(fs8obs.z);
    gsl_vector_free(fs8obs.f);
    gsl_matrix_free(fs8obs.cov);
    return _SUCCESS_;
}

