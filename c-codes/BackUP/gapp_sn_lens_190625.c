/* ===============================================
 * File Name: gapp_sn_lens.c
 * Author: ekli
 * Mail: ekli_091@mail.dlut.edu.cn  
 * Created Time: 2019-06-22 16:33:18
 * ===============================================
 */

#include "common.h"
#include "gapp.h"

typedef struct {
    int n;
    gsl_vector *zl;
    gsl_vector *zs;
    gsl_vector *dl;
    gsl_vector *ds;
    gsl_vector *obs;
    gsl_vector *obss;
    gsl_vector *sig;
    gsl_matrix *cll;
    gsl_matrix *cls;
    gsl_matrix *csl;
    gsl_matrix *css;
} REC_DZ;

REC_DZ * sl;
REC_DZ * tdsl;

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/*               from mu reconstruct dz             */
/* ------------------------------------------------ */
/* Name: rec_distance_noM                           */
/* Input: z = exp(x)                                */
/* Output: dz = 1/(1+z) exp(ln(10)/5 (mu))          */
/*         need x = log(x) to reconst               */
/*         without M_B,                             */
/*         real is dz/exp(ln(10)/5 M_B)             */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
double rec_distance_noM(double z)
{
    double x;
    x = log(z);
    return 1/(1+z)*exp(log(10)/5 *rec_mu(x));
}

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/*        covariance of reconstructed dz            */
/* ------------------------------------------------ */
/* Name: rec_cov_distance_noM                       */
/* Input: zi, zj                                    */
/* Output: cov(d(zi), d(zj))                        */
/*         real is cov / (exp(ln(10)/5 M_B)^2)      */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
double rec_cov_distance_noM(double zi, double zj)
{
    double xi, xj;
    xi = log(zi);
    xj = log(zj);
    return pow(log(10.0)/5.0, 2.0) * rec_distance_noM(zi) 
        * rec_distance_noM(zj) * rec_covariance(xi, xj);
}


/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/*               from mu reconstruct dz             */
/* ------------------------------------------------ */
/* Name: rec_distance                               */
/* Input: z = exp(x)                                */
/* Output: dz = 1/(1+z) exp(ln(10)/5 (mu-MB))       */
/*         need x = log(x) to reconst               */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
double rec_distance(double z, double MB)
{
    double x;
    x = log(z);
    return 1/(1+z)*exp(log(10)/5 *(rec_mu(x) -MB));
}

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/*        covariance of reconstructed dz            */
/* ------------------------------------------------ */
/* Name: rec_cov_distance                           */
/* Input: zi, zj, MB                                */
/* Output: cov(d(zi), d(zj))                        */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
double rec_cov_distance(double zi, double zj, double MB)
{
    double xi, xj;
    xi = log(zi);
    xj = log(zj);
    return pow(log(10.0)/5.0, 2.0) * rec_distance(zi, MB) 
        * rec_distance(zj, MB) * rec_covariance(xi, xj);
}

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/*       reconstructed dls                          */
/* ------------------------------------------------ */
/* Name: rec_distance_ls                            */
/* Input: zl, zs, MB, omk                           */
// Output: dls = ds \sqrt{1+omk dl^2} - dl \sqrt{1+omk ds^2}
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
double rec_distance_ls(
        double zl, 
        double zs, 
        double MB, 
        double omk)
{
    double dl, ds;
    dl = rec_distance(zl, MB);
    ds = rec_distance(zs, MB);
    return ds*sqrt(1+omk*pow(dl, 2.0)) - dl*sqrt(1+omk*pow(ds, 2.0));
}

// \partial dls/ \partial dl
double rec_pd_dls_over_dl(
        double zl,
        double zs,
        double MB,
        double omk
        )
{
    double dl, ds;
    dl = rec_distance(zl, MB);
    ds = rec_distance(zs, MB);
    return omk*dl*ds/sqrt(1+omk*pow(dl,2)) - sqrt(1+omk*pow(ds, 2));
}

// \partial dls/ \partial ds
double rec_pd_dls_over_ds(
        double zl,
        double zs,
        double MB,
        double omk
        )
{
    double dl, ds;
    dl = rec_distance(zl, MB);
    ds = rec_distance(zs, MB);
    return -omk*dl*ds/sqrt(1+omk*pow(ds,2)) + sqrt(1+omk*pow(dl, 2));
}

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/*            rec Einstein rings radius             */
/* ------------------------------------------------ */
/* Name: rec_Einstein_rings
 * Input: zl, zs, MB, omk
 * Output: D = dls/ds
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
double rec_Einstein_rings(
        double zl, 
        double zs,
        double MB,
        double omk)
{
    double dl, ds;
    dl = rec_distance(zl, MB);
    ds = rec_distance(zs, MB);
    return sqrt(1+omk*pow(dl, 2.0)) - dl/ds*sqrt(1+omk*pow(ds, 2.0));
}

// \partial D/\partial dl
double rec_pd_D_over_dl(
        double zl,
        double zs,
        double MB,
        double omk)
{
    double dl, ds;
    dl = rec_distance(zl, MB);
    ds = rec_distance(zs, MB);
    return omk*dl/sqrt(1+omk*pow(dl,2)) - sqrt(1+omk*pow(ds, 2))/ds;
}

// \partial D/\partial ds
double rec_pd_D_over_ds(
        double zl,
        double zs,
        double MB,
        double omk)
{
    double dl, ds;
    dl = rec_distance(zl, MB);
    ds = rec_distance(zs, MB);
    return dl/pow(ds,2)/sqrt(1+omk*pow(ds,2));
}

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/*           initial strong lensing data            */
/* ------------------------------------------------ */
/* Name: initial_SL
 * Input: N, filename
 * */
int initial_SL(int N, char * filename)
{
    int i, j;
    FILE * fp=NULL;
    double zl, zs, tE, stE, sig, ssig, Dobs, sig_D;
    double xl, xs, yl, ys, cij;

    // initial zl, zs and dl, ds
    sl = (REC_DZ *)malloc(sizeof(REC_DZ));
    sl->n = N;
    sl->zl  = gsl_vector_calloc(N);
    sl->zs  = gsl_vector_calloc(N);
    sl->dl  = gsl_vector_calloc(N);
    sl->ds  = gsl_vector_calloc(N);
    sl->obs = gsl_vector_calloc(N);
    sl->obss= gsl_vector_calloc(N);
    sl->sig = gsl_vector_calloc(N);
    sl->cll = gsl_matrix_calloc(N, N);
    sl->cls = gsl_matrix_calloc(N, N);
    sl->csl = gsl_matrix_calloc(N, N);
    sl->css = gsl_matrix_calloc(N, N);

    // open file
    fp = fopen(filename, "r");

    // read in data
    for (i=0; i<N; i++) {
        fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf\n",
                &zl, &zs, &tE, &stE, &sig, &ssig, &Dobs, &sig_D);
        gsl_vector_set(sl->zl, i, zl);
        gsl_vector_set(sl->zs, i, zs);
        gsl_vector_set(sl->obs, i, Dobs);
        gsl_vector_set(sl->sig, i, sig_D);
    }

    // initial dl ds cov
    for (i=0; i<N; i++) {
        xl = log(sl->zl->data[i]);
        xs = log(sl->zs->data[i]);
        gsl_vector_set(sl->dl, i, rec_distance_noM(xl));
        gsl_vector_set(sl->ds, i, rec_distance_noM(xs));
        for (j=0; j<=i; j++) {
            yl = log(sl->zl->data[j]);
            ys = log(sl->zs->data[j]);
            cij = rec_cov_distance_noM(xl, yl);
            gsl_matrix_set(sl->cll, i, j, cij);
            gsl_matrix_set(sl->cll, j, i, cij);
            cij = rec_cov_distance_noM(xs, ys);
            gsl_matrix_set(sl->css, i, j, cij);
            gsl_matrix_set(sl->css, j, i, cij);
        }
        for (j=0; j<N; j++) {
            yl = log(sl->zl->data[j]);
            cij = rec_cov_distance_noM(xs, yl);
            gsl_matrix_set(sl->csl, i, j, cij);
            gsl_matrix_set(sl->cls, j, i, cij);
        }
    }

    // close file
    fclose(fp);

    return _SUCCESS_;
}

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/*      initial time delay strong lensing data      */
/* ------------------------------------------------ */
/* Name: initial_TDSL
 * Input: N, filename
 * */

int initial_TDSL(int N, char * filename)
{
    int i, j;
    FILE * fp=NULL;
    double zl, zs, sig, Dobs, sig_D;
    double xl, xs, yl, ys, cij;

    // initial zl, zs and dl, ds
    tdsl = (REC_DZ *)malloc(sizeof(REC_DZ));
    tdsl->n = N;
    tdsl->zl  = gsl_vector_calloc(N);
    tdsl->zs  = gsl_vector_calloc(N);
    tdsl->dl  = gsl_vector_calloc(N);
    tdsl->ds  = gsl_vector_calloc(N);
    tdsl->obs = gsl_vector_calloc(N);
    tdsl->obss= gsl_vector_calloc(N);
    tdsl->sig = gsl_vector_calloc(N);
    tdsl->cll = gsl_matrix_calloc(N, N);
    tdsl->cls = gsl_matrix_calloc(N, N);
    tdsl->csl = gsl_matrix_calloc(N, N);
    tdsl->css = gsl_matrix_calloc(N, N);

    // open file
    fp = fopen(filename, "r");

    // read in data
    for (i=0; i<N; i++) {
        fscanf(fp, "%lf %lf %lf %lf %lf\n",
                &zl, &zs, &Dobs, &sig_D, &sig);
        gsl_vector_set(tdsl->zl, i, zl);
        gsl_vector_set(tdsl->zs, i, zs);
        gsl_vector_set(tdsl->obs, i, Dobs); // mu_D
        gsl_vector_set(tdsl->sig, i, sig_D); // sig_D
        gsl_vector_set(tdsl->obss, i, sig); // lambda_D
    }

    // initial dl ds cov
    for (i=0; i<N; i++) {
        xl = log(tdsl->zl->data[i]);
        xs = log(tdsl->zs->data[i]);
        gsl_vector_set(tdsl->dl, i, rec_distance_noM(xl));
        gsl_vector_set(tdsl->ds, i, rec_distance_noM(xs));
        for (j=0; j<=i; j++) {
            yl = log(tdsl->zl->data[j]);
            ys = log(tdsl->zs->data[j]);
            cij = rec_cov_distance_noM(xl, yl);
            gsl_matrix_set(tdsl->cll, i, j, cij);
            gsl_matrix_set(tdsl->cll, j, i, cij);
            cij = rec_cov_distance_noM(xs, ys);
            gsl_matrix_set(tdsl->css, i, j, cij);
            gsl_matrix_set(tdsl->css, j, i, cij);
        }
        for (j=0; j<N; j++) {
            yl = log(tdsl->zl->data[j]);
            cij = rec_cov_distance_noM(xs, yl);
            gsl_matrix_set(tdsl->csl, i, j, cij);
            gsl_matrix_set(tdsl->cls, j, i, cij);
        }
    }

    // close file
    fclose(fp);

    return _SUCCESS_;

}

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/*  free the malloc struct & calloc array  */
/* --------------------------------------- */
int free_SL_TDSL(void)
{
    gsl_vector_free(sl->zl);
    gsl_vector_free(sl->zs);
    gsl_vector_free(sl->dl);
    gsl_vector_free(sl->ds);
    gsl_vector_free(sl->obs);
    gsl_vector_free(sl->obss);
    gsl_vector_free(sl->sig);
    gsl_matrix_free(sl->cll);
    gsl_matrix_free(sl->cls);
    gsl_matrix_free(sl->csl);
    gsl_matrix_free(sl->css);
    free(sl);
    
    gsl_vector_free(tdsl->zl);
    gsl_vector_free(tdsl->zs);
    gsl_vector_free(tdsl->dl);
    gsl_vector_free(tdsl->ds);
    gsl_vector_free(tdsl->obs);
    gsl_vector_free(tdsl->obss);
    gsl_vector_free(tdsl->sig);
    gsl_matrix_free(tdsl->cll);
    gsl_matrix_free(tdsl->cls);
    gsl_matrix_free(tdsl->csl);
    gsl_matrix_free(tdsl->css);
    free(tdsl);

    return _SUCCESS_;
}

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* strong lensing data -> theortical distance ratio */
/*         the correspondence loglikelihood         */
/* ------------------------------------------------ */
/* Name: SL_loglikelihood
 * Input: MB, omk, fE
 * Output: -0.5*chisq
 * */
double SL_loglikelihood(double MB, double omk, double fE)
{
    int i, j;
    double dl, ds, dls, fMB, FMB, cij; //dli, dlj, dsi, dsj, cij;
    double logDet, chisq;
    gsl_vector *temp = gsl_vector_calloc(sl->n); //temp
    gsl_vector *D_Dth = gsl_vector_calloc(sl->n); //D_th - D_obs
    gsl_matrix *COV = gsl_matrix_calloc(sl->n, sl->n); // Cov matrix
    double *pdD_dl = (double *)calloc(sl->n, sizeof(double));
    double *pdD_ds = (double *)calloc(sl->n, sizeof(double));

    // exp(ln10/5 M_B)
    fMB = log(10)/5.0 * MB;
    FMB = exp(fMB);

    // initial dls/ds 
    for (i=0; i<sl->n; i++) {
        ds = gsl_vector_get(sl->ds, i)/FMB;
        dl = gsl_vector_get(sl->dl, i)/FMB;
        dls = ds * sqrt(1+omk*dl*dl) - dl * sqrt(1+omk*ds*ds);
        gsl_vector_set(D_Dth, i, dls/ds - sl->obs->data[i]);
        pdD_dl[i] = omk*dl/sqrt(1+omk*dl*dl) - sqrt(1+omk*ds*ds)/ds;
        pdD_ds[i] = dl/ds/ds/sqrt(1+omk*ds*ds);
    }

    // initial total cov
    for (i=0; i<sl->n; i++) {
        //dsi = gsl_vector_get(sl->ds, i)/FMB;
        //dli = gsl_vector_get(sl->dl, i)/FMB;
        for (j=0; j<=i; j++) {
            //dsj = gsl_vector_get(sl->ds, j)/FMB;
            //dlj = gsl_vector_get(sl->dl, j)/FMB;
            cij = pdD_dl[i] * pdD_dl[j] * gsl_matrix_get(sl->cll, i, j) 
                + pdD_dl[i] * pdD_ds[j] * gsl_matrix_get(sl->cls, i, j)
                + pdD_ds[i] * pdD_dl[j] * gsl_matrix_get(sl->csl, i, j)
                + pdD_ds[i] * pdD_ds[j] * gsl_matrix_get(sl->css, i, j);
            if(i == j){
                cij += sl->sig->data[i]* sl->sig->data[i];
            }
            gsl_matrix_set(COV, i, j, cij);
            gsl_matrix_set(COV, j, i, cij);
        }
    }

    // inverse total cov
    gsl_linalg_cholesky_decomp(COV);
    logDet = 0;
    for (i=0; i<sl->n; i++) {
        logDet += log(gsl_matrix_get(COV, i, j));
    }
    logDet = 2.0*logDet;

    gsl_linalg_cholesky_invert(COV);

    gsl_blas_dgemv(CblasNoTrans,
            1.0, COV, D_Dth,
            0.0, temp);
    gsl_blas_ddot(D_Dth, temp, &chisq);

    chisq += logDet + sl->n * log(2*_PI_);
    return -0.5*chisq;
}

