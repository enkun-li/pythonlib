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
    double *zl;
    double *zs;
    double *dl;
    double *ds;
    double *obs;
    double *sig;
    double *cll;
    double *cls;
    double *csl;
    double *css;
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
    double tE, stE, sig, ssig;
    double xl, xs, yl, ys;

    // initial zl, zs and dl, ds
    sl = (REC_DZ *)malloc(sizeof(REC_DZ));
    sl->n = N;
    sl->zl = (double *)calloc(N, sizeof(double));
    sl->zs = (double *)calloc(N, sizeof(double));
    sl->dl = (double *)calloc(N, sizeof(double));
    sl->ds = (double *)calloc(N, sizeof(double));
    sl->obs = (double *)calloc(N, sizeof(double));
    sl->sig = (double *)calloc(N, sizeof(double));
    sl->cll = (double *)calloc(N*N, sizeof(double));
    sl->cls = (double *)calloc(N*N, sizeof(double));
    sl->csl = (double *)calloc(N*N, sizeof(double));
    sl->css = (double *)calloc(N*N, sizeof(double));

    // open file
    fp = fopen(filename, "r");

    // read in data
    for (i=0; i<N; i++) {
        fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf\n",
                &sl->zl[i], &sl->zs[i], &tE, &stE, &sig, &ssig,
                &sl->obs[i], &sl->sig[i]);
    }

    // initial dl ds cov
    for (i=0; i<N; i++) {
        xl = log(sl->zl[i]);
        xs = log(sl->zs[i]);
        sl->dl[i] = rec_distance_noM(xl);
        sl->ds[i] = rec_distance_noM(xs);
        for (j=0; j<N; j++) {
            yl = log(sl->zl[j]);
            ys = log(sl->zs[j]);
            sl->cll[i*N+j] = rec_cov_distance_noM(xl, yl);
            sl->css[i*N+j] = rec_cov_distance_noM(xs, ys);
            sl->csl[i*N+j] = rec_cov_distance_noM(xs, yl);
            sl->cls[i*N+j] = rec_cov_distance_noM(xl, ys);
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

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/*  free the malloc struct & calloc array  */
/* --------------------------------------- */
int free_SL_TDSL(void)
{
    free(sl->zl);
    free(sl->zs);
    free(sl->dl);
    free(sl->ds);
    free(sl->obs);
    free(sl->sig);
    free(sl->cll);
    free(sl->cls);
    free(sl->css);
    free(sl);
    free(tdsl->zl);
    free(tdsl->zs);
    free(tdsl->dl);
    free(tdsl->ds);
    free(tdsl->obs);
    free(tdsl->sig);
    free(tdsl->cll);
    free(tdsl->cls);
    free(tdsl->css);
    free(tdsl);

    return _SUCCESS_;
}
