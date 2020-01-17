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
    double *dl;
    double *ds;
    double *obs;
    double *obss;
    double *sig;
    double *cll;
    double *cls;
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
    return 1/(1+z)*pow(10.0, 0.2*rec_mu(x));
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
    return pow(0.2*log(10.0), 2.0) * rec_distance_noM(zi) 
        * rec_distance_noM(zj) * rec_covariance(xi, xj);
}

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/*  output mean and cov of dl & ds @ zl & zs of SL  */
/* ------------------------------------------------ */
/* Name: output_mean_cov_of_SL
 * Input: N, filename
 * Output: file_of_dl_ds
 * */
int output_mean_cov_of_SL(
        int N, 
        char *filename, 
        char *out_mean, 
        char *out_cov)
{
    FILE *fp = NULL;
    FILE *fpm = NULL;
    FILE *fpc = NULL;
    int i, j;
    double zl[N], zs[N];
    double tE, stE, sig0, ssig0, Dobs, sDobs;
    double cll[N*N], cls[N*N], csl[N*N], css[N*N];

    // open file
    fp = fopen(filename, "r");
    fpm = fopen(out_mean, "w+");
    fpc = fopen(out_cov, "w+");

    // read in data & write to file
    for (i=0; i<N; i++) {
        fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf\n",
                &zl[i], &zs[i], &tE, &stE, &sig0, &ssig0, 
                &Dobs, &sDobs);
        fprintf(fpm, "%8.3f %8.3f ", zl[i], zs[i]);
        fprintf(fpm, "%20.6f %20.6f ", rec_distance_noM(zl[i]), rec_distance_noM(zs[i]));
        fprintf(fpm, "%12.6f %12.6f ", Dobs, sDobs);
        fprintf(fpm, "\n");
    }

    // cov ll, ss
    for (i=0; i<N; i++) {
# pragma omp parallel for
        for (j=0; j<=i; j++) {
            cll[i*N+j] = rec_cov_distance_noM(zl[i], zl[j]);
            css[i*N+j] = rec_cov_distance_noM(zs[i], zs[j]);
            cll[j*N+i] = cll[i*N+j];
            css[j*N+i] = css[i*N+j];
        }
    }
    // cov ls, sl
    for (i=0; i<N; i++){
# pragma omp parallel for
        for (j=0; j<N; j++) {
            cls[i*N+j] = rec_cov_distance_noM(zl[i], zs[j]);
            csl[j*N+i] = cls[i*N+j];
        }
    }

    // write
    for (i=0; i<N*N; i++) {
        fprintf(fpc, "%20.6f %20.6f %20.6f %20.6f\n", 
                cll[i], cls[i], csl[i], css[i]);
    }

    // close file
    fclose(fp);
    fclose(fpm);
    fclose(fpc);

    return _SUCCESS_;
}

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* output mean and cov of dl & ds @ zl & zs of TDSL */
/* ------------------------------------------------ */
/* Name: output_mean_cov_of_TDSL
 * Input: N, filename
 * Output: file_of_dl_ds
 * */
int output_mean_cov_of_TDSL(int N, char *filename)
{
    FILE *fp = NULL;
    FILE *fpm = NULL;
    FILE *fpc = NULL;
    int i, j;
    double zl[N], zs[N];
    double mu, sig, lam;
    double cll[N*N], cls[N*N], csl[N*N], css[N*N];

    // open file
    fp = fopen(filename, "r");
    fpm = fopen("TDSL_mean.txt", "w+");
    fpc = fopen("TDSL_cov.txt", "w+");

    // read in data & write to file
    for (i=0; i<N; i++) {
        fscanf(fp, "%lf %lf %lf %lf %lf\n",
                &zl[i], &zs[i], &mu, &sig, &lam);
        fprintf(fpm, "%8.3f %8.3f ", zl[i], zs[i]);
        fprintf(fpm, "%20.6f %20.6f ", rec_distance_noM(zl[i]), rec_distance_noM(zs[i]));
        fprintf(fpm, "%12.6f %12.6f %12.6f", mu, sig, lam);
        fprintf(fpm, "\n");
    }

    // cov ll, ss
    for (i=0; i<N; i++) {
# pragma omp parallel for
        for (j=0; j<=i; j++) {
            cll[i*N+j] = rec_cov_distance_noM(zl[i], zl[j]);
            css[i*N+j] = rec_cov_distance_noM(zs[i], zs[j]);
            cll[j*N+i] = cll[i*N+j];
            css[j*N+i] = css[i*N+j];
        }
    }
    // cov ls, sl
    for (i=0; i<N; i++){
# pragma omp parallel for
        for (j=0; j<N; j++) {
            cls[i*N+j] = rec_cov_distance_noM(zl[i], zs[j]);
            csl[j*N+i] = cls[i*N+j];
        }
    }

    // write
    for (i=0; i<N*N; i++) {
        fprintf(fpc, "%20.6f %20.6f %20.6f %20.6f\n", 
                cll[i], cls[i], csl[i], css[i]);
    }

    // close file
    fclose(fp);
    fclose(fpm);
    fclose(fpc);

    return _SUCCESS_;
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
    double *zl = (double *)calloc(N, sizeof(double));
    double *zs = (double *)calloc(N, sizeof(double));


    // initial zl, zs and dl, ds
    sl = (REC_DZ *)malloc(sizeof(REC_DZ));
    sl->n = N;
    sl->dl  = (double *)calloc(N, sizeof(double));
    sl->ds  = (double *)calloc(N, sizeof(double));
    sl->obs = (double *)calloc(N, sizeof(double));
    sl->obss= (double *)calloc(N, sizeof(double));
    sl->sig = (double *)calloc(N, sizeof(double));
    sl->cll = (double *)calloc(N*N, sizeof(double));
    sl->cls = (double *)calloc(N*N, sizeof(double));
    sl->css = (double *)calloc(N*N, sizeof(double));

    // open file
    fp = fopen(filename, "r");

    // read in data
    for (i=0; i<N; i++) {
        fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf\n",
                &zl[i], &zs[i], &tE, &stE, &sig, &ssig, 
                &sl->obs[i], &sl->sig[i]);
    }

    // close file
    fclose(fp);

    // initial dl ds cov
    for (i=0; i<N; i++) {
        sl->dl[i] = rec_distance_noM(zl[i]);
        sl->ds[i] = rec_distance_noM(zs[i]);
        for (j=0; j<=i; j++) {
            sl->cll[i*N+j] = rec_cov_distance_noM(zl[i], zl[j]);
            sl->cll[j*N+i] = sl->cll[i*N+j];
            sl->css[i*N+j] = rec_cov_distance_noM(zs[i], zs[j]);
            sl->css[j*N+i] = sl->css[i*N+j];
        }
        for (j=0; j<N; j++) {
            sl->cls[i*N+j] = rec_cov_distance_noM(zl[i], zs[j]);
        }
    }

    // write data
    fp = fopen("./sl_dl_ds.txt", "w+");
    for (i=0; i<sl->n; i++){
        fprintf(fp, "%20.8e %20.8e\n", sl->dl[i], sl->ds[i]);
    }
    // close file
    fclose(fp);

    // write data
    fp = fopen("./sl_cov_dl_ds.txt", "w+");
    for (i=0; i<sl->n*sl->n; i++){
        fprintf(fp, "%20.8e %20.8e %20.8e\n", sl->cll[i], sl->cls[i], sl->css[i]);
    }
    // close file
    fclose(fp);
    free(zl);
    free(zs);

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
    double cij;
    double *zl  = (double *)calloc(N, sizeof(double));
    double *zs  = (double *)calloc(N, sizeof(double));

    // initial zl, zs and dl, ds
    tdsl = (REC_DZ *)malloc(sizeof(REC_DZ));
    tdsl->n = N;
    tdsl->dl  = (double *)calloc(N, sizeof(double));
    tdsl->ds  = (double *)calloc(N, sizeof(double));
    tdsl->obs = (double *)calloc(N, sizeof(double));
    tdsl->obss= (double *)calloc(N, sizeof(double));
    tdsl->sig = (double *)calloc(N, sizeof(double));
    tdsl->cll = (double *)calloc(N*N, sizeof(double));
    tdsl->cls = (double *)calloc(N*N, sizeof(double));
    tdsl->css = (double *)calloc(N*N, sizeof(double));

    // open file
    fp = fopen(filename, "r");

    // read in data
    for (i=0; i<N; i++) {
        fscanf(fp, "%lf %lf %lf %lf %lf\n",
                &zl[i], &zs[i], &tdsl->obs[i], &tdsl->sig[i], &tdsl->obss[i]);
    }

    // initial dl ds cov
    for (i=0; i<N; i++) {
        tdsl->dl[i] = rec_distance_noM(zl[i]);
        tdsl->ds[i] = rec_distance_noM(zs[i]);
        for (j=0; j<=i; j++) {
            cij = rec_cov_distance_noM(zl[i], zl[j]);
            tdsl->cll[i*N+j] = cij;
            tdsl->cll[j*N+i] = cij;
            cij = rec_cov_distance_noM(zs[i], zs[j]);
            tdsl->css[i*N+j] = cij;
            tdsl->css[j*N+i] = cij;
        }
        for (j=0; j<N; j++) {
            cij = rec_cov_distance_noM(zl[i], zs[j]);
            tdsl->cls[i*N+j] = cij;
        }
    }

    // close file
    fclose(fp);

    free(zl);
    free(zs);

    return _SUCCESS_;

}

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/*  free the malloc struct & calloc array  */
/* --------------------------------------- */
int free_SL_TDSL(void)
{
    free(sl->dl);
    free(sl->ds);
    free(sl->obs);
    free(sl->obss);
    free(sl->sig);
    free(sl->cll);
    free(sl->cls);
    free(sl->css);
    free(sl);
    
    free(tdsl->dl);
    free(tdsl->ds);
    free(tdsl->obs);
    free(tdsl->obss);
    free(tdsl->sig);
    free(tdsl->cll);
    free(tdsl->cls);
    free(tdsl->css);
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
        ds = sl->ds[i]/FMB;
        dl = sl->dl[i]/FMB;
        dls = ds * sqrt(1+omk*dl*dl) - dl * sqrt(1+omk*ds*ds);
        gsl_vector_set(D_Dth, i, dls/ds - sl->obs[i]);
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
            cij = pdD_dl[i] * pdD_dl[j] * sl->cll[i*sl->n+j]
                + pdD_dl[i] * pdD_ds[j] * sl->cls[i*sl->n+j]
                + pdD_ds[i] * pdD_dl[j] * sl->cls[j*sl->n+i]
                + pdD_ds[i] * pdD_ds[j] * sl->css[i*sl->n+j];
            if(i == j){
                cij += sl->sig[i]* sl->sig[i];
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
    
    free(pdD_dl);
    free(pdD_ds);
    gsl_vector_free(temp);
    gsl_vector_free(D_Dth);
    gsl_matrix_free(COV);

    return -0.5*chisq;
}

