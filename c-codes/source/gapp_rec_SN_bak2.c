/* ===============================================
 * File Name: gapp_rec_SN.c
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-06-28 14:20:55
 * =============================================== 
 */

#include "common.h"
#include "gapp.h"
#include "gapp_sn_lens.h"
#include "utils.h"

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
int sl_use_num = 0;

int sl_n[] = {205, 173, 167, 138, 120};
char *sl_dat[] = {
    "/home/ekli/myworks/cosmodata/strong_lens_obs.txt",
    "/home/ekli/myworks/cosmodata/strong_lens_obs_real.txt",
    "/home/ekli/myworks/cosmodata/strong_lens_obs_SN.txt",
    "/home/ekli/myworks/cosmodata/strong_lens_obs_real_SN.txt",
    "/home/ekli/myworks/cosmodata/strong_lens_obs_real_1.55.txt"
};

char *sl_out[] = {
    "./data/sl_dldscov.txt",
    "./data/sl_dldscov_real.txt",
    "./data/sl_dldscov_SN.txt",
    "./data/sl_dldscov_real_SN.txt",
    "./data/sl_dldscov_real_1.55.txt"
};

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
int flag_TDSL_N = 4;
char * flag_TDSL_dat = "/home/ekli/myworks/cosmodata/time_delay_obs.txt";
char * flag_TDSL_cov = "tdsl_dldscov.txt";

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// type data
typedef struct {
    int n;
    double *dl;
    double *ds;
    double *obs;
    double *sig;
    double *lam;
    double *cll;
    double *cls;
    double *csl;
    double *css;
} DATA;

DATA *sl;
DATA *tdsl;

int flag_SL = 0;
int flag_TDSL = 0;


// dz*10^{0.2*(Mb-25)} = 1/(1+z) 10^{0.2*(m-25)}
double dz_star(double z)
{
    double x = log(z);
    return 1.0/(1.0+z)*pow(10.0, 0.2*(rec_mu(x)-25));
}

// cov dz*10^{0.2*(Mb-25)}
double cov_dz_star(double zi, double zj)
{
    return pow(0.2*log(10), 2.0) * dz_star(zi) * dz_star(zj) *
        rec_covariance(zi, zj);
}

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/*             out put used data                 */
/* --------------------------------------------- */
int initial_SL_usedat(int num)
{
    int N = sl_n[num];
    char * file_dat = sl_dat[num];
    char * file_cov = sl_out[num];

    int i, status;
    FILE *fp = NULL;
    FILE *fps = NULL;
    int nsn = 1048;
    char * filecov = "/home/ekli/myworks/cosmodata/sn_full_cov.txt";
    char * filename = "/home/ekli/myworks/cosmodata/sn_full_lnz_dat.txt";
    double sigf=39.21533362, lenf=12.97267722;
    double zl, zs, tE, stE, sig0, ssig0, obs, sig;
    double dl, ds, cll, cls, csl, css;
    
    // initial pantheon
    if((status = initial_gapp(nsn, filename)) == _SUCCESS_
            && (status = initial_gapp_cov(N, filecov)) == _SUCCESS_
            && (status = setup_gapp(sigf, lenf)) == _SUCCESS_ 
            )
        printf("\nSuccessful initial gapp\n");
    else
        printf("\nCan not setup covariance of !\n");

    // open file
    fp = fopen(file_dat, "r");
    fps = fopen(file_cov, "w");

    // read in means and obs data
    // read in data & write to file
    for (i=0; i<N; i++) {
        fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf\n",
                &zl, &zs, &tE, &stE, &sig0, &ssig0, &obs, &sig);
        dl = dz_star(zl);
        ds = dz_star(zs);
        cll = cov_dz_star(zl,zl);
        cls = cov_dz_star(zl,zs);
        csl = cov_dz_star(zs,zl);
        css = cov_dz_star(zs,zs);
        fprintf(fps, "%12.8f %12.8f %20.8e %20.8e %20.8e %20.8e %12.8f %20.8e\n",
                dl, ds, cll, cls, csl, css, obs, sig);
    }

    // close file
    fclose(fp);
    fclose(fps);

    // free gapp
    free_gapp();

    return 0;

}

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/*                 SL_loglike                     */
/* ---------------------------------------------- */
int initial_SL_loglike(
        int N, 
        char * file_dat, 
        char * file_cov)
{
    int i, status;
    FILE *fp = NULL;
    FILE *fps = NULL;
    int nsn = 1048;
    char * filecov = "/home/ekli/myworks/cosmodata/sn_full_cov.txt";
    char * filename = "/home/ekli/myworks/cosmodata/sn_full_lnz_dat.txt";
    double sigf=39.21533362, lenf=12.97267722;
    //double zl, zs, tE, stE, sig0, ssig0;
    
    // initial pantheon
    if((status = initial_gapp(nsn, filename)) == _SUCCESS_
            && (status = initial_gapp_cov(N, filecov)) == _SUCCESS_
            && (status = setup_gapp(sigf, lenf)) == _SUCCESS_ 
            )
        printf("\nSuccessful initial gapp\n");
    else
        printf("\nCan not setup covariance of !\n");

    
    // initial zl, zs and dl, ds
    sl = (DATA *)malloc(sizeof(DATA));
    sl->n = N;
    sl->dl  = (double *)calloc(N, sizeof(double));
    sl->ds  = (double *)calloc(N, sizeof(double));
    sl->obs = (double *)calloc(N, sizeof(double));
    sl->sig = (double *)calloc(N, sizeof(double));
    sl->lam = (double *)calloc(N, sizeof(double));
    sl->cll = (double *)calloc(N, sizeof(double));
    sl->cls = (double *)calloc(N, sizeof(double));
    sl->csl = (double *)calloc(N, sizeof(double));
    sl->css = (double *)calloc(N, sizeof(double));

    // open file
    fp = fopen(file_dat, "r");
    fps = fopen(file_cov, "r");

    // read in means and obs data
    // read in data & write to file
    //for (i=0; i<N; i++) {
    //    fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf\n",
    //            &zl, &zs, &tE, &stE, &sig0, &ssig0, 
    //            &sl->obs[i], &sl->sig[i]);
    //    sl->dl[i] = dz_star(zl);
    //    sl->ds[i] = dz_star(zs);
    //    sl->cll[i] = cov_dz_star(zl,zl);
    //    sl->cls[i] = cov_dz_star(zl,zs);
    //    sl->csl[i] = cov_dz_star(zs,zl);
    //    sl->css[i] = cov_dz_star(zs,zs);
    //}

    // write cov data
    //for (i=0; i<N; i++) {
    //    fprintf(fps, "%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n",
    //            sl->dl[i], sl->ds[i],
    //            sl->cll[i], sl->cls[i], sl->csl[i], sl->css[i]);
    //}

    //
    for (i=0; i<N; i++) {
        fscanf(fps, "%lf %lf %le %le %le %le %lf %le\n",
                &sl->dl[i], &sl->ds[i], 
                &sl->cll[i], &sl->cls[i], &sl->csl[i], 
                &sl->css[i], &sl->obs[i], &sl->sig[i]);
    }
    // close file
    fclose(fp);
    fclose(fps);

    // free gapp
    free_gapp();

    return 0;

}

int free_SL(void)
{
    free(sl);
    flag_SL = 0;
    return _SUCCESS_;
}

double SL_loglike(
        double MB, 
        double omk, 
        double fE)
{
    int N = sl->n;
    int i;
    double chisq;
    double FMB = pow(10.0, 0.2*(25-MB));
    double FMBsq = pow(FMB,2.0);
    double dl, ds;
    double Dth[N], pdDdl[N], pdDds[N],
           VD[N], Cov[N];

    // give Dth - Dobs vector
//# pragma omp parallel for
    for (i=0; i<N; i++) {
        dl = sl->dl[i]*FMB;
        ds = sl->ds[i]*FMB;
        Dth[i] = sqrt(1+omk*pow(dl, 2.0)) 
            - dl/ds*sqrt(1+omk*pow(ds,2.0));
        
        pdDdl[i] = omk * dl/sqrt(1+omk*pow(dl, 2.0)) 
            - sqrt(1+omk*pow(ds, 2))/ds;
        pdDds[i] = dl/pow(ds, 2.0)/sqrt(1+omk*pow(ds, 2.0));
        
        VD[i] = Dth[i] * pow(fE, 2.0) - sl->obs[i];
    }

    chisq = 0.0;
    for (i=0; i<N; i++) {
        Cov[i] = pdDdl[i] * pdDdl[i] * sl->cll[i]*FMBsq
            + pdDdl[i] * pdDds[i] * sl->cls[i]*FMBsq
            + pdDds[i] * pdDdl[i] * sl->csl[i]*FMBsq
            + pdDds[i] * pdDds[i] * sl->css[i]*FMBsq 
            + pow(sl->sig[i], 2.0);
        chisq += pow(VD[i], 2.0)/Cov[i]; //+ log(Cov[i]);
    }

    if(gsl_isnan(chisq))
        return -3e30;
    return -0.5*chisq;
}

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/*                 TDSL_loglike                  */
/* --------------------------------------------- */
int initial_TDSL_loglike(
        int N, 
        char * file_dat, 
        char * file_cov)
{
    int i, status;
    FILE *fp = NULL;
    FILE *fps = NULL;
    int nsn = 1048;
    char * filecov = "/home/ekli/myworks/cosmodata/sn_full_cov.txt";
    char * filename = "/home/ekli/myworks/cosmodata/sn_full_lnz_dat.txt";
    double sigf=39.21533362, lenf=12.97267722;
    double zl, zs;
    
    // initial pantheon
    if((status = initial_gapp(nsn, filename)) == _SUCCESS_
            && (status = initial_gapp_cov(N, filecov)) == _SUCCESS_
            && (status = setup_gapp(sigf, lenf)) == _SUCCESS_ 
            )
        printf("\nSuccessful initial gapp\n");
    else
        printf("\nCan not setup covariance of !\n");

    
    // initial zl, zs and dl, ds
    tdsl = (DATA *)malloc(sizeof(DATA));
    tdsl->n = N;
    tdsl->dl  = (double *)calloc(N, sizeof(double));
    tdsl->ds  = (double *)calloc(N, sizeof(double));
    tdsl->obs = (double *)calloc(N, sizeof(double));
    tdsl->sig = (double *)calloc(N, sizeof(double));
    tdsl->lam = (double *)calloc(N, sizeof(double));
    tdsl->cll = (double *)calloc(N, sizeof(double));
    tdsl->cls = (double *)calloc(N, sizeof(double));
    tdsl->csl = (double *)calloc(N, sizeof(double));
    tdsl->css = (double *)calloc(N, sizeof(double));

    // open file
    fp = fopen(file_dat, "r");
    fps = fopen(file_cov, "w");

    // read in means and obs data
    // read in data & write to file
    for (i=0; i<N; i++) {
        fscanf(fp, "%lf %lf %lf %lf %lf\n",
                &zl, &zs,
                &tdsl->obs[i], &tdsl->sig[i], &tdsl->lam[i]);
        tdsl->dl[i] = dz_star(zl);
        tdsl->ds[i] = dz_star(zs);
        tdsl->cll[i] = cov_dz_star(zl,zl);
        tdsl->cls[i] = cov_dz_star(zl,zs);
        tdsl->csl[i] = cov_dz_star(zs,zl);
        tdsl->css[i] = cov_dz_star(zs,zs);
    }

    // write cov data
    for (i=0; i<N; i++) {
        fprintf(fps, "%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n",
                tdsl->dl[i], tdsl->ds[i],
                tdsl->cll[i], tdsl->cls[i], tdsl->csl[i], tdsl->css[i]);
    }

    // close file
    fclose(fp);
    fclose(fps);

    // free gapp
    free_gapp();

    //
    printf("\nSuccessful initial_TDSL_loglike\n");

    return 0;

}

int free_TDSL(void)
{
    free(tdsl);
    flag_TDSL = 0;
    return _SUCCESS_;
}

double TDSL_loglike(
        double MB, 
        double omk, 
        double H0)
{
    int N = tdsl->n;
    int i;
    double chisq;
    double FMB = pow(10.0, 0.2*(25-MB));
    double FMBsq = pow(FMB,2.0);
    double dl, ds, dls;
    double Dth[N], pdDdl[N], pdDds[N],
           VD[N], Cov[N], lCov[N];

    // give Dth - Dobs vector
//# pragma omp parallel for
    for (i=0; i<N; i++) {
        dl = tdsl->dl[i]*FMB;
        ds = tdsl->ds[i]*FMB;
        
        dls =  ds*sqrt(1+omk*pow(dl, 2.0))
            - dl*sqrt(1+omk*pow(ds,2.0));
        
        Dth[i] = (_C_/1000.0)/H0 * dl*ds / dls;
        
        pdDdl[i] = (_C_/1000.0)/H0 * pow(ds, 2)/pow(dls, 2) /
            sqrt(1+omk*pow(dl, 2.0));
        pdDds[i] = (_C_/1000.0)/H0 * pow(dl, 2)/pow(dls, 2) /
            sqrt(1+omk*pow(ds, 2.0));
        VD[i] = log(Dth[i] - tdsl->lam[i])  - tdsl->obs[i];
    }

    chisq = 0.0;
    for (i=0; i<N; i++) {
        lCov[i] = pdDdl[i] * pdDdl[i] * tdsl->cll[i]*FMBsq
            + pdDdl[i] * pdDds[i] * tdsl->cls[i]*FMBsq
            + pdDds[i] * pdDdl[i] * tdsl->csl[i]*FMBsq
            + pdDds[i] * pdDds[i] * tdsl->css[i]*FMBsq;
        Cov[i] = lCov[i]/pow(Dth[i] -tdsl->lam[i], 2.0)
            + pow(tdsl->sig[i], 2.0);
        chisq += pow(VD[i], 2.0)/Cov[i] 
            + log(Cov[i] * pow(Dth[i] -tdsl->lam[i], 2.0));
    }

    if(gsl_isnan(chisq))
        return -3e30;
    return -0.5*chisq;
}

/* @@@@@@@@@@@@@@@@@@@@@@ */
/*  output loglike  */
/* ---------------------- */
double return_SL_loglike(
        double MB, 
        double omk,
        double fE)
{
    int flag_SL_N = sl_n[sl_use_num];
    char * flag_SL_dat = sl_dat[sl_use_num];
    char * flag_SL_cov = sl_out[sl_use_num];

    if(flag_SL == 0) {
        initial_SL_loglike(flag_SL_N, flag_SL_dat, flag_SL_cov);
        flag_SL = 1;
    }
    return SL_loglike(MB,omk,fE);
}

// margin over MB
double int_SL_loglike(double MB, double * pars)
{
    return exp( SL_loglike(MB,*pars, *(pars+1)) );
}

double margin_of_MB_SL_loglike(
        double omk,
        double fE)
{
    int flag_SL_N = sl_n[sl_use_num];
    char * flag_SL_dat = sl_dat[sl_use_num];
    char * flag_SL_cov = sl_out[sl_use_num];

    double pars[] = {omk,fE};
    double chisq;

    if(flag_SL == 0) {
        initial_SL_loglike(flag_SL_N,flag_SL_dat,flag_SL_cov);
        flag_SL = 1;
    }

    chisq = log( qgausleg(&int_SL_loglike, 16.0, 32.0, &pars[0]) );
    //chisq = SL_loglike(23.81, omk, fE);
    if (gsl_isinf(chisq) == -1 || gsl_isinf(chisq) == 1)
        return -3e30;
    return chisq;
}

// @@@@@@@@@@@@@@@@@@@@@@
// TDSL loglike
double return_TDSL_loglike(
        double MB, 
        double omk,
        double H0)
{
    if(flag_TDSL == 0) {
        initial_TDSL_loglike(flag_TDSL_N, flag_TDSL_dat, flag_TDSL_cov);
        flag_TDSL = 1;
    }
    return TDSL_loglike(MB,H0,omk);
}

// margin over MB
double int_TDSL_loglike(double MB, double * pars)
{
    return exp( TDSL_loglike(MB,*pars, *(pars+1)) );
}

double margin_of_MB_TDSL_loglike(
        double omk,
        double H0)
{
    double pars[] = {H0,omk};
    double chisq;
      
    if(flag_TDSL == 0) {
        printf("\nInitial TDSL lolike\n");
        initial_TDSL_loglike(flag_TDSL_N, flag_TDSL_dat, flag_TDSL_cov);
        flag_TDSL = 1;
    }

    chisq = log( qgausleg(&int_TDSL_loglike, 16.0, 32.0, &pars[0]) );
    //chisq = TDSL_loglike(23.81, omk, H0);
    if (gsl_isinf(chisq) == -1 || gsl_isinf(chisq) == 1)
        return -3e30;
    return chisq;
}

