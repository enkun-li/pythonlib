/* ===============================================
 * File Name: stronglens.c
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-06-27 20:41:09
 * =============================================== 
 */

#include "common.h"
#include "cosmo_bk.h"
#include "cosmography.h"

// typedef data struct
typedef struct {
    int n;
    double *zl;
    double *zs;
    double *obs;
    double *sig;
} DATA;

DATA *sl;

// define a flag to determine whether or not initial_SL_data
int flag = 0;

// initial_SL_data
int initial_SL_data(int N, char *filename)
{
    int i;
    FILE *fp=NULL;
    double tE, stE, sig, ssig;

    // initial data array
    sl = (DATA *)malloc(sizeof(DATA));
    sl->n = N;
    sl->zl = (double *)calloc(N, sizeof(double));
    sl->zs = (double *)calloc(N, sizeof(double));
    sl->obs = (double *)calloc(N, sizeof(double));
    sl->sig = (double *)calloc(N, sizeof(double));

    // initial sl data
    fp = fopen(filename, "r");

    // read in data
    for (i=0; i<N; i++) {
        fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf\n",
                &sl->zl[i], &sl->zs[i], &tE, &stE, &sig, &ssig,
                &sl->obs[i], &sl->sig[i]);
    }
    return _SUCCESS_;
}

// sl loglikelihood
double sl_loglikelihood(double (*func)(double, double), double fE)
{
    int i;
    double difobs;
    double chisq = 0;

    // get th - obs and chisq
    for (i=0; i<sl->n; i++) {
        difobs = pow(fE,2.0) * (*func)(sl->zl[i], sl->zs[i]) - sl->obs[i];
        chisq += pow(difobs, 2)/pow(sl->sig[i], 2);
    }

    // check whether is nan
    if(gsl_isnan(chisq))
        return -3.0e30;    
    return -0.5*chisq;
}

// free sl
int free_stronglens(void)
{
    free(sl);
    return _SUCCESS_;
}

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* output LCDM model loglikelihood with Stong Lens data */
/* ---------------------------------------------------- */
double lcdm_stronglens(double h0, double omm, double omk, 
        double fE)
{
    //int N = 205;
    //char * file_SL = "/home/ekli/myworks/cosmodata/strong_lens_obs.txt";
    //int N = 173;
    //char * file_SL = "/home/ekli/myworks/cosmodata/strong_lens_obs_real.txt";
    int N = 167;
    char * file_SL = "/home/ekli/myworks/cosmodata/strong_lens_obs_SN.txt";

    if(flag == 0){
        initial_SL_data(N, file_SL);
        flag = 1;
    }

    // initial LCDM model
    initial_lcdm(omm,omk,h0);

    return sl_loglikelihood(&lcdm_Dobs_SL, fE);
}

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* output cosmography model loglikelihood with Stong Lens data */
/* ----------------------------------------------------------- */
double cosmography_stronglens(double omk, double a1, double a2,
        double fE)
{
    //int N = 205;
    //char * file_SL = "/home/ekli/myworks/cosmodata/strong_lens_obs.txt";
    //int N = 173;
    //char * file_SL = "/home/ekli/myworks/cosmodata/strong_lens_obs_real.txt";
    int N = 167;
    char * file_SL = "/home/ekli/myworks/cosmodata/strong_lens_obs_SN.txt";

    if(flag == 0){
        initial_SL_data(N, file_SL);
        flag = 1;
    }

    // initial cosmography model
    initial_cosmography(omk,a1,a2);

    return sl_loglikelihood(&cosmography_Dobs_SL, fE);
}


