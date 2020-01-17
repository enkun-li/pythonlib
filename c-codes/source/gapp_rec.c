/* ===============================================
 * File Name: gapp_rec.c
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-06-26 19:54:08
 * =============================================== 
 */

#include "common.h"
#include "gapp.h"
#include "utils.h"
#include "cosmo_bk.h"

// rec mu0/mu
double rec_mu0_over_mu(double x, double * args)
{
    return pow(rec_mu(0.0)/rec_mu(x), *args);
}

// rec cov(x,y)/mu(x)/mu(y)
double rec_covxy_over_muxy(double x, double y)
{
    return rec_covariance(x,y)/rec_mu(x)/rec_mu(y);
}

// rec cov mu0/mu
double rec_cov_mu0_over_mu(double x, double y)
{
    double args[] = {1};
    return rec_mu0_over_mu(x, &args[0]) *
        rec_mu0_over_mu(y, &args[0]) *
        (rec_covxy_over_muxy(0.0, 0.0)
        - rec_covxy_over_muxy(0.0,y)
        - rec_covxy_over_muxy(x,0.0)
        + rec_covxy_over_muxy(x,y) );
}

// rec cov(Ei,Ej)/(Ei*Ej)^arg
double rec_covmu0mu_over_fmu0mu(double x, double y, double * args)
{
    return rec_cov_mu0_over_mu(x,y)/
        (rec_mu0_over_mu(x, &args[0]) * rec_mu0_over_mu(y, &args[0]));
}

// rec int_0^x mu0/mu(x) dx
double int_mu0_over_mu(double x)
{
    double args[] = {1.0};
    return qgausleg(&rec_mu0_over_mu, 0.0, x, &args[0]);
}

// rec cov int_0^x mu0/mu(x)
double cov_int_mu0_over_mu(double x, double y)
{
    double args[] = {0.0};
    return qgausleg2d_fix(&rec_covmu0mu_over_fmu0mu, 0, x, 0, y, &args[0]);
}

/* @@@@@@@@@@@@@@@@@@ */
/* output Pantheon dz */
/* ------------------ */
int return_Pantheon_dz(void)
{
    int N =1048;
    char * filename = "/home/ekli/myworks/cosmodata/sn_full_dat.txt";
    FILE *fp = NULL;
    FILE *fps = NULL;
    int i;
    double z, m, sig, dz, sigdz;

    // open file
    fp = fopen(filename, "r");
    fps = fopen("Pantheon_gapp_dz.txt", "w");

    initial_lcdm(0.3,0.0,70);

    // read write file
    for (i=0; i<N; i++) {
        fscanf(fp, "%lf %lf %lf\n", &z, &m, &sig);
        dz = int_mu0_over_mu(z);
        sigdz = sqrt( cov_int_mu0_over_mu(z, z) );
        fprintf(fps, "%8.3f %12.6f %20.8e %12.6f %12.6f",
                z, dz, sigdz, m, sig);
        fprintf(fps, "%12.6f %20.8e\n",
                //lcdm_como_ang_diam_distance(z)/(_C_/1000)*70,
                5*log10((1+z)*dz) ,//-5*log10(rec_mu(0.0)/(_C_/1000) +25),
                5/log(10.0)*sigdz/dz);
    }

    return _SUCCESS_;
}
