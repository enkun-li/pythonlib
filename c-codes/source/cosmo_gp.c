/* ===============================================
 * File Name: cosmo_gp.c
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-07-01 15:11:16
 * =============================================== 
 */

#include "common.h"
#include "gapp.h"
#include "gapp_rec.h"
#include "utils.h"

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/*            general initial defination        */
/* -------------------------------------------- */
int which_ohd = 0; // using which ohd data-set

int nhz[] = {51, 31, 23, 54};
char * file_ohd[] = {
    "/home/ekli/myworks/cosmodata/OHD_51.txt",
    "/home/ekli/myworks/cosmodata/OHD_31_CC.txt",
    "/home/ekli/myworks/cosmodata/OHD_23_BAO.txt",
    "/home/ekli/myworks/cosmodata/OHD_54_CC_BAO_P18.txt"
};
double hyper_sigf[] = {170.75260325, 133.80863893, 220.30220737, 157.85759895};
double hyper_lenf[] = {2.61947963, 1.93610328, 3.6377493, 2.35798007};

double omegak = 0.0;

// define a struct
typedef struct {
    int n; // number
    double * x; // x
    double * y; // f(x)
    double * s; // sigma
    double * cov; // covariance
} DATA;

DATA * sn;

typedef struct {
    double *obs;
    double *cov;
} GP;

GP * gp;
GP * thgp;

// a global parameter to determine wherther to call again
int call_again = 1; // 0: yes; 1 no;

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/*                  initial gapp                */
/* -------------------------------------------- */
int initial_cosmo_gp(int num)
{
    int N = nhz[num];
    char * filename = file_ohd[num];
    double sigf = hyper_sigf[num];
    double lenf = hyper_lenf[num];

    printf("\nInitial GP with: %s\n", filename);
    printf("sigf = % .6f\n", sigf);
    printf("lenf = % .6f\n", lenf);
    
    initial_gapp(N, filename);
    setup_gapp(sigf, lenf);

    printf("\nInitial GP Done\n");

    return _SUCCESS_;
}

int setup_cosmo_gp(double omk, int ohdtype)
{
    omegak = omk;
    which_ohd = ohdtype;
    return _SUCCESS_;
}

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/*          reconst cosmology distance          */
/* -------------------------------------------- */
// reconst hubble
double gp_hz(double z)
{
    return rec_mu(z);
}

double gp_cov_hz(double zi, double zj)
{
    return rec_covariance(zi, zj);
}

// reconst dimensionless comoving distance
// dc = \int_0^z dz/(H/H0)
double gp_dimles_comdistance(double z)
{
    double pars = 1.0;
    return qgausleg(&rec_mu0_over_mu, 0.0, z, &pars);
    //return int_mu0_over_mu(z);
}

double gp_cov_dimles_comdistance(double zi, double zj)
{
    double pars = 0.0;
    return qgausleg2d_fix(&rec_covmu0mu_over_fmu0mu, 0.0, zi, 0.0, zj, &pars);
}

// reconst dimensionless comoving angular diameter distance
double gp_dimles_comangdistance(double z)
{
    if(omegak > 1e-8) {
        return 1/sqrt(omegak) * sinh(sqrt(omegak) * gp_dimles_comdistance(z));
    } else if (omegak < -1e-8) {
        return 1/sqrt(-omegak) * sin(sqrt(-omegak) * gp_dimles_comdistance(z));
    } else {
        return gp_dimles_comdistance(z);
    }
}

double gp_cov_dimles_comangdistance(double zi, double zj)
{
    double csni = gp_dimles_comdistance(zi);
    double csnj = gp_dimles_comdistance(zj);
    if(omegak > 1e-8) {
        return cosh(sqrt(omegak) * csni) * cosh(sqrt(omegak) * csnj) *
            gp_cov_dimles_comdistance(zi,zj);
    } else if (omegak < -1e-8) {
        return cos(sqrt(-omegak) * csni) * cos(sqrt(-omegak) * csnj) *
            gp_cov_dimles_comdistance(zi,zj);
    } else {
        return gp_cov_dimles_comdistance(zi,zj);
    }
}

// reconst comoving distance in [Mpc]
double gp_com_distance(double z)
{
    return (_C_/1000)/rec_mu(0.0) * gp_dimles_comdistance(z);
}

// reconst angular diameter distanc in [Mpc]
double gp_ang_distance(double z)
{
    return (_C_/1000)/rec_mu(0.0)/(1+z) * gp_dimles_comangdistance(z);
}

// reconst lumnisty distance
double gp_lum_distance(double z)
{
    return (_C_/1000)/rec_mu(0.0)*(1+z) * gp_dimles_comangdistance(z);
}

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/*          reconstruct distance modulus        */
/* -------------------------------------------- */
// mu = 5 log10( DL/10pc) = 5 log10(DL/Mpc) +25 
//    = 5log10( dL ) - 5log10(H0/c) +25
// mu = 5 log10( (1+z) dz) - 5 log10(H0/c) +25
double gp_DistanceModulus_star(double z)
{
    return 5*log10((1+z)*gp_dimles_comangdistance(z));
}

double gp_cov_DistanceModulus_star(double zi, double zj)
{
    double dzi = gp_dimles_comangdistance(zi);
    double dzj = gp_dimles_comangdistance(zj);
    return pow(5.0/log(10.0), 2) / dzi /dzj * 
        gp_cov_dimles_comangdistance(zi, zj);
}

double gp_DistanceModulus(double z)
{
    return gp_DistanceModulus_star(z) - 5*log10(rec_mu(0.0)/(_C_/1000)) + 25;
}

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/*               distance sum role              */
/* -------------------------------------------- */
// dls/ds = \sqrt{1+omegak dl^2} -dl/ds \sqrt{1+omegak ds^2}
double gp_DistanceSumRole(double zl, double zs)
{
    double dl = gp_dimles_comangdistance(zl);
    double ds = gp_dimles_comangdistance(zs);

    return sqrt(1+omegak * pow(dl, 2)) - dl/ds * sqrt(1+omegak * pow(ds, 2));
}

double gp_cov_DistanceSumRole(double zli, double zsi, double zlj, double zsj)
{
    double dli = gp_dimles_comangdistance(zli);
    double dsi = gp_dimles_comangdistance(zsi);
    double dlj = gp_dimles_comangdistance(zlj);
    double dsj = gp_dimles_comangdistance(zsj);

    double pDdli = omegak * dli/sqrt(1+omegak*pow(dli, 2.0)) 
            - sqrt(1+omegak*pow(dsi, 2))/dsi;
    double pDdsi = dli/pow(dsi, 2.0)/sqrt(1+omegak*pow(dsi, 2.0));
    double pDdlj = omegak * dlj/sqrt(1+omegak*pow(dlj, 2.0)) 
            - sqrt(1+omegak*pow(dsj, 2))/dsj;
    double pDdsj = dlj/pow(dsj, 2.0)/sqrt(1+omegak*pow(dsj, 2.0));

    return pDdli * pDdlj * gp_cov_dimles_comangdistance(zli, zlj)
        + pDdli * pDdsj * gp_cov_dimles_comangdistance(zli, zsj)
        + pDdsi * pDdlj * gp_cov_dimles_comangdistance(zsi, zlj)
        + pDdsi * pDdsj * gp_cov_dimles_comangdistance(zsi, zsj);
}

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/*          out put comoving distance           */
/* -------------------------------------------- */
int GP_OHD_Pantheon_output_dc(int num)
{
    int i, N=1048;
    FILE * fp=NULL;
    FILE *fp_dc = NULL;
    //FILE *fpcov = NULL;
    double z[N], sn, sig;
    //double *cov;

    // initial gp ohd with case: num
    initial_cosmo_gp(num);

    // open data file
    fp = fopen("/home/ekli/myworks/cosmodata/sn_full_dat.txt", "r");
    fp_dc = fopen("./data/Pantheon_out_dc.txt", "w");
    //fpcov = fopen("./data/Pantheon_out_cov.txt", "w");

    // read in data
    printf("\nRead in Pantheon data\n");
    for (i=0; i<N; i++) {
        printf("%5d-th\n", i);
        fscanf(fp, "%lf %lf %lf\n", &z[i], &sn, &sig);
        fprintf(fp_dc, "%12.8f %12.8f %12.8f\n",
                z[i], gp_dimles_comdistance(z[i]),
                sqrt( gp_cov_dimles_comdistance(z[i], z[i])) );
    }

    // calculation cov obs 
    //cov = (double *)calloc(N*N, sizeof(double));
    //for (i=0; i<N; i++) {
    //    printf("%5d-th\n", i);
    //    for (j=0; j<=i; j++) {
    //        cov[i*N+j] = gp_cov_dimles_comdistance(z[i], z[j]);
    //        cov[j*N+i] = cov[i*N+j];
    //    }
    //}

    //// write out cov
    //for (i=0; i<N*N; i++) {
    //    fprintf(fpcov, "%16.12f\n", cov[i]);
    //}

    // close data file
    fclose(fp);
    fclose(fp_dc);
    //fclose(fpcov);

    free_gapp();

    return _SUCCESS_;
}

int GP_OHD_StrongLens_output_dc(int num, int N, char * file_sl, char *file_slout)
{
    int i;
    FILE * fp=NULL;
    FILE *fp_dc = NULL;
    double zl, zs, tE, stE, sc, ssc, obs, sobs;
    double theta, sigtheta, Dobs, sigDobs;

    // initial gp ohd with case: num
    initial_cosmo_gp(num);

    // open data file
    fp = fopen(file_sl, "r");
    fp_dc = fopen(file_slout, "w");

    // read in data
    printf("\nRead in Pantheon data\n");
    for (i=0; i<N; i++) {
        printf("%5d-th\n", i);
        fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf\n", 
                &zl, &zs, &tE, &stE, &sc, &ssc, &obs, &sobs);
        theta = tE/3600*_PI_/180;
        sigtheta = stE/3600*_PI_/180;
        Dobs = pow(_C_/1000, 2) * theta/(4*_PI_*pow(sc, 2));
        sigDobs = Dobs * sqrt(pow(sigtheta/theta, 2) + 4*pow(ssc/sc, 2));
        fprintf(fp_dc, "%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n",
                zl, zs, 
                gp_dimles_comdistance(zl),
                gp_dimles_comdistance(zs),
                sqrt( gp_cov_dimles_comdistance(zl, zl)),
                sqrt( gp_cov_dimles_comdistance(zl, zs)),
                sqrt( gp_cov_dimles_comdistance(zs, zl)),
                sqrt( gp_cov_dimles_comdistance(zs, zs)),
                Dobs, sigDobs
                );
    }
    
    // close data file
    fclose(fp);
    fclose(fp_dc);
    //fclose(fpcov);

    free_gapp();

    return _SUCCESS_;
}

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/*         do the calculation in array          */
/* -------------------------------------------- */
// initial Pantheon data
int GP_OHD_initial_Pantheon(void)
{
    int i, j, N=1048;
    FILE * fp=NULL;
    FILE *fps=NULL;
    double z, s;

    // initial data array
    // initial obs of theory
    sn = (DATA *)malloc(sizeof(DATA));
    gp = (GP *)malloc(sizeof(GP));
    thgp = (GP *)malloc(sizeof(GP));
    sn->n = N;
    sn->x = (double *)calloc(N, sizeof(double));
    sn->y = (double *)calloc(N, sizeof(double));
    sn->s = (double *)calloc(N, sizeof(double));
    sn->cov = (double *)calloc(N*N, sizeof(double));
    gp->obs = (double *)calloc(N, sizeof(double));
    gp->cov = (double *)calloc(N*N, sizeof(double));
    thgp->obs = (double *)calloc(N, sizeof(double));
    thgp->cov = (double *)calloc(N*N, sizeof(double));

    // open data file
    fp = fopen("/home/ekli/myworks/cosmodata/sn_full_dat.txt", "r");
    fps = fopen("./data/Pantheon_out_dc.txt", "r");

    // read in data
    printf("\nRead in Pantheon data\n");
    for (i=0; i<N; i++) { 
        fscanf(fp, "%lf %lf %lf\n", &sn->x[i], &sn->y[i], &sn->s[i]);
        fscanf(fps, "%lf %lf %lf", &z, &gp->obs[i], &s);
        gp->cov[i*N+i] = s*s;
    }

    // close data file
    fclose(fp);
    fclose(fps);

    // open cov file
    fp = fopen("/home/ekli/myworks/cosmodata/sn_full_cov.txt", "r");

    // read in cov
    printf("\nRead in Pantheon sys cov\n");
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            fscanf(fp, "%lf\n", &sn->cov[i*N+j]);
            if(i==j)
                sn->cov[i*N+j] += pow(sn->s[i], 2.0);
            //printf("%8d %16.8f\n", i*N+j, sn->cov[i*N+j]);
        }
    }

    // close cov file
    fclose(fp);
    
    return _SUCCESS_;
}

// sn_loglikelihood
double GP_OHD_sn_loglikelihood(void)
{
    int i;
    double magA, magB, magC, chisq;
    double difmu[sn->n];
    //double logDet;

# pragma omp parallel for
    for (i=0; i<sn->n; i++) {
        difmu[i] = thgp->obs[i] - sn->y[i];
    }
    
    // get inv Cov
    gsl_matrix_view M = gsl_matrix_view_array(thgp->cov, sn->n, sn->n);
    gsl_linalg_cholesky_decomp(&M.matrix);
    gsl_linalg_cholesky_invert(&M.matrix);

    gsl_vector_view V = gsl_vector_view_array(difmu, sn->n);

    magC = 0.0;
//# pragma omp parallel for
    for (i=0; i<sn->n*sn->n; i++) {
        magC += thgp->cov[i];
    }

    gsl_vector *C = gsl_vector_calloc(sn->n);
    gsl_blas_dgemv(CblasNoTrans,
            1.0, &M.matrix, &V.vector,
            0.0, C);

    magB = 0.0;
//# pragma omp parallel for
    for (i=0; i<sn->n; i++) {
        magB += gsl_vector_get(C, i);
    }

    gsl_blas_ddot(&V.vector, C, &magA);

    gsl_vector_free(C);
    chisq = magA - pow(magB, 2.0)/magC + log(magC/_PI_);

    // check whether is nan
    if(gsl_isnan(chisq))
        return -3.0e30;
    return -0.5*chisq;
}

double GP_OHD_rofchi(double chi, double omk)
{
    if(omk > 1e-8) {
        return 1/sqrt(omk) * sinh(sqrt(omk) * chi);
    } else if (omk < -1e-8) {
        return 1/sqrt(-omk) * sin(sqrt(-omk) * chi);
    } else {
        return chi;
    }
}

double return_gp_ohd_loglike(double omk)
{
    int i, j, N;
    double dz[1048];
    
    // initial GP_OHD_initial_Pantheon
    if (call_again == 1) {
        printf("\nInitial Pantheon data\n");
        GP_OHD_initial_Pantheon();
        call_again = 0;
        printf("\nInitial Pantheon data Done!\n");
    }

    // initial mu^th
# pragma omp parallel for
    for (i=0; i<sn->n; i++) {
        dz[i] = GP_OHD_rofchi(gp->obs[i], omk);
        thgp->obs[i] = 5*log10((1+sn->x[i])*dz[i]);
    }

    N = sn->n;

# pragma omp parallel for
    for (i=0; i<N; i++) {
# pragma omp parallel for
        for (j=0; j<N; j++) {
            thgp->cov[i*N+j] = pow(5/log(10), 2.0) *
                sqrt(1+omk*pow(dz[i], 2)) / (dz[i]) *
                sqrt(1+omk*pow(dz[j], 2)) / (dz[j]) *
                gp->cov[i*N+j] +
                sn->cov[i*N+j];
        }
    }

    return GP_OHD_sn_loglikelihood();
}

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/*         return loglike of strong lens        */
/* -------------------------------------------- */
int GP_OHD_initial_StrongLens(void)
{
    return _SUCCESS_;
}
