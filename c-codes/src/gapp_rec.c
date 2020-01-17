/* ===============================================
 * File Name: gapp_rec.c
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-06-10 13:19:49
 * =============================================== 
 */

#include "common.h"
#include "utils.h"
#include "gapp.h"

/* @@@@@@@@@@@@@@@@@@@ */
/*   base equations    */
/* @@@@@@@@@@@@@@@@@@@ */

// func: 1/mu
double func_one_over_mu(
        double x, // x position
        double *args // pars of func
        )
{
    return pow(1+x, *(args))/pow(rec_mu(x), *(args+1));
}


/* covariance of integrate_one_over_mu */

double func_covxy_over_muxy(double x, double y)
{
    return rec_covariance(x,y)/rec_mu(x)/rec_mu(y);
}

double func_fxy_covxy_over_muxy(double x, double y, double *args)
{
    return pow(1+x, *(args))
        *pow(1+y, *(args+1))
        /pow(rec_mu(x), *(args+2))
        /pow(rec_mu(y), *(args+3))
        * func_covxy_over_muxy(x,y);
}

double func_fx_covxj_over_muxj(double x, double *args)
{
    double y = *(args);
    return pow(1+x, *(args+1))
        /pow(rec_mu(x), *(args+2))
        * func_covxy_over_muxy(x,y);
}

// func: covxy/mu_x mu_y
double func_covxy_over_mux_muy(
        double x, // x
        double y, // y
        double *args // pars
        )
{
    return pow(1+x, *(args))
       * pow(1+y, *(args))
       * rec_covariance(x, y)
        / pow(rec_mu(x), *(args+1)) 
        /pow(rec_mu(y), *(args+1));
}

// func: (1+x)^arg[1] * cov[x, arg[0]]/mu_1^arg[2] mu(arg[0])^arg[3]
double func_covxj_over_mux_muj(
        double x, // x
        double *args // pars
        )
{
    double y = *(args);

    return pow(1+x, *(args+1))
       * rec_covariance(x, y)
        / pow(rec_mu(x), *(args+2)) 
        / pow(rec_mu(y), *(args+3));
}


/* @@@@@@@@@@@@@@@@@@@@@@@@ */
/* some integrate functions */
/* @@@@@@@@@@@@@@@@@@@@@@@@ */

// integrate (1+x)^al / mu^be
double integrate_fx_over_fmu(
        double a, double b, // lower and upper limit
        double *args        // some pars of func
        )
{
    return qgausleg(&func_one_over_mu, a, b, &args[0]);
}


// integrate (1+x)^al / mu^be *cov(x, j)
double integrate_fxcovxj_over_fmux(
        double a, double b, // lower and upper limit
        double *args        // some pars of func
        )
{
    return qgausleg(&func_covxj_over_mux_muj, a, b, &args[0]);
}

// integrate (1+x)^al (1+y)^al / mux^be / muy^be * cov(x,y)
double integrate_fxcovxy_over_fmuxy(
        double a, // lower limit of i
        double b, // upper limit of i
        double c, // lower limit of j
        double d, // upper limit of j
        double *args // some pars of func
        )
{
    return qgausleg2d_fix(&func_covxy_over_mux_muy, a, b, c, d, &args[0]);
}

/* @@@@@@@@@@@@@@@@@@@ */
/* some test functions */
/* @@@@@@@@@@@@@@@@@@@ */

// integrate func: 1/mu
double integrate_one_over_mu(
        double a, // lower limit
        double b  // upper limit
        )
{
    double pars[] = {0.0, 1.0};
    //return qgausleg(&func_one_over_mu, a, b, &pars[0]);
    return integrate_fx_over_fmu(a,b,&pars[0]);
}

// cov of integrate func: 1/mu
double integrate_cov_one_over_mu(
        double a, // lower limit of i
        double b, // upper limit of i
        double c, // lower limit of j
        double d  // upper limit of j
        )
{
    double pars[] = {0.0, 2.0};
    // return qgausleg2d_fix(&func_covxy_over_mux_muy, a, b, c, d, &pars[0]);
    return integrate_fxcovxy_over_fmuxy(a, b, c, d, &pars[0]);
}

/* @@@@@@@@@@@@@@@@@@@@@@@@@ */
/* some special function rec */
/* @@@@@@@@@@@@@@@@@@@@@@@@@ */

// func: mu0/mu
double func_mu0_over_mu(
        double x, // x position
        double *args // pars of func
        )
{
    return pow(1+x, *(args)) * pow(rec_mu(0.0)/rec_mu(x), *(args+1));
}

// rec mu/mu0
double rec_mu_over_mu0(double x)
{
    double mu0;
    mu0 = rec_mu(0.0);
    return rec_mu(x)/mu0;
}

double rec_cov_mu_over_mu0(double x, double y)
{
    double mux, muy, mu0;
    double cov00, covx0, cov0y, covxy;
    double pars[] = {0.0, 1.0};

    mu0 = rec_mu(0.0);
    mux = rec_mu(x);
    muy = rec_mu(y);

    cov00 = func_covxy_over_mux_muy(0,0,&pars[0]);
    covx0 = func_covxy_over_mux_muy(x,0,&pars[0]);
    cov0y = func_covxy_over_mux_muy(0,y,&pars[0]);
    covxy = func_covxy_over_mux_muy(x,y,&pars[0]);

    return mux/mu0 * muy/mu0 * (covxy - covx0 - cov0y + cov00);
}

/* @@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* some special base function */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@ */

// func: covxy/mu/mu0_x mu/mu0_y
double func_mu0_mu0_covxy_over_mux_muy(
        double x, // x
        double y, // y
        double *args // pars
        )
{
    return pow(1+x, *(args))
       * pow(1+y, *(args))
       * rec_cov_mu_over_mu0(x, y)
        / pow(rec_mu_over_mu0(x), *(args+1)) 
        /pow(rec_mu_over_mu0(y), *(args+1));
}

// func: (1+x)^arg[1] * cov[x, arg[0]]/mu_1^arg[2] mu(arg[0])^arg[3]
double func_mu0_mu0_covxj_over_mux_muj(
        double x, // x
        double *args // pars
        )
{
    double y = *(args);

    return pow(1+x, *(args+1))
       * rec_cov_mu_over_mu0(x, y)
        / pow(rec_mu_over_mu0(x), *(args+2)) 
        / pow(rec_mu_over_mu0(y), *(args+3));
}

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* some special integrate function */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

// integrate (1+x)^al / mu^be *mu0^be
double integrate_fxmu0_over_fmu(
        double a, double b, // lower and upper limit
        double *args        // some pars of func
        )
{
    return qgausleg(&func_mu0_over_mu, a, b, &args[0]);
}

// integrate (1+x)^al / mu^be *cov(x, j)
double integrate_fxmu0covxj_over_fmux(
        double a, double b, // lower and upper limit
        double *args        // some pars of func
        )
{
    return qgausleg(&func_mu0_mu0_covxj_over_mux_muj, a, b, &args[0]);
}

// integrate (1+x)^al (1+y)^al / mux^be / muy^be * cov(x,y)
double integrate_fxmu0covxy_over_fmuxy(
        double a, // lower limit of i
        double b, // upper limit of i
        double c, // lower limit of j
        double d, // upper limit of j
        double *args // some pars of func
        )
{
    return qgausleg2d_fix(&func_mu0_mu0_covxy_over_mux_muy, a, b, c, d, &args[0]);
}


/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* still some integrate function */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

// rec mu * int 1/mu dx
double rec_mu_int_one_over_mu(double x)
{
    return rec_mu(x) * integrate_one_over_mu(0.0, x);
}

double rec_cov_mu_int_one_over_mu(double x, double y)
{
    double fx, fy, mux, muy;
    double cov00, covx0, cov0y, covxy;
    double pars00[] = {0.0, 1.0};
    double parsx0[] = {y, 0.0, 1.0+1.0, 1.0};
    double pars0y[] = {x, 0.0, 1.0+1.0, 1.0};
    double parsxy[] = {0.0, 2.0};

    mux = rec_mu(x);
    muy = rec_mu(y);

    fx = integrate_one_over_mu(0.0,x);
    fy = integrate_one_over_mu(0.0,y);

    cov00 = func_covxy_over_mux_muy(x,y,&pars00[0]);
    covx0 = integrate_fxcovxj_over_fmux(0.0, x, &parsx0[0]);
    cov0y = integrate_fxcovxj_over_fmux(0.0, y, &pars0y[0]);
    covxy = integrate_fxcovxy_over_fmuxy(0.0, x, 0.0, y, &parsxy[0]);

    return fx*fy*cov00 - mux * fy *covx0 - fx *muy * cov0y + mux *muy *covxy;
}
