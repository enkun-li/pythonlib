/* ===============================================
 * File Name: cosmography.c
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-06-28 09:35:16
 * =============================================== 
 */

#include "common.h"
#include "utils.h"

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/*       cosmography model       */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

double a1, a2;
double omk;

// initial cosmography
int initial_cosmography(double omegak, double alpha, double beta)
{
    omk = omegak;
    a2 = alpha;
    a2 = beta;
    return _SUCCESS_;
}

// dimensionless comoving angular distance
double cosmography_dimless_distance(double z)
{
    return z + a1 * pow(z, 2.0) + a2 * pow(z, 3.0);
}

// d(zl, zs)
double cosmography_dimless_distance2(double zl, double zs)
{
    double dl, ds;
    dl = cosmography_dimless_distance(zl);
    ds = cosmography_dimless_distance(zs);
    return ds * sqrt(1+omk*pow(dl,2)) - dl * sqrt(1+omk*pow(ds,2));
}

double cosmography_Dobs_SL(double zl, double zs)
{
    return cosmography_dimless_distance2(zl, zs)/
        cosmography_dimless_distance(zs);
}
