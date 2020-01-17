/* ===============================================
 * File Name: cosmography.h
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-06-28 09:42:03
 * =============================================== 
 */

#ifndef __COSMOGRAPHY__
#define __COSMOGRAPHY__

int initial_cosmography(double omegak, double alpha, double beta);
double cosmography_dimless_distance(double z);
double cosmography_dimless_distance2(double zl, double zs);

double cosmography_Dobs_SL(double zl, double zs);


#endif /* __COSMOGRAPHY__ */

