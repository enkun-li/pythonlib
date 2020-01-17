/* ===============================================
 * File Name: cosmo_gp.h
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-07-01 16:31:02
 * =============================================== 
 */

#ifndef __COSMO_GP__
#define __COSMO_GP__

int initial_cosmo_gp(int num);
int setup_cosmo_gp(double omegak, int ohdtype);

double gp_hz(double z);
double gp_cov_hz(double zi, double zj);

double gp_dimles_comdistance(double z);
double gp_cov_dimles_comdistance(double zi, double zj);

double gp_dimles_comangdistance(double z);
double gp_cov_dimles_comangdistance(double zi, double zj);

double gp_com_distance(double z);
double gp_ang_distance(double z);
double gp_lum_distance(double z);

double gp_DistanceModulus_star(double z);
double gp_cov_DistanceModulus_star(double zi, double zj);

double gp_DistanceModulus(double z);

double gp_DistanceSumRole(double zl, double zs);
double gp_cov_DistanceSumRole(double zli, double zsi, double zlj, double zsj);

int GP_OHD_Pantheon_output_dc(int num);
int GP_OHD_StrongLens_output_dc(int num, int N, char * file_sl, char *file_slout);

int GP_OHD_initial_Pantheon(void);
double GP_OHD_sn_loglikelihood(void);
double GP_OHD_rofchi(double chi, double omk);
double return_gp_ohd_loglike(double omk);

#endif /* __COSMO_GP__ */
