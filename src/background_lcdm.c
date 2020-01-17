/* ===============================================
 * File Name: CosmoTheory.c
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-10-21 18:11:24
 * =============================================== 
 */

#include "common.h"
#include "models.h"
#include "utils.h"

double lcdm_Hofz(struct CosmoTheory this, double z)
{
    double H0 = this.param[0];
    double omm = this.param[1];
    return H0*sqrt(omm*pow(1+z, 3) + 1-omm);
}

double lcdm_one_over_Hz(struct CosmoTheory this, double z)
{
    return 1.0/lcdm_Hofz(this, z);
}

double lcdm_D_C(double c)
{
    double pars[] = {0.0}
}

CosmoTheory *LCDM_model(CosmoTheory *this, double *pars)
{
    this->model = "LCDM";
    this->param = pars;
    this->initial_model = initial_model;
    this->Hofz = lcdm_Hofz;
    this->D_A = lcdm_D_A;

    return 0;
}

double test_func(double z, void * model)
{
    CosmoTheory *cosmo = (CosmoTheory*)model;
    return cosmo->Hofz(*cosmo, z);
}

int main()
{
    int i;
    double *pars = (double *)malloc(2*sizeof(double));
    pars[0] = 67.27;
    pars[1] = 0.3;
    double x;

    CosmoTheory lcdm;

    LCDM_model(&lcdm, pars);

    printf("%s\n", lcdm.model);
    printf("%d\n", lcdm.initial_model(lcdm));

    for (i=0; i<10; i++) {
        x = i*0.3;
        printf("%d : %f %f\n", i, lcdm.Hofz(lcdm, (double)i*0.3), 
                test_func(x, &lcdm));
    }
    return 0;
}
