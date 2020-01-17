/* ===============================================
 * File Name: Likelihood_HST.c
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-10-21 15:25:29
 * =============================================== 
 */

#include <stdio.h>
#include "common.h"

// HST data
static Obs_Data HST;

// ===============================
// initial HST data
// ===============================
int initial_HST(char * dat_file)
{
    FILE *fp;

    HST.info = "HST";
    HST.n = 1;
    HST.x = (double *)malloc(sizeof(double));
    HST.y = (double *)malloc(sizeof(double));
    HST.z = (double *)malloc(sizeof(double));

    if(NULL == (fp = fopen(dat_file, "r"))) {
        perror("fopen");
        return -1;
    }

    fscanf(fp, "%lf %lf %lf\n", &HST.x[0], &HST.y[0], &HST.z[0]);

    fclose(fp);
    return 0;
}

// ====================================
// HST likelihood
// ====================================
double Likelihood_HST(double *param,
                      void * model)
{
    double chisq = 0;
    chisq = pow((param[0] - HST.y[0])/HST.z[0], 2);
    return -0.5*chisq;
}
