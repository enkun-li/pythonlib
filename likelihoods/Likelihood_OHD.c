/* ===============================================
 * File Name: Likelihood_OHD.c
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-10-21 18:03:02
 * =============================================== 
 */

#include <stdio.h>
#include "common.h"

// OHD data
static Obs_Data OHD;

// ===============================
// initial OHD data
// ===============================
int initial_OHD(int N, char * dat_file)
{
    int i;
    FILE *fp;

    OHD.info = "OHD";
    OHD.n = N;
    OHD.x = (double *)malloc(N*sizeof(double));
    OHD.y = (double *)malloc(N*sizeof(double));
    OHD.z = (double *)malloc(N*sizeof(double));

    if(NULL == (fp = fopen(dat_file, "r"))) {
        perror("fopen");
        return -1;
    }

    for (i=0; i<N; i++) {
        fscanf(fp, "%lf %lf %lf\n", &OHD.x[i], &OHD.y[i], &OHD.z[i]);
    }

    fclose(fp);
    return 0;
}

// ====================================
// OHD likelihood
// ====================================
double Likelihood_OHD(double *param,
                      void * model)
{
    int i;
    double zi, Hzi;
    double chisq = 0;

    for (i=0; i<OHD.n; i++) {
        zi = OHD.x[i];
        Hzi = model->Hofz(zi);
        chisq += pow((param[i] - OHD.y[i])/OHD.z[i], 2);
    }
    return -0.5*chisq;
}
