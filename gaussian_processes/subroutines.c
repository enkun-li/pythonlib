/* ===============================================
 * File Name: subroutines.c
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-05-19 10:18:37
 * =============================================== 
 */

#include <stdio.h>
#include <math.h>

typedef struct params{
    double sigf;
    double lenf;
} Pars;

double kernel(double x, double y, double sigf, double lenf) 
{
    return pow(sigf, 2) * exp(-pow(x-y, 2)/2/pow(lenf, 2));
}

double rec_mu_x(double x, )
