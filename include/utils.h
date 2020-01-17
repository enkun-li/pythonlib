/* ===============================================
 * File Name: utils.h
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-10-21 20:51:37
 * =============================================== 
 */

#ifndef __UTILS__
#define __UTILS__

#include <stdio.h>

/* gaus leg integrate */
int gauleg(double x1, double x2, double x[], double w[]);
double qgausleg(double (*func)(double, double *), double a, double b, double *args );
double qgausleg2d_fix(double (*func)(double, double, double *), double a, double b, double c, double d, double *args );

/* */
double qgaus(double (*func)(double, double *), double a, double b, double *args);
double qgaus2d_fix(double (*func)(double, double, double *), double a, double b, double c, double d, double *args);

/* inverse matrix */
int gsl_matrix_inv(gsl_matrix *M);

/* gsl integrate */
double quad2d(double (*func)(double, double),
        double x1,
        double x2,
        double y1, 
        double y2);

double quad1d(double (*func)(double),
        double x,
        double y);

#endif /* __UTILS__ */
