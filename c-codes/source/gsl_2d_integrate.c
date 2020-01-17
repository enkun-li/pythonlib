/* ===============================================
 * File Name: gsl_2d_integrate.c
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-06-30 23:33:21
 * =============================================== 
 */

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

#define EPS 1e-7

// define a func1d to integrate
static double (*func1d)(double);

// define a func2d to integrate
static double (*func2d)(double, double);

// general gsl qags integrate
int gsl_qags(double (*func)(double, void *), 
        double x, double y, 
        double *result, double *error,
        void * args)
{
    gsl_integration_workspace * w
        = gsl_integration_workspace_alloc (1000);

    gsl_function F;
    F.function = (*func);
    F.params = args;

    gsl_integration_qags(&F, x, y, 0, EPS, 1000,
            w, result, error);

    gsl_integration_workspace_free(w);
    
    return 0;
}

// \int f(x) dx
double func_1d(double x, void *pars)
{
    return func1d(x);
}

double quad1d(double (*func)(double),
        double x,
        double y)
{
    double res, err;
    double pars = 0;
    func1d = func;
    gsl_qags(&func_1d, x, y, &res, &err, &pars);
    return res;
}

// \int \int h(x,y) dy dx
double func_f(double y, void *pars)
{
    double x = *(double*)pars;
    return func2d(x,y);
}

double func_h(double x, void *pars)
{
    double res, err;
    double *y = (double *)pars;

    gsl_qags(&func_f, *y, *(y+1), &res, &err, &x);
    return res;
}

double quad2d(double (*func)(double, double),
        double x1,
        double x2,
        double y1, 
        double y2)
{
    double res, err;
    double pars[] = {y1, y2};

    func2d = func;
    gsl_qags(&func_h, x1, x2, &res, &err, pars);

    return res;
}

//double f(double x, double y)
//{
//    return x*y+x+ y* sin(y) + 1/sqrt(1+x+y)*exp(-x*y);
//}
//
//int main(void)
//{
//    double res = quad2d(&f, 0,1,0,1);
//    printf("result is % .8f\n", res);
//    return 0;
//}
