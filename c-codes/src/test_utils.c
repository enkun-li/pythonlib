/* ===============================================
 * File Name: test_utils.c
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-06-10 10:38:38
 * =============================================== 
 */

#include "common.h"
#include "utils.h"

/* ============================== */
double func(double x, double * args)
{
    return *(args) * sin(*(args+1)*x);
}

int test_qgaus(void)
{
    int i;
    double x, ff, ff2, rr;
    double pars[] = {1.0, 1.0};

    for(i=0; i<51; i++){
        x = _PI_*2.0/50*i;
        ff = qgaus(&func, 0.0, x, &pars[0]);
        ff2 = qgausleg(&func, 0.0, x, &pars[0]); 
        rr = 1-cos(x);
        printf("int sinx dx = %12.9f  real = %12.9f  Delta = %12.9e  ", ff, rr, ff-rr );
        printf("leg = %12.9f  Delta = %12.9e\n", ff2, ff2-rr);
    }
    return _SUCCESS_;
}

/* ============================== */
double func_2d(double x, double y, double * args)
{
    return pow(x, *(args)) + pow(y, *(args+1)) +x*y;
}

int test_qgaus2d_fix(void)
{
    int i;
    double x, ff, ff2, rr;
    double pars[] = {2.0, 2.0};

    for(i=0; i<51; i++){
        x = 2*_PI_/50*i;
        ff = qgaus2d_fix(&func_2d, 0.0, x, 0.0, x, &pars[0]);
        rr = 11.0/12*pow(x, 4.0);
        // rr = pow(x, 3); // 2.0/3.0*pow(x,4);
        printf("int sinx dx = %15.9f  real = %15.9f  Delta = %15.9e  ", ff, rr, ff-rr );
        ff2 = qgausleg2d_fix(&func_2d, 0.0, x, 0.0, x, &pars[0]);
        printf("leg = %15.9f  Delta = %15.9e\n", ff2, ff2-rr);
    }
    return _SUCCESS_;
}

/* ============================== */
int main(void)
{
    int status;

    status = test_qgaus();
    status = test_qgaus2d_fix();

    if(status == _FAILURE_)
        printf("Failed");

    return _SUCCESS_;
}
