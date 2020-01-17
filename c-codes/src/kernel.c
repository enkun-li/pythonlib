/* ===============================================
 * File Name: kernel.c
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-06-10 11:50:14
 * =============================================== 
 */

#include "common.h"

double kernel(double x, double y, double sigf, double lenf)
{
    return pow(sigf, 2.0) * exp(- pow(x-y, 2.0)/2.0/pow(lenf, 2.0));
}

