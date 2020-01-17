/* ===============================================
 * File Name: kernel.c
 * Author: ekli
 * Mail: ekli_091@mail.dlut.edu.cn  
 * Created Time: 2019-06-21 15:24:57
 * ===============================================
 */

#include "common.h"
#include "gapp.h"

double kernel(double x, double y, double sigf, double lenf)
{
    return sigf*sigf * exp(-(x-y)*(x-y)/2.0/lenf/lenf);
}

