/* ===============================================
 * File Name: models.h
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-10-21 20:23:46
 * =============================================== 
 */

#ifndef __MODELS__
#define __MODELS__

#include <stdio.h>

// define a Cosmo Theory class
typedef struct CosmoTheory{
    char *model;
    double *param;
    double (*Hofz)(struct CosmoTheory This, double z);
    double (*D_A)(struct CosmoTheory This, double z);
    double (*D_L)(struct CosmoTheory This, double z);
} CosmoTheory;

#endif /* __MODELS__ */
