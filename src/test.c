/* ===============================================
 * File Name: test.c
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-10-21 17:25:35
 * =============================================== 
 */

#include <stdio.h>
#include "common.h"
#include "likelihoods.h"

void model()
{
}

int main()
{
    char * files = "HST2019.txt";
    double param[] = {70.0};

    initial_HST(files);
    printf("%f\n", Likelihood_HST(param, model));

    return 0;
}
