/* ===============================================
 * File Name: test.c
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-05-16 16:21:13
 * =============================================== 
 */

#include <stdio.h>
#include "common.h"

struct DATA {
    double * x;
    double * y;
    double * sig;
};

int main()
{
    struct DATA data;
    int n, i;

    printf("Input size of data:\n");
    scanf("%d", &n);

    data.x = class_malloc(double, n);
    data.y = class_malloc(double, n);
    data.sig = class_malloc(double, n);

    for (i=0;i<n;i++){
        printf("%d : %6.3f\n", i, data.x[i]);
    }
    return 0;
}
