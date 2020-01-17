/* ===============================================
 * File Name: utils.c
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-10-21 16:42:17
 * =============================================== 
 */

#include <stdio.h>
#include <string.h>

/* ===============================
 * get ini key string:
 *
 * ===============================
 * */
int ini_file_read_string(char* key, char* filename)
{
    FILE * fp;
    int flag = 0;
    char sLine[1024];
    char *wTmp;

    if(NULL == (fp = fopen(filename, "r"))) {
        perror("fopen");
        return -1;
    }

    while (NULL != fgets(sLine, 1024, fp)) {
        if(0 == strncmp("//", sLine, 2)) continue;
        if('#' == sLine[0]) {
            wTmp = strchr(sLine, "=");
        }
        if ((NULL != wTmp) && (1 == flag)) {
            if(0 == strcmp(key, sLine, strlen(key))) {
                sLine[strlen(sLine) -1] = '\0';
                fclose(fp);
                while (*(wTmp +1) == ' ') {
                    wTmp++;
                }
                strcpy(buf, wTmp+1);
                return 0;
            }
        } else {
            if (0 == strncmp(sT))
        }
    }
}

