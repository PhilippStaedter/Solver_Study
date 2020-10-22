#include "sundials/sundials_types.h"

void dwdx_colptrs_model2_sarma7(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 4;
    colptrs[3] = 6;
    colptrs[4] = 7;
    colptrs[5] = 8;
    colptrs[6] = 10;
    colptrs[7] = 14;
    colptrs[8] = 18;
    colptrs[9] = 20;
    colptrs[10] = 24;
    colptrs[11] = 30;
    colptrs[12] = 31;
    colptrs[13] = 33;
}