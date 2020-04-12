#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model3_naresh2(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 3;
    colptrs[3] = 4;
    colptrs[4] = 5;
    colptrs[5] = 6;
    colptrs[6] = 8;
    colptrs[7] = 10;
    colptrs[8] = 11;
    colptrs[9] = 13;
    colptrs[10] = 15;
    colptrs[11] = 17;
    colptrs[12] = 18;
    colptrs[13] = 20;
}