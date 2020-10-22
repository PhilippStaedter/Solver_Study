#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_sarma3(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 4;
    colptrs[2] = 8;
    colptrs[3] = 9;
    colptrs[4] = 14;
    colptrs[5] = 21;
    colptrs[6] = 28;
    colptrs[7] = 34;
    colptrs[8] = 40;
    colptrs[9] = 46;
    colptrs[10] = 50;
    colptrs[11] = 54;
}