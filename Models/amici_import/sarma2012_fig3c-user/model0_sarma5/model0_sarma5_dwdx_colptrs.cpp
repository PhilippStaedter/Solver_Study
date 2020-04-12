#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_sarma5(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 4;
    colptrs[3] = 6;
    colptrs[4] = 8;
    colptrs[5] = 12;
    colptrs[6] = 16;
    colptrs[7] = 18;
    colptrs[8] = 22;
    colptrs[9] = 27;
    colptrs[10] = 28;
    colptrs[11] = 30;
}