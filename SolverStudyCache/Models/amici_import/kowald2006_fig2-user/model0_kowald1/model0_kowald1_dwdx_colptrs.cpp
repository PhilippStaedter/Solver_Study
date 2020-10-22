#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_kowald1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 6;
    colptrs[2] = 11;
    colptrs[3] = 14;
    colptrs[4] = 17;
    colptrs[5] = 19;
    colptrs[6] = 20;
    colptrs[7] = 21;
    colptrs[8] = 23;
    colptrs[9] = 24;
}