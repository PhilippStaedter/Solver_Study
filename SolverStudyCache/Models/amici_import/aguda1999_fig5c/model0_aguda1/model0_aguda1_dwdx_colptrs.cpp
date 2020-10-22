#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_aguda1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 4;
    colptrs[2] = 6;
    colptrs[3] = 7;
    colptrs[4] = 11;
    colptrs[5] = 17;
    colptrs[6] = 20;
    colptrs[7] = 21;
    colptrs[8] = 23;
    colptrs[9] = 27;
    colptrs[10] = 30;
    colptrs[11] = 31;
    colptrs[12] = 31;
}