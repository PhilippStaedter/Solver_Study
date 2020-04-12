#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_essunger1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 7;
    colptrs[2] = 13;
    colptrs[3] = 15;
    colptrs[4] = 17;
    colptrs[5] = 21;
    colptrs[6] = 23;
    colptrs[7] = 26;
}