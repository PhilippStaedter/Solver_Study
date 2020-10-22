#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_essunger6(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 5;
    colptrs[2] = 9;
    colptrs[3] = 11;
    colptrs[4] = 13;
    colptrs[5] = 15;
    colptrs[6] = 17;
    colptrs[7] = 21;
}