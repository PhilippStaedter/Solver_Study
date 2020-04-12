#include "sundials/sundials_types.h"

void dwdx_colptrs_model2_essunger3(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 5;
    colptrs[2] = 9;
    colptrs[3] = 11;
    colptrs[4] = 13;
    colptrs[5] = 17;
    colptrs[6] = 19;
    colptrs[7] = 23;
}