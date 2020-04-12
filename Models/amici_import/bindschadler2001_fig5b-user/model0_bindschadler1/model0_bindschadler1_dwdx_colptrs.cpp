#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_bindschadler1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 5;
    colptrs[2] = 10;
    colptrs[3] = 13;
    colptrs[4] = 16;
}