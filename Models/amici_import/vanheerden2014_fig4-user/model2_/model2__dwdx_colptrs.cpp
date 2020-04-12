#include "sundials/sundials_types.h"

void dwdx_colptrs_model2_(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 0;
    colptrs[2] = 0;
    colptrs[3] = 0;
    colptrs[4] = 3;
    colptrs[5] = 5;
    colptrs[6] = 5;
    colptrs[7] = 7;
}