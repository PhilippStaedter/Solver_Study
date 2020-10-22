#include "sundials/sundials_types.h"

void dwdx_colptrs_model0_teusink2(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 4;
    colptrs[2] = 5;
    colptrs[3] = 6;
    colptrs[4] = 8;
}