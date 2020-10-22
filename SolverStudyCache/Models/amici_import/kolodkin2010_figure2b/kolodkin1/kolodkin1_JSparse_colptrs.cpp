#include "sundials/sundials_types.h"

void JSparse_colptrs_kolodkin1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 3;
    colptrs[2] = 6;
    colptrs[3] = 9;
}