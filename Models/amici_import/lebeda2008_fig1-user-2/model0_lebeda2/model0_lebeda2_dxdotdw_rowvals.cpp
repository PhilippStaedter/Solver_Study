#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model0_lebeda2(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 1;
    rowvals[2] = 0;
    rowvals[3] = 3;
    rowvals[4] = 2;
    rowvals[5] = 3;
}