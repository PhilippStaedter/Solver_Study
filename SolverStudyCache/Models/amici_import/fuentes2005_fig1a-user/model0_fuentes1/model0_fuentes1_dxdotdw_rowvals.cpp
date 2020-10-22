#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model0_fuentes1(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 2;
    rowvals[2] = 3;
    rowvals[3] = 0;
    rowvals[4] = 1;
    rowvals[5] = 3;
    rowvals[6] = 0;
    rowvals[7] = 1;
    rowvals[8] = 2;
}