#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model0_marhl(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 3;
    rowvals[2] = 2;
    rowvals[3] = 3;
    rowvals[4] = 4;
    rowvals[5] = 2;
    rowvals[6] = 3;
    rowvals[7] = 4;
    rowvals[8] = 0;
    rowvals[9] = 3;
    rowvals[10] = 0;
    rowvals[11] = 3;
    rowvals[12] = 1;
    rowvals[13] = 3;
    rowvals[14] = 1;
    rowvals[15] = 3;
}