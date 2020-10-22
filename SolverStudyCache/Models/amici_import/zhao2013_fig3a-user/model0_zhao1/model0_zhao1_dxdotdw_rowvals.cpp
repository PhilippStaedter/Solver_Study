#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model0_zhao1(sunindextype *rowvals){
    rowvals[0] = 4;
    rowvals[1] = 1;
    rowvals[2] = 2;
    rowvals[3] = 4;
    rowvals[4] = 4;
    rowvals[5] = 3;
    rowvals[6] = 1;
    rowvals[7] = 3;
    rowvals[8] = 3;
    rowvals[9] = 2;
    rowvals[10] = 2;
    rowvals[11] = 1;
}