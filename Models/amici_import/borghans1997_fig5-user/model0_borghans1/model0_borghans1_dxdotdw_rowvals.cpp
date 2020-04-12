#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model0_borghans1(sunindextype *rowvals){
    rowvals[0] = 3;
    rowvals[1] = 4;
    rowvals[2] = 3;
    rowvals[3] = 4;
    rowvals[4] = 3;
    rowvals[5] = 4;
    rowvals[6] = 0;
    rowvals[7] = 4;
    rowvals[8] = 1;
    rowvals[9] = 2;
    rowvals[10] = 1;
    rowvals[11] = 2;
    rowvals[12] = 0;
    rowvals[13] = 4;
}