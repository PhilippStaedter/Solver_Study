#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model0_ma(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 6;
    rowvals[2] = 5;
    rowvals[3] = 5;
    rowvals[4] = 1;
    rowvals[5] = 1;
    rowvals[6] = 0;
    rowvals[7] = 3;
    rowvals[8] = 3;
    rowvals[9] = 2;
    rowvals[10] = 2;
    rowvals[11] = 4;
    rowvals[12] = 4;
    rowvals[13] = 6;
}