#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model0_kouril7(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 1;
    rowvals[2] = 2;
    rowvals[3] = 7;
    rowvals[4] = 0;
    rowvals[5] = 1;
    rowvals[6] = 9;
    rowvals[7] = 10;
    rowvals[8] = 2;
    rowvals[9] = 5;
    rowvals[10] = 6;
    rowvals[11] = 8;
    rowvals[12] = 4;
    rowvals[13] = 5;
    rowvals[14] = 6;
}