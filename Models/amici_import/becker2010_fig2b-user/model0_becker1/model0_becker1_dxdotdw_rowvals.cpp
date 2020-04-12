#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model0_becker1(sunindextype *rowvals){
    rowvals[0] = 1;
    rowvals[1] = 1;
    rowvals[2] = 0;
    rowvals[3] = 1;
    rowvals[4] = 2;
    rowvals[5] = 0;
    rowvals[6] = 1;
    rowvals[7] = 2;
    rowvals[8] = 2;
    rowvals[9] = 3;
    rowvals[10] = 0;
    rowvals[11] = 1;
    rowvals[12] = 3;
    rowvals[13] = 3;
    rowvals[14] = 5;
    rowvals[15] = 3;
    rowvals[16] = 4;
}