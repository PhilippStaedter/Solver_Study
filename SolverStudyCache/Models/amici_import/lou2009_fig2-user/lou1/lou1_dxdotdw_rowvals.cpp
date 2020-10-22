#include "sundials/sundials_types.h"

void dxdotdw_rowvals_lou1(sunindextype *rowvals){
    rowvals[0] = 4;
    rowvals[1] = 0;
    rowvals[2] = 3;
    rowvals[3] = 3;
    rowvals[4] = 0;
    rowvals[5] = 1;
    rowvals[6] = 4;
    rowvals[7] = 4;
    rowvals[8] = 1;
    rowvals[9] = 5;
    rowvals[10] = 2;
    rowvals[11] = 5;
    rowvals[12] = 5;
    rowvals[13] = 2;
    rowvals[14] = 3;
}