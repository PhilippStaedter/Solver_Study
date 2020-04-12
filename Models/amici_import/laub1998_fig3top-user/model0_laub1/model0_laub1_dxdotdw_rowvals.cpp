#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model0_laub1(sunindextype *rowvals){
    rowvals[0] = 4;
    rowvals[1] = 4;
    rowvals[2] = 0;
    rowvals[3] = 0;
    rowvals[4] = 5;
    rowvals[5] = 5;
    rowvals[6] = 2;
    rowvals[7] = 2;
    rowvals[8] = 6;
    rowvals[9] = 6;
    rowvals[10] = 3;
    rowvals[11] = 3;
    rowvals[12] = 1;
    rowvals[13] = 1;
}