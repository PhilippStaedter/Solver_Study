#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model0_perelson1(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 0;
    rowvals[2] = 1;
    rowvals[3] = 2;
    rowvals[4] = 1;
    rowvals[5] = 3;
    rowvals[6] = 1;
    rowvals[7] = 3;
}