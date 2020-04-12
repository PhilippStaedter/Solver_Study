#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model0_borghans2(sunindextype *rowvals){
    rowvals[0] = 2;
    rowvals[1] = 3;
    rowvals[2] = 2;
    rowvals[3] = 3;
    rowvals[4] = 2;
    rowvals[5] = 3;
    rowvals[6] = 1;
    rowvals[7] = 3;
    rowvals[8] = 0;
    rowvals[9] = 0;
    rowvals[10] = 0;
    rowvals[11] = 1;
    rowvals[12] = 3;
}