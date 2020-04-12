#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model0_miao2008(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 0;
    rowvals[2] = 1;
    rowvals[3] = 0;
    rowvals[4] = 3;
    rowvals[5] = 0;
    rowvals[6] = 1;
    rowvals[7] = 2;
    rowvals[8] = 3;
    rowvals[9] = 1;
    rowvals[10] = 1;
    rowvals[11] = 2;
    rowvals[12] = 3;
    rowvals[13] = 2;
    rowvals[14] = 3;
    rowvals[15] = 2;
}