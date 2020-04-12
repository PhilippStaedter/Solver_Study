#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model0_bier2(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 2;
    rowvals[2] = 1;
    rowvals[3] = 3;
    rowvals[4] = 0;
    rowvals[5] = 1;
    rowvals[6] = 2;
    rowvals[7] = 3;
    rowvals[8] = 2;
    rowvals[9] = 3;
}