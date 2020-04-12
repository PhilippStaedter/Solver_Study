#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model1_alexander1(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 0;
    rowvals[2] = 3;
    rowvals[3] = 0;
    rowvals[4] = 3;
    rowvals[5] = 4;
    rowvals[6] = 4;
    rowvals[7] = 2;
    rowvals[8] = 0;
    rowvals[9] = 4;
    rowvals[10] = 2;
    rowvals[11] = 3;
}