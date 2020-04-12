#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model0_bachar1(sunindextype *rowvals){
    rowvals[0] = 3;
    rowvals[1] = 1;
    rowvals[2] = 3;
    rowvals[3] = 3;
    rowvals[4] = 1;
    rowvals[5] = 1;
    rowvals[6] = 2;
    rowvals[7] = 1;
    rowvals[8] = 2;
    rowvals[9] = 2;
    rowvals[10] = 2;
}