#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model0_karin5(sunindextype *rowvals){
    rowvals[0] = 1;
    rowvals[1] = 2;
    rowvals[2] = 0;
    rowvals[3] = 3;
    rowvals[4] = 4;
}