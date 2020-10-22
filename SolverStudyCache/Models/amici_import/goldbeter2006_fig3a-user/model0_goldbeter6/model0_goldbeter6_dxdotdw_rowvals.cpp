#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model0_goldbeter6(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 0;
    rowvals[2] = 1;
    rowvals[3] = 1;
    rowvals[4] = 2;
    rowvals[5] = 2;
}