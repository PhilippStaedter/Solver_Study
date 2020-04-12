#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model0_goldbeter2(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 0;
    rowvals[2] = 0;
    rowvals[3] = 1;
    rowvals[4] = 1;
    rowvals[5] = 2;
    rowvals[6] = 2;
}