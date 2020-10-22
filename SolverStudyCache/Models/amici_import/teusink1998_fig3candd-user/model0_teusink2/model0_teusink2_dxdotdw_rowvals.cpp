#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model0_teusink2(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 3;
    rowvals[2] = 0;
    rowvals[3] = 1;
    rowvals[4] = 3;
    rowvals[5] = 0;
    rowvals[6] = 1;
    rowvals[7] = 0;
}