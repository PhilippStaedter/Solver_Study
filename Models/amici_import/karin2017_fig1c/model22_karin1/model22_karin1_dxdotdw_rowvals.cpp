#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model22_karin1(sunindextype *rowvals){
    rowvals[0] = 1;
    rowvals[1] = 1;
    rowvals[2] = 0;
}