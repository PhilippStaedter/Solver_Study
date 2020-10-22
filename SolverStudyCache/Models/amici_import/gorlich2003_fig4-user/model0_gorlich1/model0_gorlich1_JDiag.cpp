#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_model0_gorlich1(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[2] = -83333333333.333328*dwdx2 + 83333333333.333328*dwdx3;
    JDiag[3] = 83333333333.333328*dwdx4 - 83333333333.333328*dwdx5;
    JDiag[4] = -83333333333.333328*dwdx6 + 83333333333.333328*dwdx7;
    JDiag[5] = 83333333333.333328*dwdx8 - 83333333333.333328*dwdx9;
    JDiag[6] = -55555555555.555557*dwdx10;
    JDiag[8] = 55555555555.555557*dwdx13;
    JDiag[9] = -83333333333.333328*dwdx14 - 83333333333.333328*dwdx15;
    JDiag[10] = -55555555555.555557*dwdx16 + 55555555555.555557*dwdx17;
    JDiag[11] = 55555555555.555557*dwdx18 - 55555555555.555557*dwdx19 - 55555555555.555557*dwdx20;
    JDiag[12] = -83333333333.333328*dwdx21 + 83333333333.333328*dwdx22;
}