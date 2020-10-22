#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_Pritchard2002(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = GLCo;
    x_solver[1] = GLCi;
    x_solver[2] = ATP;
    x_solver[3] = G6P;
    x_solver[4] = ADP;
    x_solver[5] = F6P;
    x_solver[6] = F16bP;
    x_solver[7] = AMP;
    x_solver[8] = F26bP;
    x_solver[9] = DHAP;
    x_solver[10] = GAP;
    x_solver[11] = NAD;
    x_solver[12] = BPG;
    x_solver[13] = NADH;
    x_solver[14] = P3G;
    x_solver[15] = P2G;
    x_solver[16] = PEP;
    x_solver[17] = PYR;
    x_solver[18] = AcAld;
    x_solver[19] = CO2;
    x_solver[20] = EtOH;
    x_solver[21] = Glycerol;
    x_solver[22] = Glycogen;
    x_solver[23] = Trehalose;
    x_solver[24] = Succinate;
}