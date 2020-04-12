#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_BIOMD0000000037(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 1.0*Glucose_sensor_inactivation_kG*Ya;
    dwdx[1] = 1.0*Photoreceptor_activation_IfrSfrPfr;
    dwdx[2] = 1.0*Photoreceptor_decay_kd;
    dwdx[3] = 1.0*Photoreceptor_inactivation_IrSrPr;
    dwdx[4] = 1.0*Transducer_activation_kia*Xi;
    dwdx[5] = 1.0*S_degradation_kd_s;
    dwdx[6] = -3.0*pow(S, 2)*V_formation_alpha2/pow(pow(S, 3) + 1, 2);
    dwdx[7] = -3.0*S_formation_alpha1*pow(V, 2)/pow(pow(V, 3) + 1, 2);
    dwdx[8] = 1.0*V_degradation_kd_v;
    dwdx[9] = 1.0*Transducer_inactivation_kai;
    dwdx[10] = 1.0*preS_formation_kx*prepreS;
    dwdx[11] = 1.0*Pr*Transducer_activation_kia;
    dwdx[12] = 1.0*Gluc*Glucose_sensor_inactivation_kG;
    dwdx[13] = 1.0*S_generation_ky*preS;
    dwdx[14] = 1.0*S_generation_ky*Ya;
    dwdx[15] = 1.0*Xa*preS_formation_kx;
}