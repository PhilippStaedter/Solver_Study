#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_hald(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 0.5*ACA0*kVap;
    w[1] = 2*V1adh1*V2adh1*sc*(-ACA*NADH/Keq + EtOH*NAD)/(ACA*EtOH*NAD*V2adh1/Kiadh1ACA + ACA*EtOH*NADH*V1adh1/(Keq*Kiadh1EtOH) + ACA*Kadh1NADH*V1adh1/Keq + ACA*Kadh1NADH*NAD*V1adh1/(Keq*Kiadh1NAD) + ACA*NADH*V1adh1/Keq + EtOH*Kadh1NAD*V2adh1 + EtOH*Kadh1NAD*NADH*V2adh1/Kiadh1NADH + EtOH*NAD*V2adh1 + Kadh1ACA*NADH*V1adh1/Keq + Kadh1EtOH*Kiadh1NAD*V2adh1 + Kadh1EtOH*NAD*V2adh1);
    w[2] = sc*(-pow(ADP, 2)*kakr + AMP*ATP*kakf);
    w[3] = sc*(-DHAP*GAP*Valdf/(Kaldeq*(DHAP*GAP*Valdf/(Kaldeq*Valdr) + DHAP*KaldGAP*Valdf/(Kaldeq*Valdr) + FBP*GAP/KaldIGAP + FBP + GAP*KaldDHAP*Valdf/(Kaldeq*Valdr) + KaldFBP)) + FBP*Valdf/(DHAP*GAP*Valdf/(Kaldeq*Valdr) + DHAP*KaldGAP*Valdf/(Kaldeq*Valdr) + FBP*GAP/KaldIGAP + FBP + GAP*KaldDHAP*Valdf/(Kaldeq*Valdr) + KaldFBP));
    w[4] = ACA*NAD*kOAc*sc;
    w[5] = 0.83199999999999996*ATP*Vatpconm*sc/(ATP + KatpconATP);
    w[6] = DHAP*HCN*kdhapf - DHAPCN*kdhapb;
    w[7] = sc*(-DPG*NADH*Vgapdh/(KgapdhGAP*KgapdhNAD*Kgapdheq*(1 + NADH/KgapdhNADH + NAD/KgapdhNAD)*(DPG/KgapdhDPG + GAP/KgapdhGAP + 1)) + GAP*NAD*Vgapdh/(KgapdhGAP*KgapdhNAD*(1 + NADH/KgapdhNADH + NAD/KgapdhNAD)*(DPG/KgapdhDPG + GAP/KgapdhGAP + 1)));
    w[8] = sc*(-Glc*Vglctr/(KglctrGlc*Vratio*(G6P*Glc/(KglctrGlc*KglctrIIG6P) + G6P/KglctrIG6P + Glc/KglctrGlc + (Glc0/KglctrGlc + 1)*(Glc*Pglctr/KglctrGlc + 1)/(Glc0*Pglctr/KglctrGlc + 1) + 1)) + Glc0*Vglctr/(KglctrGlc*Vratio*(Glc0/KglctrGlc + 1 + (Glc0*Pglctr/KglctrGlc + 1)*(G6P*Glc/(KglctrGlc*KglctrIIG6P) + G6P/KglctrIG6P + Glc/KglctrGlc + 1)/(Glc*Pglctr/KglctrGlc + 1))));
    w[9] = HCN0*kVap;
    w[10] = ATP*Glc*Vhkm*sc/(ATP*Glc + ATP*KhkGlc + Glc*KhkATP + KhkATP*KhkDGlc);
    w[11] = ACA0*HCN0*kacaf - kacab*lacto;
    w[12] = 0.69999999999999996*Pyr*V1PDC*sc/(KaPDC*(G6P/KxPDC + 1)/(G6P*bet/(KxPDC*alp) + 1) + Pyr*(G6P/(KxPDC*alp) + 1)/(G6P*bet/(KxPDC*alp) + 1));
    w[13] = ATP*pow(F6P, 2)*Vpfkm*sc/((ATP + KpfkmATP)*(pow(F6P, 2) + Kpfk*(1 + pow(ATP, 2)*kappapfk/pow(AMP, 2))));
    w[14] = sc*(-F6P*Vpgi/(Kpgieq*(F6P*KpgiG6P/KpgiF6P + G6P + KpgiG6P)) + G6P*Vpgi/(F6P*KpgiG6P/KpgiF6P + G6P + KpgiG6P));
    w[15] = ADP*PEP*Vpkm*sc/((ADP + KpkADP)*(KpkPEP + PEP));
    w[16] = HCN*Pyr*kpyrf - PyrCN*kpyrb;
    w[17] = 1.1499999999999999*Pyr*VPYRd*sc/(KmPYRd + Pyr);
    w[18] = ATP*G6P*kg6pst*sc;
    w[19] = sc*(DHAP*Vtimm/(DHAP + GAP*KtimDHAP/KtimGAP + KtimDHAP) - GAP*Vtimm/(Ktimeq*(DHAP + GAP*KtimDHAP/KtimGAP + KtimDHAP)));
    w[20] = sc*(ACA*dACA/Vratio - ACA0*dACA/Vratio);
    w[21] = sc*(EtOH*dEtOH/Vratio - EtOH0*dEtOH/Vratio);
    w[22] = sc*(Glyc*dGlyc/Vratio - Glyc0*dGlyc/Vratio);
    w[23] = sc*(HCN*dHCN/Vratio - HCN0*dHCN/Vratio);
    w[24] = sc*(OAc*dOAc/Vratio - OAc0*dOAc/Vratio);
    w[25] = DHAP*Vlpglycm*sc/(DHAP*(KlpglycNADH*(1 + NAD/KlpglycINAD)/NADH + 1) + KlpglycDHAP*(KlpglycINADH*(1 + NAD/KlpglycINAD)/NADH + 1));
    w[26] = sc*(ADP*DPG*klppepf - ATP*PEP*klppepr);
}