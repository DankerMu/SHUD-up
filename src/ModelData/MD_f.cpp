//  MD_f.cpp
//
//  Created by Lele Shu on 1/27/19.
//  Copyright © 2019 Lele Shu. All rights reserved.
//

#include "Model_Data.hpp"
#include "WaterBalanceDiag.hpp"
void Model_Data:: f_loop(double t){
    int i;
    for (i = 0; i < NumEle; i++) {
        if(lakeon && Ele[i].iLake > 0){
            /* Lake elements */
            Ele[i].updateLakeElement();
            fun_Ele_lakeVertical(i, t);
            qLakeEvap[Ele[i].iLake - 1] += qEleEvapo[i] / lake[Ele[i].iLake - 1].NumEleLake;
            qLakePrcp[Ele[i].iLake - 1] += qElePrep[i] / lake[Ele[i].iLake - 1].NumEleLake;
        }else{
            f_etFlux(i, t);
            /*DO INFILTRATION FRIST, then do LATERAL FLOW.*/
            /*========infiltration/Recharge Function==============*/
            Ele[i].updateElement(uYsf[i] , uYus[i] , uYgw[i] ); // step 1 update the kinf, kh, etc. for elements.
            fun_Ele_Infiltraion(i, t); // step 2 calculate the infiltration.
            fun_Ele_Recharge(i, t); // step 3 calculate the recharge.
        }
    }
    for (i = 0; i < NumEle; i++) {
        if(lakeon && Ele[i].iLake > 0){
            /* Lake elements */
            fun_Ele_lakeHorizon(i, t);
        }else{
            /*========surf/gw flow Function==============*/
            fun_Ele_surface(i, t);  // AFTER infiltration, do the lateral flux. ESP for overland flow.
            fun_Ele_sub(i, t);
        }
    } //end of for loop.
    for (i = 0; i < NumSegmt; i++) {
        fun_Seg_surface(RivSeg[i].iEle-1, RivSeg[i].iRiv-1, i);
        fun_Seg_sub(RivSeg[i].iEle-1, RivSeg[i].iRiv-1, i);
    }
    for (i = 0; i < NumRiv; i++) {
        Flux_RiverDown(t, i);
    }
    for (i = 0; i < NumLake; i++) {
        qLakeEvap[i] = min(qLakeEvap[i], qLakePrcp[i] + yLakeStg[i]);
        qLakeEvap[i] = max(0, qLakeEvap[i]);
    }
    /* Shared for both OpenMP and Serial, to update */
    PassValue();
}

void Model_Data::f_applyDY(double *DY, double t){
    double area;
    int isf, ius, igw;

    WaterBalanceDiag *wb = wbdiag;
    const bool want_quad = (wb != nullptr && wb->quadEnabled());
    const double *ic_raw = want_quad ? wb->icRawData() : nullptr;
    double p_rate_m3min = 0.0;
    double et_rate_m3min = 0.0;
    double qbc_rate_m3min = 0.0;
    double qss_rate_m3min = 0.0;
    double qedge_rate_m3min = 0.0;
    double noncons_edge_rate_m3min = 0.0;
    for (int i = 0; i < NumEle; i++) {
        isf = iSF; ius = iUS; igw = iGW;
        area = Ele[i].area;
        QeleSurfTot[i] = Qe2r_Surf[i];
        QeleSubTot[i] = Qe2r_Sub[i];
        for (int j = 0; j < 3; j++) {
            QeleSurfTot[i] += QeleSurf[i][j];
            QeleSubTot[i] += QeleSub[i][j];
            CheckNANij(QeleSurf[i][j], i, "QeleSurf[i][j]");
            CheckNANij(QeleSub[i][j], i, "QeleSub[i][j]");

            if (want_quad) {
                const int inabr = Ele[i].nabr[j] - 1;
                const int ilake = Ele[i].lakenabr[j] - 1;
                const double q_edge = QeleSurf[i][j] + QeleSub[i][j];
                if (inabr < 0 && ilake < 0 && !CS.CloseBoundary) {
                    qedge_rate_m3min += q_edge;
                }
                if (inabr >= 0) {
                    noncons_edge_rate_m3min += q_edge;
                }
            }
        }
        DY[i] = qEleNetPrep[i] - qEleInfil[i] + qEleExfil[i] - QeleSurfTot[i] / area - qEs[i];
        DY[ius] = qEleInfil[i] - qEleRecharge[i] - qEu[i] - qTu[i];
        DY[igw] = qEleRecharge[i] - qEleExfil[i] - QeleSubTot[i] / area - qEg[i] - qTg[i];
//        if(uYgw[i] > Ele[i].WetlandLevel ){ /* IF GW above surface, exfil from GW, OR, from UNSAT*/
//            DY[igw] += - qEleExfil[i];
//        }else{
//            if(uYus[i] > 0.){
//                DY[ius] += - qEleExfil[i];
//            }else{
//                DY[igw] += - qEleExfil[i];
//            }
//        }
        /* Boundary condition and Source/Sink */
        if(Ele[i].iBC == 0){
        }else if(Ele[i].iBC > 0){ // Fix head of GW.
            DY[igw] = 0;
        }else if(Ele[i].iBC < 0){ // Fix flux in GW
            DY[igw] += Ele[i].QBC / area;
        }

        if(Ele[i].iSS == 0){
        }else if(Ele[i].iSS > 0){ // SS in Landusrface
            DY[isf] += Ele[i].QSS / area;
        }else if(Ele[i].iSS < 0){ // SS in GW
            DY[igw] += Ele[i].QSS / area;
        }

        if (want_quad) {
            p_rate_m3min += qElePrep[i] * area;
            if (ic_raw != nullptr) {
                et_rate_m3min += (ic_raw[i] + qEs[i] + qEu[i] + qEg[i] + qTu[i] + qTg[i]) * area;
            } else {
                et_rate_m3min += (qEleE_IC[i] + qEs[i] + qEu[i] + qEg[i] + qTu[i] + qTg[i]) * area;
            }
            if (Ele[i].iBC < 0) {
                qbc_rate_m3min += Ele[i].QBC;
            }
            if (Ele[i].iSS != 0) {
                qss_rate_m3min += Ele[i].QSS;
            }
        }
        /* Convert with specific yield */
        DY[ius] /= Ele[i].Sy;
        DY[igw] /= Ele[i].Sy;
        
//        if(i+1==ID_ELE && DY[i]+uYsf[i] > 0.1 ){  // debug only
//            printf("%.3f, %d: %.2e, %.2e, %.2e | (%.2e, %.2e, %.2e, %.2e, %.2e )\n",
//                   t, i+1,
//                   uYsf[i], uYus[i], uYgw[i], qEleInfil[i], - qEleRecharge[i], - qEu[i],  - qTu[i], QeleSurfTot[i] / area);
//            printf("\n");
//        }
//        if(i +1== ID_ELE){  // debug only
//            printf("%.3f, %d: %f, %f, %f | (%f, %f, %f), %.2e, %.2e\n",
//                   t, i+1,
//                   DY[i], DY[ius], DY[igw], QeleSurf[i][0] / area, QeleSurf[i][1] / area, QeleSurf[i][2] / area,
//                   -QeleSurfTot[i] / area, qEleRecharge[i]);
//            printf("\n");
//        }
        if(Ele[i].iLake > 0){
            DY[i] = 0.;
            DY[ius] = 0.;
            DY[igw] = 0.;
        }
#ifdef DEBUG
        CheckNANi(DY[i], i, "DY[i] (Model_Data::f_applyDY)");
        CheckNANi(DY[ius], i, "DY[ius] (Model_Data::f_applyDY)");
        CheckNANi(DY[igw], i, "DY[igw] (Model_Data::f_applyDY)");
#endif
    }
    for (int i = 0; i < NumRiv; i++) {
        if(Riv[i].BC > 0){
//            Newmann condition.
            DY[iRIV] = 0.;
        }else{
            DY[iRIV] = (- QrivUp[i] - QrivSurf[i] - QrivSub[i] - QrivDown[i] + Riv[i].qBC) / Riv[i].Length; // dA on CS
            if(DY[iRIV] < -1. * Riv[i].u_CSarea){ /* The negative dA cannot larger then Availalbe Area. */
                DY[iRIV] = -1. * Riv[i].u_CSarea;
            }
            DY[iRIV] = fun_dAtodY(DY[iRIV], Riv[i].u_topWidth, Riv[i].bankslope);
            
//            if(i+1 == ID_RIV && fabs(DY[iRIV])> 0.001){
//                printf("%d, up:%f, down:%f, sub:%f, surf:%f, dy=%f\n", i+1,
//                       - QrivUp[i] / Riv[i].u_TopArea, - QrivDown[i] / Riv[i].u_TopArea,
//                       - QrivSub[i] / Riv[i].u_TopArea, - QrivSurf[i] / Riv[i].u_TopArea,
//                       DY[iRIV]);
//                i=i;
//            }
        }
#ifdef DEBUG
        CheckNANi(DY[i + 3 * NumEle], i, "DY[i] of river (Model_Data::f_applyDY)");
#endif
    }
    for(int i = 0; i < NumLake; i++){
//        DY[i + 3 * NumEle + NumRiv]
        DY[iLAKE] = qLakePrcp[i] - qLakeEvap[i]  +
                    (QLakeRivIn[i] - QLakeRivOut[i] + QLakeSub[i] + QLakeSurf[i] ) / y2LakeArea[i] ;        
//        if(fabs(DY[iLAKE]) > 1.0e-4){
//            printf("%f: %g + %g\n",t, yLakeStg[i], DY[iLAKE]);
//            i=i;
//        }
#ifdef DEBUG
        CheckNANi(DY[iLAKE], i, "DY[i] of LAKE (Model_Data::f_applyDY)");
#endif
    }

    if (want_quad) {
        double qout_rate_m3min = 0.0;
        for (int i = 0; i < NumRiv; i++) {
            if (Riv[i].toLake >= 0) {
                continue;
            }
            if (Riv[i].down < 0) {
                qout_rate_m3min += QrivDown[i];
            }
            if (Riv[i].BC < 0) {
                qbc_rate_m3min += Riv[i].qBC;
            }
        }
        wb->updateQuadRates(t,
                            p_rate_m3min,
                            et_rate_m3min,
                            qout_rate_m3min,
                            qedge_rate_m3min,
                            qbc_rate_m3min,
                            qss_rate_m3min,
                            noncons_edge_rate_m3min);
    }
}

void Model_Data::PassValue(){
    int i, ie, ir;
    for (i = 0; i < NumRiv; i++) {
        QrivSurf[i] = 0.;
        QrivSub[i] = 0.;
        QrivUp[i] = 0.;
    }
    for (i = 0; i < NumEle; i++) {
        Qe2r_Surf[i] = 0.;
        Qe2r_Sub[i] = 0.;
    }
    for (i = 0; i < NumSegmt; i++) {
        ie = RivSeg[i].iEle-1;
        ir = RivSeg[i].iRiv-1;
        QrivSurf[ir] += QsegSurf[i]; // Positive from River to Element
        QrivSub[ir] += QsegSub[i];
        Qe2r_Surf[ie] += -QsegSurf[i]; // Positive from Element to River
        Qe2r_Sub[ie] += -QsegSub[i];
    }
    for (i = 0; i < NumRiv; i++) {
        if(iDownStrm >= 0 && Riv[i].toLake <= 0){
            QrivUp[iDownStrm] += - QrivDown[i];
        }
    }
    //    for (i = 0; i < NumEle; i++) { /*Check flux A->B  = Flux B->A*/
    //        for (j = 0; j < 3; j++) {
    //            inabr = Ele[i].nabr[j] - 1;
    //            if (inabr >= 0) {
    //                jnabr = Ele[i].nabrToMe[j];
    //                if(Ele[inabr].iupdSF[jnabr]){
    //                    //void
    //                }else{
    //                    QeleSurf[inabr][jnabr] = - QeleSurf[i][j];
    //                    Ele[inabr].iupdSF[jnabr] = 1;
    //                    QeleSub[inabr][jnabr] = - QeleSub[i][j];
    //                    Ele[inabr].iupdGW[jnabr] = 1;
    //                }
    //            }
    //        }
    //    }
}

void Model_Data::applyBCSS(double *DY, int i){
    if(Ele[i].iBC > 0){ // Fix head of GW.
        DY[iGW] = 0;
    }else if(Ele[i].iBC < 0){ // Fix flux in GW
        DY[iGW] += Ele[i].QBC / Ele[i].area;
    }else{}
    
    if(Ele[i].iSS > 0){ // SS in Landusrface
        DY[iSF] += Ele[i].QSS / Ele[i].area;
    }else if(Ele[i].iSS < 0){ // SS in GW
        DY[iGW] += Ele[i].QSS / Ele[i].area;
    }else{}
}
