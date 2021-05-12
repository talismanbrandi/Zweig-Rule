/*
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <limits>
#include "myObservables.h"


myObservables::myObservables(const StandardModel& SM_i)
: ThObservable(SM_i),
  VCKM(3,3,0.),
  my_model(static_cast<const myModel*> (&SM_i))
{}

myObservables::~myObservables()
{}

void myObservables::updateParameters()
{
    T = my_model->getT();
    C = my_model->getC();
    Tp = T;
    Cp = C;
    TD = T;
    CD = C;
    kappa = my_model->getkappa();
    kappa_p = my_model->getkappa_p();
    // kappa_p = kappa;
    D = my_model->getD();
    K = my_model->getK();
    Kp = my_model->getKp();
    delta_0 = my_model->getdelta_0();
    delta_0_p = my_model->getdelta_0_p();
    eps_d = my_model->geteps_d();
    delta_12 = my_model->getdelta_12();
    delta_12s = my_model->getdelta_12()*(1. - eps_d);
    delta_1 = my_model->getdelta_1();
    delta_1s = my_model->getdelta_1()*(1. - eps_d);
    phi = my_model->getphi();

    cos2phi = cos(2. * phi);
    sin2phi = sin(2. * phi);

    ii = gslpp::complex(0., 1., 0);

    eIdelta_0 = exp(ii * delta_0);
    eIdelta_0_p = exp(ii * delta_0_p);
    eIdelta_1 = exp(ii * delta_1);
    eIdelta_12 = exp(ii * delta_12);
    eIdelta_1s = exp(ii * delta_1s);
    eIdelta_12s = exp(ii * delta_12s);

    /* CPV Parameters */

    Delta_3 = my_model->getDelta_3_T() * T;
    Delta_4 = my_model->getDelta_4_T() * T;
    P = my_model->getP_T() * T;

    GF = SM.getGF();
    VCKM = SM.getVCKM();

    V_ud = VCKM(0,0);
    V_us = VCKM(0,1);
    V_cd = VCKM(1,0);
    V_cs = VCKM(1,1);
    V_ub = VCKM(0,2);
    V_cb = VCKM(1,2);

    ckm_scs = (V_us * (V_cs).conjugate() - V_ud * (V_cd).conjugate())/2.;
    ckm_scs_bar = (ckm_scs).conjugate();

    ckm_cf = V_ud * (V_cs).conjugate();
    ckm_cf_bar = (ckm_cf).conjugate();

    ckm_dcs = - V_us * (V_cd).conjugate(); /*-ve sign for Buccella.*/
    ckm_dcs_bar = (ckm_dcs).conjugate();

    ckm_b = (V_ub * (V_cb).conjugate());
    ckm_b_bar = (ckm_b).conjugate();

    ckm_CPV = (V_us * (V_cs).conjugate() + V_ud * (V_cd).conjugate())/2.;
    ckm_CPV_bar = ckm_CPV.conjugate();

    //std::cout << " ckm_scs " << ckm_scs << " ckm_scs_bar " << ckm_scs_bar << " ckm_CPV " << ckm_CPV << " ckm_CPV_bar " << ckm_CPV_bar << " R1 " << ckm_CPV/ckm_scs << " R2 " <<  ckm_CPV_bar/ckm_scs_bar << std::endl;

    MD0 = 1.86484;
    MDP = 1.86961;
    MDS = 1.96830;
    MPIP = 0.13957018;
    MPIM = 0.13957018;
    MKP = 0.493677;
    MKM = 0.493677;
    MPI0 = 0.1349766;
    MKS = 0.497614;
    MKL = 0.497614;
    MET = 0.547862;
    METP = 0.95778;

    GD0 = 1./0.4101 * 6.58211928e-13;
    GDP = 1./1.040 * 6.58211928e-13;
    GDS = 1./0.500 * 6.58211928e-13;

    sqrt_10 = sqrt(10.);
    sqrt_30 = sqrt(30.);
    sqrt_3_5 = sqrt(3./5.);
}

double myObservables::PS(const double M, const double m1, const double m2) const
{
    double pc = sqrt( (M*M - (m1 + m2)*(m1 + m2))*(M*M - (m1 - m2)*(m1 - m2)) )/2.0/M;
    return ( pc/(8.*M_PI*M*M)*GF*GF/2.0 );
}

double myObservables::BR_bar(double amp, double amp_bar){
    return (amp + amp_bar)/2.;
}

double myObservables::aCP(double amp, double amp_bar){
    return (amp - amp_bar)/(amp + amp_bar) * 1.e2;
}

/*******************************************************************************
 * Observables                                                                 *
 * ****************************************************************************/

 PENG_T::PENG_T(const StandardModel& SM_i)
 : myObservables(SM_i)
 {}

 double PENG_T::computeThValue()
 {
    updateParameters();
    return (P + T)/T;
 }

 /************************************************
 *  SCS Decays
 ************************************************/

BRpppm::BRpppm(const StandardModel& SM_i, int obstype_i)
: myObservables(SM_i)
{   obstype = obstype_i; }

double BRpppm::computeThValue()
{
    updateParameters();

    gslpp::complex Apppm = (T - 2./3. * C) * (-3./10. * (eIdelta_0 + eIdelta_0_p) + (-3./10. * cos2phi + 3./4./sqrt_10 * sin2phi) * (eIdelta_0_p - eIdelta_0)) - (T + C) * 2./5.;
    gslpp::complex ApppmB = ckm_scs_bar * Apppm;
    Apppm = ckm_scs * Apppm;

    gslpp::complex Bpppm = (P + T + Delta_3) * (1./2. * (eIdelta_0_p + eIdelta_0) + (eIdelta_0_p - eIdelta_0) * (-1./6. * cos2phi - 7./4./sqrt_10 * sin2phi)) + (T + C) * (-3./20. * (eIdelta_0_p + eIdelta_0) + 3./10. + (1./60. * cos2phi + 1./2./sqrt_10 * sin2phi) * (eIdelta_0_p - eIdelta_0)) + Delta_4 * (eIdelta_0_p - eIdelta_0) * (-1./3. * cos2phi - 1./4./sqrt_10 * sin2phi);
    gslpp::complex BpppmB = ckm_CPV_bar * Bpppm;
    Bpppm = ckm_CPV * Bpppm;

    if (obstype == 1) return PS(MD0,MPIP,MPIM) * BR_bar((Apppm + Bpppm).abs2(), (ApppmB + BpppmB).abs2()) / GD0;
    else if (obstype == 2) return aCP((Apppm + Bpppm).abs2(), (ApppmB + BpppmB).abs2());
    else return 0.;
}

BRp0p0::BRp0p0(const StandardModel& SM_i, int obstype_i)
: myObservables(SM_i)
{   obstype = obstype_i; }

double BRp0p0::computeThValue()
{
    updateParameters();

    gslpp::complex Ap0p0 = (T - 2./3. * C) * (-3./10. * (eIdelta_0 + eIdelta_0_p) + (-3./10. * cos2phi + 3./4./sqrt_10 * sin2phi) * (eIdelta_0_p - eIdelta_0)) + (T + C) * 3./5.;
    gslpp::complex Ap0p0B = ckm_scs_bar * Ap0p0;
    Ap0p0 = ckm_scs * Ap0p0;

    gslpp::complex Bpppm = (P + T + Delta_3) * (1./2. * (eIdelta_0_p + eIdelta_0) + (eIdelta_0_p - eIdelta_0) * (-1./6. * cos2phi - 7./4./sqrt_10 * sin2phi)) + (T + C) * (-3./20. * (eIdelta_0_p + eIdelta_0) + 3./10. + (1./60. * cos2phi + 1./2./sqrt_10 * sin2phi) * (eIdelta_0_p - eIdelta_0)) + Delta_4 * (eIdelta_0_p - eIdelta_0) * (-1./3. * cos2phi - 1./4./sqrt_10 * sin2phi);
    gslpp::complex Bp0p0 = Bpppm - (T + C);
    gslpp::complex Bp0p0B = ckm_CPV_bar * Bp0p0;
    Bp0p0 = ckm_CPV * Bp0p0;

    if (obstype == 1) return PS(MD0,MPI0,MPI0) * BR_bar((Ap0p0 + Bp0p0).abs2(), (Ap0p0B + Bp0p0B).abs2())/2. / GD0;
    else if (obstype == 2) return aCP((Ap0p0 + Bp0p0).abs2(), (Ap0p0B + Bp0p0B).abs2());
    else return 0.;
}

BRkpkm::BRkpkm(const StandardModel& SM_i, int obstype_i)
: myObservables(SM_i)
{   obstype = obstype_i; }

double BRkpkm::computeThValue()
{
    updateParameters();

    gslpp::complex Akpkm = (T - 2./3. * C) * (3./20. * (eIdelta_0 + eIdelta_0_p) + (3./20. * cos2phi + 3./4./sqrt_10 * sin2phi) * (eIdelta_0_p - eIdelta_0) + 3./10. * eIdelta_1) + (T + C) * 2./5.;
    gslpp::complex AkpkmB = ckm_scs_bar * Akpkm;
    Akpkm = ckm_scs * Akpkm;

    gslpp::complex Bkpkm = (P + T + Delta_3) * (1./4. * (eIdelta_0_p + eIdelta_0) + (eIdelta_0_p - eIdelta_0) * (-5./12. * cos2phi + 1./4./sqrt_10 * sin2phi) + 1./2. * eIdelta_1) + (T + C) * (-1./20. * (eIdelta_0_p + eIdelta_0) + 3./10. + (7./60. * cos2phi) * (eIdelta_0_p - eIdelta_0) - 1./5. * eIdelta_1) + Delta_4 * (1./4. * (eIdelta_0_p + eIdelta_0) + (eIdelta_0_p - eIdelta_0) * (-1./12. * cos2phi + 3./4./sqrt_10 * sin2phi) - 1./2. * eIdelta_1);
    gslpp::complex BkpkmB = ckm_CPV_bar * Bkpkm;
    Bkpkm = ckm_CPV * Bkpkm;

    if (obstype == 1) return PS(MD0,MKP,MKM) * BR_bar((Akpkm + Bkpkm).abs2(), (AkpkmB + BkpkmB).abs2()) / GD0;
    else if (obstype == 2) return aCP((Akpkm + Bkpkm).abs2(), (AkpkmB + BkpkmB).abs2());
    else return 0.;
}

BRkSkS::BRkSkS(const StandardModel& SM_i, int obstype_i)
: myObservables(SM_i)
{   obstype = obstype_i; }

double BRkSkS::computeThValue()
{
    updateParameters();

    gslpp::complex Ak0k0bar = (T - 2./3. * C) * (3./20. * (eIdelta_0 + eIdelta_0_p) + (3./20. * cos2phi + 3./4./sqrt_10 * sin2phi) * (eIdelta_0_p - eIdelta_0) - 3./10. * eIdelta_1);
    gslpp::complex Ak0k0barB = ckm_scs_bar * Ak0k0bar;
    Ak0k0bar = ckm_scs * Ak0k0bar;

    gslpp::complex Bk0k0bar = (P + T + Delta_3) * (1./4. * (eIdelta_0_p + eIdelta_0) + (eIdelta_0_p - eIdelta_0) * (-5./12. * cos2phi + 1./4./sqrt_10 * sin2phi) - 1./2. * eIdelta_1) + (T + C) * (-1./20. * (eIdelta_0_p + eIdelta_0) - 1./10. + (7./60. * cos2phi) * (eIdelta_0_p - eIdelta_0) + 1./5. * eIdelta_1) + Delta_4 * (1./4. * (eIdelta_0_p + eIdelta_0) + (eIdelta_0_p - eIdelta_0) * (-1./12. * cos2phi + 3./4./sqrt_10 * sin2phi) + 1./2. * eIdelta_1);
    gslpp::complex Bk0k0barB = ckm_CPV_bar * Bk0k0bar;
    Bk0k0bar = ckm_CPV * Bk0k0bar;

    if (obstype == 1) return PS(MD0,MKS,MKS) * BR_bar((Ak0k0bar + Bk0k0bar).abs2(), (Ak0k0barB + Bk0k0barB).abs2())/2. / GD0;
    else if (obstype == 2) return aCP((Ak0k0bar + Bk0k0bar).abs2(), (Ak0k0barB + Bk0k0barB).abs2());
    else return 0.;
}

BRetet::BRetet(const StandardModel& SM_i, int obstype_i)
: myObservables(SM_i)
{   obstype = obstype_i; }

double BRetet::computeThValue()
{
    updateParameters();

    gslpp::complex Aet8et8 = (T - 2./3. * C) * (3./10. * (eIdelta_0 + eIdelta_0_p) + (3./10. * cos2phi + 3./4./sqrt_10 * sin2phi) * (eIdelta_0_p - eIdelta_0)) - (T + C) * 3./5.;
    gslpp::complex AetetB = ckm_scs_bar * Aet8et8;
    gslpp::complex Aetet = ckm_scs * Aet8et8;

    if (obstype == 1) return PS(MD0,MET,MET) * BR_bar((Aetet).abs2(), (AetetB).abs2())/2. / GD0;
    else return 0.;
}

BRp0et::BRp0et(const StandardModel& SM_i, int obstype_i)
: myObservables(SM_i)
{   obstype = obstype_i; }

double BRp0et::computeThValue()
{
    updateParameters();

    gslpp::complex Ap0et8 = sqrt(3.)/5. * ((T - 2./3. * C) * eIdelta_1 - (T + C));
    gslpp::complex Ap0etB = ckm_scs_bar * Ap0et8;
    gslpp::complex Ap0et = ckm_scs * Ap0et8;

    gslpp::complex Bp0et8 = 1./sqrt(3.) * ((P + T + Delta_3 - Delta_4) * eIdelta_1 - (T + C) * (2./5. * eIdelta_1 + 3./5.));
    gslpp::complex Bp0etB = ckm_CPV_bar * Bp0et8;
    gslpp::complex Bp0et = ckm_CPV * Bp0et8;

    if (obstype == 1) return PS(MD0,MPI0,MET) * BR_bar((Ap0et + Bp0et).abs2(), (Ap0etB + Bp0etB).abs2()) / GD0;
    else if (obstype == 2) return aCP((Ap0et + Bp0et).abs2(), (Ap0etB + Bp0etB).abs2());
    else return 0.;
}

BRkpkS::BRkpkS(const StandardModel& SM_i, int obstype_i)
: myObservables(SM_i)
{   obstype = obstype_i; }

double BRkpkS::computeThValue()
{
    updateParameters();

    // gslpp::complex Akpk0bar = 1./5. * (2. * T - 3. * C + D) * eIdelta_1 + 3./5. * (T + C);
    gslpp::complex Akpk0bar = 1./5. * (2. * T - 3. * C + D) * eIdelta_1 + 3./5. * (T + C);
    gslpp::complex Akpk0barB = ckm_scs_bar * Akpk0bar;
    Akpk0bar = ckm_scs * Akpk0bar;

    gslpp::complex Bkpk0bar = (P + T + Delta_3 - Delta_4 - 1./5. * D) * eIdelta_1 + 1./5. * (T + C) * (1. - 1. * eIdelta_1);
    gslpp::complex Bkpk0barB = ckm_CPV_bar * Bkpk0bar;
    Bkpk0bar = ckm_CPV * Bkpk0bar;

    if (obstype == 1) return PS(MDP,MKP,MKS) * BR_bar((Akpk0bar + Bkpk0bar).abs2(), (Akpk0barB + Bkpk0barB).abs2())/2. / GDP;
    else if (obstype == 2) return aCP((Akpk0bar + Bkpk0bar).abs2(), (Akpk0barB + Bkpk0barB).abs2());
    else return 0.;
}

BRppp0::BRppp0(const StandardModel& SM_i, int obstype_i)
: myObservables(SM_i)
{   obstype = obstype_i; }

double BRppp0::computeThValue()
{
    updateParameters();

    gslpp::complex Appp0 = (T + C) / sqrt(2.);
    gslpp::complex Appp0B = ckm_scs_bar * Appp0;
    Appp0 = ckm_scs * Appp0;

    gslpp::complex Bppp0 = -Appp0;
    gslpp::complex Bppp0B = ckm_CPV_bar * Bppp0;
    Bppp0 = ckm_CPV * Bppp0;

    if (obstype == 1) return PS(MDP,MPIP,MPI0) * BR_bar((Appp0 + Bppp0).abs2(), (Appp0B + Bppp0B).abs2()) / GDP;
    else return 0.;
}

BRppet::BRppet(const StandardModel& SM_i, int obstype_i)
: myObservables(SM_i)
{   obstype = obstype_i; }

double BRppet::computeThValue()
{
    updateParameters();

    //gslpp::complex Appet8 = sqrt(2./3.)/5. * (2. * T - 3. * C + D) * eIdelta_1 - 3./5. * sqrt(3./2.) * (T + C);
    gslpp::complex Appet8 = sqrt(2./3.)/5. * (2. * T - 3. * C + D) * eIdelta_1 - 3./5. * sqrt(3./2.) * (T + C);
    gslpp::complex AppetB = ckm_scs_bar * Appet8;
    gslpp::complex Appet = ckm_scs *Appet8;

    gslpp::complex Bppet8 = sqrt(2./3.) * (P + T + Delta_3 - Delta_4 - 1./5. * D - 1./5. * (T + C)) * eIdelta_1 - sqrt(6.)/10. * (T + C);
    gslpp::complex BppetB = ckm_CPV_bar * Bppet8;
    gslpp::complex Bppet = ckm_CPV * Bppet8;

    if (obstype == 1) return PS(MDP,MPIP,MET) * BR_bar((Appet + Bppet).abs2(), (AppetB + BppetB).abs2()) / GDP;
    else if (obstype == 2) return aCP((Appet + Bppet).abs2(), (AppetB + BppetB).abs2());
    else return 0.;
}

BRDsppkS::BRDsppkS(const StandardModel& SM_i, int obstype_i)
: myObservables(SM_i)
{   obstype = obstype_i; }

double BRDsppkS::computeThValue()
{
    updateParameters();

    gslpp::complex ADsppk0 = -1./5. * (2. * T - 3. * C + D - Kp) * eIdelta_12s - 3./5. * (T + C);
    gslpp::complex ADsppk0B = ckm_scs_bar * ADsppk0;
    ADsppk0 = ckm_scs * ADsppk0;

    // gslpp::complex BDsppk0 = (P + T + Delta_3 - Delta_4 - 1./5. * (D + T + C)) * eIdelta_12s + 1./5. * (T + C);
    gslpp::complex BDsppk0 = -(P + T + Delta_3 - Delta_4 - 1./5. * (D + T + C)) * eIdelta_12s - 1./5. * (T + C);
    gslpp::complex BDsppk0B = ckm_CPV_bar * BDsppk0;
    BDsppk0 = ckm_CPV * BDsppk0;

    if (obstype == 1) return PS(MDS,MPIP,MKS) * BR_bar((ADsppk0 + BDsppk0).abs2(), (ADsppk0B + BDsppk0B).abs2())/2. / GDS;
    else if (obstype == 2) return aCP((ADsppk0 + BDsppk0).abs2(), (ADsppk0B + BDsppk0B).abs2());
    else return 0.;
}

BRDsp0kp::BRDsp0kp(const StandardModel& SM_i, int obstype_i)
: myObservables(SM_i)
{   obstype = obstype_i; }

double BRDsp0kp::computeThValue()
{
    updateParameters();

    gslpp::complex ADsp0kp = -1./5./sqrt(2.) * (2. * T - 3. * C + D - Kp) * eIdelta_12s + 2./5./sqrt(2.) * (T + C);
    gslpp::complex ADsp0kpB = ckm_scs_bar * ADsp0kp;
    ADsp0kp = ckm_scs * ADsp0kp;

    // gslpp::complex BDsp0kp = sqrt(1./2.) * (P + T + Delta_3 - Delta_4 - 1./5. * (D + T + C)) * eIdelta_12s - 2./5. * sqrt(2.) * (T + C);
    gslpp::complex BDsp0kp = - sqrt(1./2.) * (P + T + Delta_3 - Delta_4 - 1./5. * (D + T + C)) * eIdelta_12s + 2./5. * sqrt(2.) * (T + C);
    gslpp::complex BDsp0kpB = ckm_CPV_bar * BDsp0kp;
    BDsp0kp = ckm_CPV * BDsp0kp;

    if (obstype == 1) return PS(MDS,MPI0,MKP) * BR_bar((ADsp0kp + BDsp0kp).abs2(), (ADsp0kpB + BDsp0kpB).abs2()) / GDS;
    else if (obstype == 2) return aCP((ADsp0kp + BDsp0kp).abs2(), (ADsp0kpB + BDsp0kpB).abs2());
    else return 0.;
}

BRDskpet::BRDskpet(const StandardModel& SM_i, int obstype_i)
: myObservables(SM_i)
{   obstype = obstype_i; }

double BRDskpet::computeThValue()
{
    updateParameters();

    gslpp::complex ADskpet8 = 1./5./sqrt(6.) * (2. * T - 3. * C + D - Kp) * eIdelta_12s - 2. * sqrt(6.)/5. * (T + C);
    gslpp::complex ADskpetB = ckm_scs_bar * ADskpet8;
    gslpp::complex ADskpet = ckm_scs * ADskpet8;

    if (obstype == 1) return PS(MDS,MKP,MET) * BR_bar((ADskpet).abs2(), (ADskpetB).abs2()) / GDS;
    else return 0.;
}

/************************************************
*  CF & DCS Decays
************************************************/

BRppkm::BRppkm(const StandardModel& SM_i, int obstype_i)
: myObservables(SM_i)
{   obstype = obstype_i; }

double BRppkm::computeThValue()
{
    updateParameters();

    gslpp::complex Appkm = 1./5. * (3. * Tp - 2. * Cp - K) * eIdelta_12 + 2./5. * (Tp + Cp + kappa);
    gslpp::complex AppkmB = ckm_cf_bar * Appkm;
    Appkm = ckm_cf * Appkm;

    if (obstype == 1) return PS(MD0,MPIP,MKM) * BR_bar((Appkm).abs2(), (AppkmB).abs2()) / GD0;
    else return 0.;
}

BRp0kS::BRp0kS(const StandardModel& SM_i, int obstype_i)
: myObservables(SM_i)
{   obstype = obstype_i; }

double BRp0kS::computeThValue()
{
    updateParameters();

    gslpp::complex Ap0k0bar = -1./5./sqrt(2.) * (3. * Tp - 2. * Cp - K) * eIdelta_12 + 3./5./sqrt(2.) * (Tp + Cp + kappa);
    gslpp::complex Ap0k0barB = ckm_cf_bar * Ap0k0bar;
    Ap0k0bar = ckm_cf * Ap0k0bar;

    // gslpp::complex Ap0k0 = 1./5./sqrt(2.) * (3. * TD - 2. * CD + K) * eIdelta_12 - 3./5./sqrt(2.) * (TD + CD + kappa);/* ERROR */
    gslpp::complex Ap0k0 = 1./5./sqrt(2.) * (3. * TD - 2. * CD + K) * eIdelta_12 - 3./5./sqrt(2.) * (TD + CD + kappa_p);/* ERROR */
    gslpp::complex Ap0k0B = ckm_dcs_bar * Ap0k0;
    Ap0k0 = ckm_dcs * Ap0k0;

    if (obstype == 1) return PS(MD0,MPI0,MKS) * BR_bar((Ap0k0 - Ap0k0bar).abs2(), (Ap0k0B - Ap0k0barB).abs2())/2. / GD0;
    else return 0.;
}

BRp0kL::BRp0kL(const StandardModel& SM_i, int obstype_i)
: myObservables(SM_i)
{   obstype = obstype_i; }

double BRp0kL::computeThValue()
{
    updateParameters();

    gslpp::complex Ap0k0bar = -1./5./sqrt(2.) * (3. * Tp - 2. * Cp - K) * eIdelta_12 + 3./5./sqrt(2.) * (Tp + Cp + kappa);
    gslpp::complex Ap0k0barB = ckm_cf_bar * Ap0k0bar;
    Ap0k0bar = ckm_cf * Ap0k0bar;

    // gslpp::complex Ap0k0 = 1./5./sqrt(2.) * (3. * TD - 2. * CD + K) * eIdelta_12 - 3./5./sqrt(2.) * (TD + CD + kappa);
    gslpp::complex Ap0k0 = 1./5./sqrt(2.) * (3. * TD - 2. * CD + K) * eIdelta_12 - 3./5./sqrt(2.) * (TD + CD + kappa_p);
    gslpp::complex Ap0k0B = ckm_dcs_bar * Ap0k0;
    Ap0k0 = ckm_dcs * Ap0k0;

    if (obstype == 1) return PS(MD0,MPI0,MKS) * BR_bar((Ap0k0 + Ap0k0bar).abs2(), (Ap0k0B + Ap0k0barB).abs2())/2. / GD0;
    else return 0.;
}

BRkSet::BRkSet(const StandardModel& SM_i, int obstype_i)
: myObservables(SM_i)
{   obstype = obstype_i; }

double BRkSet::computeThValue()
{
    updateParameters();

    gslpp::complex Ak0baret8 = -1./5./sqrt(6.) * (3. * Tp - 2. * Cp - K) * eIdelta_12 + 3./5./sqrt(6.) * (Tp + Cp + kappa);
    gslpp::complex Ak0baretB = ckm_cf_bar * Ak0baret8;
    gslpp::complex Ak0baret = ckm_cf * Ak0baret8;

    // gslpp::complex Ak0et8 = 1./5./sqrt(6.) * (3. * TD - 2. * CD + K) * eIdelta_12 - sqrt(3.)/5./sqrt(2.) * (TD + CD + kappa);
    gslpp::complex Ak0et8 = 1./5./sqrt(6.) * (3. * TD - 2. * CD + K) * eIdelta_12 - sqrt(3.)/5./sqrt(2.) * (TD + CD + kappa_p);
    gslpp::complex Ak0etB = ckm_dcs_bar * Ak0et8;
    gslpp::complex Ak0et = ckm_dcs * Ak0et8;

    if (obstype == 1) return PS(MD0,MKS,MET) * BR_bar((Ak0et - Ak0baret).abs2(), (Ak0etB - Ak0baretB).abs2())/2. / GD0;
    else return 0.;
}

BRppkS::BRppkS(const StandardModel& SM_i, int obstype_i)
: myObservables(SM_i)
{   obstype = obstype_i; }

double BRppkS::computeThValue()
{
    updateParameters();

    gslpp::complex Appk0bar = (Tp + Cp + kappa);
    gslpp::complex Appk0barB = ckm_cf_bar * Appk0bar;
    Appk0bar = ckm_cf * Appk0bar;

    // gslpp::complex Appk0 = 1./5. * (2. * TD - 3. * CD + D - Kp) * eIdelta_12 - 2./5. * (TD + CD + kappa);
    gslpp::complex Appk0 = 1./5. * (2. * TD - 3. * CD + D - Kp) * eIdelta_12 - 2./5. * (TD + CD + kappa_p);
    gslpp::complex Appk0B = ckm_dcs_bar * Appk0;
    Appk0 = ckm_dcs * Appk0;

    if (obstype == 1) return PS(MDP,MPIP,MKS) * BR_bar((Appk0 - Appk0bar).abs2(), (Appk0B - Appk0barB).abs2())/2. / GDP;
    else return 0.;
}

BRppkL::BRppkL(const StandardModel& SM_i, int obstype_i)
: myObservables(SM_i)
{   obstype = obstype_i; }

double BRppkL::computeThValue()
{
    updateParameters();

    gslpp::complex Appk0bar = (Tp + Cp + kappa);
    gslpp::complex Appk0barB = ckm_cf_bar * Appk0bar;
    Appk0bar = ckm_cf * Appk0bar;

    // gslpp::complex Appk0 = 1./5. * (2. * TD - 3. * CD + D - Kp) * eIdelta_12 - 2./5. * (TD + CD + kappa);
    gslpp::complex Appk0 = 1./5. * (2. * TD - 3. * CD + D - Kp) * eIdelta_12 - 2./5. * (TD + CD + kappa_p);
    gslpp::complex Appk0B = ckm_dcs_bar * Appk0;
    Appk0 = ckm_dcs * Appk0;

    if (obstype == 1) return PS(MDP,MPIP,MKL) * BR_bar((Appk0 + Appk0bar).abs2(), (Appk0B + Appk0barB).abs2())/2. / GDP;
    else return 0.;
}

BRDskpkS::BRDskpkS(const StandardModel& SM_i, int obstype_i)
: myObservables(SM_i)
{   obstype = obstype_i; }

double BRDskpkS::computeThValue()
{
    updateParameters();

    gslpp::complex ADskpk0bar = -1./5. * (2. * Tp - 3. * Cp + D) * eIdelta_1s + 2./5. * (Tp + Cp + kappa);
    gslpp::complex ADskpk0barB = ckm_cf_bar * ADskpk0bar;
    ADskpk0bar = ckm_cf * ADskpk0bar;

    // gslpp::complex ADskpk0 = -(TD + CD + kappa);
    gslpp::complex ADskpk0 = -(TD + CD + kappa_p);
    gslpp::complex ADskpk0B = ckm_dcs_bar * ADskpk0;
    ADskpk0 = ckm_dcs * ADskpk0;

    if (obstype == 1) return PS(MDS,MKP,MKS) * BR_bar((ADskpk0 - ADskpk0bar).abs2(), (ADskpk0B - ADskpk0barB).abs2())/2. / GDS;
    else return 0.;
}

BRDskpkL::BRDskpkL(const StandardModel& SM_i, int obstype_i)
: myObservables(SM_i)
{   obstype = obstype_i; }

double BRDskpkL::computeThValue()
{
    updateParameters();

    gslpp::complex ADskpk0bar = -1./5. * (2. * Tp - 3. * Cp + D) * eIdelta_1s + 2./5. * (Tp + Cp + kappa);
    gslpp::complex ADskpk0barB = ckm_cf_bar * ADskpk0bar;
    ADskpk0bar = ckm_cf * ADskpk0bar;

    // gslpp::complex ADskpk0 =  -(TD + CD + kappa);
    gslpp::complex ADskpk0 =  -(TD + CD + kappa_p);
    gslpp::complex ADskpk0B = ckm_dcs_bar * ADskpk0;
    ADskpk0 = ckm_dcs * ADskpk0;

    if (obstype == 1) return PS(MDS,MKP,MKS) * BR_bar((ADskpk0 + ADskpk0bar).abs2(), (ADskpk0B + ADskpk0barB).abs2())/2. / GDS;
    else return 0.;
}

BRDskpkS_kpkL::BRDskpkS_kpkL(const StandardModel& SM_i)
: myObservables(SM_i), myBRDskpkS(SM_i, 1), myBRDskpkL(SM_i, 1)
{}

double BRDskpkS_kpkL::computeThValue()
{
    return myBRDskpkS.computeThValue() + myBRDskpkL.computeThValue();

}

BRDsppet::BRDsppet(const StandardModel& SM_i, int obstype_i)
: myObservables(SM_i)
{   obstype = obstype_i; }

double BRDsppet::computeThValue()
{
    updateParameters();

    gslpp::complex ADsppet8 = -sqrt(2./3.)/5. * (2. * Tp - 3. * Cp + D) * eIdelta_1s - sqrt(6.)/5. * (Tp + Cp + kappa);
    gslpp::complex ADsppetB = ckm_cf_bar * ADsppet8;
    gslpp::complex ADsppet = ckm_cf * ADsppet8;

    if (obstype == 1) return PS(MDS,MPIP,MET) * BR_bar((ADsppet).abs2(), (ADsppetB).abs2()) / GDS;
    else return 0.;
}

BRpmkp::BRpmkp(const StandardModel& SM_i, int obstype_i)
: myObservables(SM_i)
{   obstype = obstype_i; }

double BRpmkp::computeThValue()
{
    updateParameters();

    // gslpp::complex Apmkp = -1./5. * (3. * TD - 2. * CD + K) * eIdelta_12 - 2./5. * (TD + CD + kappa);
    gslpp::complex Apmkp = -1./5. * (3. * TD - 2. * CD + K) * eIdelta_12 - 2./5. * (TD + CD + kappa_p);
    gslpp::complex ApmkpB = ckm_dcs_bar * Apmkp;
    Apmkp = ckm_dcs * Apmkp;

    if (obstype == 1) return PS(MD0,MPIM,MKP) * BR_bar((Apmkp).abs2(), (ApmkpB).abs2()) / GD0;
    else return 0.;
}

BRp0kp::BRp0kp(const StandardModel& SM_i, int obstype_i)
: myObservables(SM_i)
{   obstype = obstype_i; }

double BRp0kp::computeThValue()
{
    updateParameters();

    // gslpp::complex Ap0kp = 1./5./sqrt(2.) * (2. * TD - 3. * CD + D - Kp) * eIdelta_12 + 3./5./sqrt(2.) * (TD + CD + kappa);
    gslpp::complex Ap0kp = 1./5./sqrt(2.) * (2. * TD - 3. * CD + D - Kp) * eIdelta_12 + 3./5./sqrt(2.) * (TD + CD + kappa_p);
    gslpp::complex Ap0kpB = ckm_dcs_bar * Ap0kp;
    Ap0kp = ckm_dcs * Ap0kp;

    if (obstype == 1) return PS(MDP,MPI0,MKP) * BR_bar((Ap0kp).abs2(), (Ap0kpB).abs2()) / GDP;
    else return 0.;
}

BRkpet::BRkpet(const StandardModel& SM_i, int obstype_i)
: myObservables(SM_i)
{   obstype = obstype_i; }

double BRkpet::computeThValue()
{
    updateParameters();

    // gslpp::complex Akpet8 = -1./5./sqrt(6.) * (2. * TD - 3. * CD + D - Kp) * eIdelta_12 - 3./5./sqrt(6.) * (TD + CD + kappa);
    gslpp::complex Akpet8 = -1./5./sqrt(6.) * (2. * TD - 3. * CD + D - Kp) * eIdelta_12 - 3./5./sqrt(6.) * (TD + CD + kappa_p);
    gslpp::complex AkpetB = ckm_dcs_bar * Akpet8;
    gslpp::complex Akpet = ckm_dcs * Akpet8;

    if (obstype == 1) return PS(MDP,MKP,MET) * BR_bar((Akpet).abs2(), (AkpetB).abs2()) / GDP;
    else return 0.;
}

BRppkm_pmkp::BRppkm_pmkp(const StandardModel& SM_i)
: myObservables(SM_i), myBRppkm(SM_i, 1), myBRpmkp(SM_i, 1)
{}

double BRppkm_pmkp::computeThValue()
{
    return myBRppkm.computeThValue() + myBRpmkp.computeThValue();

}

/* The strong phase difference for D0 to K pi. */

delta_kpi::delta_kpi(const StandardModel& SM_i)
: myObservables(SM_i)
{}

double delta_kpi::computeThValue()
{
    updateParameters();

    gslpp::complex Apmkp = -1./5. * (3. * TD - 2. * CD + K) * eIdelta_12 - 2./5. * (TD + CD + kappa_p);
    gslpp::complex Appkm = 1./5. * (3. * Tp - 2. * Cp - K) * eIdelta_12 + 2./5. * (Tp + Cp + kappa);

    return (-Apmkp/Appkm).arg()*180./M_PI;
}

/* Implementation of the Gronou Relations: arXiv:1501.03272*/

GR_R1::GR_R1(const StandardModel& SM_i)
: myObservables(SM_i)
{}

double GR_R1::computeThValue()
{
    updateParameters();

    gslpp::complex Apmkp = -1./5. * (3. * TD - 2. * CD + K) * eIdelta_12 - 2./5. * (TD + CD + kappa);
    gslpp::complex Appkm = 1./5. * (3. * Tp - 2. * Cp - K) * eIdelta_12 + 2./5. * (Tp + Cp + kappa);

    return Apmkp.abs()/Appkm.abs();

}

GR_R2::GR_R2(const StandardModel& SM_i)
: myObservables(SM_i)
{}

double GR_R2::computeThValue()
{
    updateParameters();

    gslpp::complex Akpkm = (T - 2./3. * C) * (3./20. * (eIdelta_0 + eIdelta_0_p) + (3./20. * cos2phi + 3./4./sqrt_10 * sin2phi) * (eIdelta_0_p - eIdelta_0) + 3./10. * eIdelta_1) + (T + C) * 2./5.;
    gslpp::complex Apppm = (T - 2./3. * C) * (-3./10. * (eIdelta_0 + eIdelta_0_p) + (-3./10. * cos2phi + 3./4./sqrt_10 * sin2phi) * (eIdelta_0_p - eIdelta_0)) - (T + C) * 2./5.;

    return Akpkm.abs()/Apppm.abs();

}

GR_R3::GR_R3(const StandardModel& SM_i)
: myObservables(SM_i)
{}

double GR_R3::computeThValue()
{
    updateParameters();

    gslpp::complex Apmkp = -1./5. * (3. * TD - 2. * CD + K) * eIdelta_12 - 2./5. * (TD + CD + kappa);
    gslpp::complex Appkm = 1./5. * (3. * Tp - 2. * Cp - K) * eIdelta_12 + 2./5. * (Tp + Cp + kappa);

    gslpp::complex Akpkm = (T - 2./3. * C) * (3./20. * (eIdelta_0 + eIdelta_0_p) + (3./20. * cos2phi + 3./4./sqrt_10 * sin2phi) * (eIdelta_0_p - eIdelta_0) + 3./10. * eIdelta_1) + (T + C) * 2./5.;
    gslpp::complex Apppm = (T - 2./3. * C) * (-3./10. * (eIdelta_0 + eIdelta_0_p) + (-3./10. * cos2phi + 3./4./sqrt_10 * sin2phi) * (eIdelta_0_p - eIdelta_0)) - (T + C) * 2./5.;

    return (Akpkm.abs() + Apppm.abs())/(Appkm.abs() + Apmkp.abs());

}

GR_R4::GR_R4(const StandardModel& SM_i)
: myObservables(SM_i)
{}

double GR_R4::computeThValue()
{
    updateParameters();

    gslpp::complex Apmkp = -1./5. * (3. * TD - 2. * CD + K) * eIdelta_12 - 2./5. * (TD + CD + kappa);
    gslpp::complex Appkm = 1./5. * (3. * Tp - 2. * Cp - K) * eIdelta_12 + 2./5. * (Tp + Cp + kappa);

    gslpp::complex Akpkm = (T - 2./3. * C) * (3./20. * (eIdelta_0 + eIdelta_0_p) + (3./20. * cos2phi + 3./4./sqrt_10 * sin2phi) * (eIdelta_0_p - eIdelta_0) + 3./10. * eIdelta_1) + (T + C) * 2./5.;
    gslpp::complex Apppm = (T - 2./3. * C) * (-3./10. * (eIdelta_0 + eIdelta_0_p) + (-3./10. * cos2phi + 3./4./sqrt_10 * sin2phi) * (eIdelta_0_p - eIdelta_0)) - (T + C) * 2./5.;

    return sqrt((Akpkm.abs() * Apppm.abs())/(Appkm.abs() * Apmkp.abs()));

}

GR_DR::GR_DR(const StandardModel& SM_i)
: myObservables(SM_i)
{}

double GR_DR::computeThValue()
{
    updateParameters();

    gslpp::complex Apmkp = -1./5. * (3. * TD - 2. * CD + K) * eIdelta_12 - 2./5. * (TD + CD + kappa);
    gslpp::complex Appkm = 1./5. * (3. * Tp - 2. * Cp - K) * eIdelta_12 + 2./5. * (Tp + Cp + kappa);

    gslpp::complex Akpkm = (T - 2./3. * C) * (3./20. * (eIdelta_0 + eIdelta_0_p) + (3./20. * cos2phi + 3./4./sqrt_10 * sin2phi) * (eIdelta_0_p - eIdelta_0) + 3./10. * eIdelta_1) + (T + C) * 2./5.;
    gslpp::complex Apppm = (T - 2./3. * C) * (-3./10. * (eIdelta_0 + eIdelta_0_p) + (-3./10. * cos2phi + 3./4./sqrt_10 * sin2phi) * (eIdelta_0_p - eIdelta_0)) - (T + C) * 2./5.;

    double R1 = Apmkp.abs()/Appkm.abs();
    double R2 = Akpkm.abs()/Apppm.abs();
    double R3 = (Akpkm.abs() + Apppm.abs())/(Appkm.abs() + Apmkp.abs());
    double R4 = sqrt((Akpkm.abs() * Apppm.abs())/(Appkm.abs() * Apmkp.abs()));

    return R3 - R4 + 1./8. * ((sqrt(2.*R1 - 1.) - 1.)*(sqrt(2.*R1 - 1.) - 1.) - (sqrt(2.*R2 - 1.) - 1.)*(sqrt(2.*R2 - 1.) - 1.));

}

GR_ReEps1::GR_ReEps1(const StandardModel& SM_i)
: myObservables(SM_i)
{}

double GR_ReEps1::computeThValue()
{
    updateParameters();

    gslpp::complex Apmkp = -1./5. * (3. * TD - 2. * CD + K) * eIdelta_12 - 2./5. * (TD + CD + kappa);
    gslpp::complex Appkm = 1./5. * (3. * Tp - 2. * Cp - K) * eIdelta_12 + 2./5. * (Tp + Cp + kappa);

    double R1 = Apmkp.abs()/Appkm.abs();

    return 1./2. * (sqrt(2.*R1 - 1.) - 1.);

}

GR_ReEps2::GR_ReEps2(const StandardModel& SM_i)
: myObservables(SM_i)
{}

double GR_ReEps2::computeThValue()
{
    updateParameters();

    gslpp::complex Akpkm = (T - 2./3. * C) * (3./20. * (eIdelta_0 + eIdelta_0_p) + (3./20. * cos2phi + 3./4./sqrt_10 * sin2phi) * (eIdelta_0_p - eIdelta_0) + 3./10. * eIdelta_1) + (T + C) * 2./5.;
    gslpp::complex Apppm = (T - 2./3. * C) * (-3./10. * (eIdelta_0 + eIdelta_0_p) + (-3./10. * cos2phi + 3./4./sqrt_10 * sin2phi) * (eIdelta_0_p - eIdelta_0)) - (T + C) * 2./5.;

    double R2 = Akpkm.abs()/Apppm.abs();

    return 1./2. * (sqrt(2.*R2 - 1.) - 1.);

}

GR_ReEps0p::GR_ReEps0p(const StandardModel& SM_i)
: myObservables(SM_i), myBRp0kS(SM_i, 1), myBRp0kL(SM_i, 1)
{}

double GR_ReEps0p::computeThValue()
{
    updateParameters();

    double tanThetaC2 = (-V_us*V_cd/(V_ud*V_cs)).abs();

    return (myBRp0kS.computeThValue() - myBRp0kL.computeThValue())/(myBRp0kS.computeThValue() + myBRp0kL.computeThValue())/(4.*tanThetaC2) - 1./2.;

}

GR_RD0::GR_RD0(const StandardModel& SM_i)
: myObservables(SM_i), myBRp0kS(SM_i, 1), myBRp0kL(SM_i, 1)
{}

double GR_RD0::computeThValue()
{
    updateParameters();

    return (myBRp0kS.computeThValue() - myBRp0kL.computeThValue())/(myBRp0kS.computeThValue() + myBRp0kL.computeThValue());

}

GR_RDp::GR_RDp(const StandardModel& SM_i)
: myObservables(SM_i), myBRppkS(SM_i, 1), myBRppkL(SM_i, 1)
{}

double GR_RDp::computeThValue()
{
    updateParameters();

    return (myBRppkS.computeThValue() - myBRppkL.computeThValue())/(myBRppkS.computeThValue() + myBRppkL.computeThValue());

}

GR_RDsp::GR_RDsp(const StandardModel& SM_i)
: myObservables(SM_i), myBRDskpkS(SM_i, 1), myBRDskpkL(SM_i, 1)
{}

double GR_RDsp::computeThValue()
{
    updateParameters();

    return (myBRDskpkS.computeThValue() - myBRDskpkL.computeThValue())/(myBRDskpkS.computeThValue() + myBRDskpkL.computeThValue());

}

DACP::DACP(const StandardModel& SM_i)
: myObservables(SM_i), myBRpppm(SM_i, 2), myBRkpkm(SM_i, 2)
{}

double DACP::computeThValue()
{
    updateParameters();

    return (myBRkpkm.computeThValue() - myBRpppm.computeThValue());

}

delta_1_s::delta_1_s(const StandardModel& SM_i)
: myObservables(SM_i)
{}

double delta_1_s::computeThValue()
{
    updateParameters();

    return delta_1s;

}

delta_12_s::delta_12_s(const StandardModel& SM_i)
: myObservables(SM_i)
{}

double delta_12_s::computeThValue()
{
    updateParameters();

    return delta_12s;

}

A1_T::A1_T(const StandardModel& SM_i)
: myObservables(SM_i)
{}

double A1_T::computeThValue()
{
    updateParameters();

    return Kp/5./T;

}

A2_T::A2_T(const StandardModel& SM_i)
: myObservables(SM_i)
{}

double A2_T::computeThValue()
{
    updateParameters();

    return K/5./T;

}

/* The Reduced Matrix Elements in the SU(3) formalism */
/* CF Sector */

RD_8_6_1::RD_8_6_1(const StandardModel& SM_i, int i)
: myObservables(SM_i)
{ type = i; }

double RD_8_6_1::computeThValue()
{
    updateParameters();

    if (type == 1) return (((2.*Tp - 3.*Cp + D)*eIdelta_1s + (3.*Tp - 2.*Cp - K)*eIdelta_12)/sqrt_30).abs();
    if (type == 2) return (((2.*Tp - 3.*Cp + D)*eIdelta_1s + (3.*Tp - 2.*Cp - K)*eIdelta_12)/sqrt_30).arg();
    else return 0.;

}

RD_8_15_1_CF::RD_8_15_1_CF(const StandardModel& SM_i, int i)
: myObservables(SM_i)
{ type = i; }

double RD_8_15_1_CF::computeThValue()
{
    updateParameters();

    if (type == 1) return (((2.*Tp - 3.*Cp + D)*eIdelta_1s - (3.*Tp - 2.*Cp - K)*eIdelta_12)/sqrt_30).abs();
    if (type == 2) return (((2.*Tp - 3.*Cp + D)*eIdelta_1s - (3.*Tp - 2.*Cp - K)*eIdelta_12)/sqrt_30).arg();
    else return 0.;

}

RD_27_15_1_CF::RD_27_15_1_CF(const StandardModel& SM_i)
: myObservables(SM_i)
{}

double RD_27_15_1_CF::computeThValue()
{
    updateParameters();

    return 3./sqrt(5.) * (Tp + Cp + kappa);

}

/* DCS Sector */

RD_8_6_0::RD_8_6_0(const StandardModel& SM_i, int i)
: myObservables(SM_i)
{ type = i; }

double RD_8_6_0::computeThValue()
{
    updateParameters();

    if (type == 1) return (-(5.*TD - 5.*CD + D + K - Kp)*eIdelta_12/sqrt_30).abs();
    if (type == 2) return (-(5.*TD - 5.*CD + D + K - Kp)*eIdelta_12/sqrt_30).arg();
    else return 0.;

}

RD_8_15_1_DCS::RD_8_15_1_DCS(const StandardModel& SM_i, int i)
: myObservables(SM_i)
{ type = i; }

double RD_8_15_1_DCS::computeThValue()
{
    updateParameters();

    if (type == 1) return (-(TD + CD - D + K + Kp)*eIdelta_12/sqrt_30).abs();
    if (type == 2) return (-(TD + CD - D + K + Kp)*eIdelta_12/sqrt_30).arg();
    else return 0.;

}

RD_27_15_1_DCS::RD_27_15_1_DCS(const StandardModel& SM_i)
: myObservables(SM_i)
{}

double RD_27_15_1_DCS::computeThValue()
{
    updateParameters();

    return 3./sqrt(5.) * (TD + CD + kappa);

}

/* SCS Sector */

RD_1_3_12::RD_1_3_12(const StandardModel& SM_i, int i)
: myObservables(SM_i)
{ type = i; }

double RD_1_3_12::computeThValue()
{
    updateParameters();

    if (type == 1) return (-3./2./sqrt_10 * (T - 2./3.*C)*(eIdelta_0 - eIdelta_0_p)*sin2phi).abs();
    if (type == 2) return (-3./2./sqrt_10 * (T - 2./3.*C)*(eIdelta_0 - eIdelta_0_p)*sin2phi).arg();
    else return 0.;

}

RD_8_3_12::RD_8_3_12(const StandardModel& SM_i, int i)
: myObservables(SM_i)
{ type = i; }

double RD_8_3_12::computeThValue()
{
    updateParameters();

    if (type == 1) return (-3./8./sqrt_10 * (T - 2./3.*C) * ((eIdelta_0 + eIdelta_0_p) - cos2phi * (eIdelta_0 - eIdelta_0_p)) + 1./4./sqrt_10 * (7.*T - 8.*C + 2.*D)*eIdelta_1 - 1./2./sqrt_10 * (2.*T - 3.*C + D - Kp)*eIdelta_12s).abs();
    if (type == 2) return (-3./8./sqrt_10 * (T - 2./3.*C) * ((eIdelta_0 + eIdelta_0_p) - cos2phi * (eIdelta_0 - eIdelta_0_p)) + 1./4./sqrt_10 * (7.*T - 8.*C + 2.*D)*eIdelta_1 - 1./2./sqrt_10 * (2.*T - 3.*C + D - Kp)*eIdelta_12s).arg();
    else return 0.;

}

RD_8_6_12::RD_8_6_12(const StandardModel& SM_i, int i)
: myObservables(SM_i)
{ type = i; }

double RD_8_6_12::computeThValue()
{
    updateParameters();

    if (type == 1) return (-9./8./sqrt(15.) * (T - 2./3.*C) * ((eIdelta_0 + eIdelta_0_p) - cos2phi * (eIdelta_0 - eIdelta_0_p)) - 1./4./sqrt(15.) * (7.*T - 8.*C + 2.*D)*eIdelta_1 - 1./2./sqrt(15.) * (2.*T - 3.*C + D - Kp)*eIdelta_12s).abs();
    if (type == 2) return (-9./8./sqrt(15.) * (T - 2./3.*C) * ((eIdelta_0 + eIdelta_0_p) - cos2phi * (eIdelta_0 - eIdelta_0_p)) - 1./4./sqrt(15.) * (7.*T - 8.*C + 2.*D)*eIdelta_1 - 1./2./sqrt(15.) * (2.*T - 3.*C + D - Kp)*eIdelta_12s).arg();
    else return 0.;

}

RD_8_15_12::RD_8_15_12(const StandardModel& SM_i, int i)
: myObservables(SM_i)
{ type = i; }

double RD_8_15_12::computeThValue()
{
    updateParameters();

    if (type == 1) return (9./8./sqrt_10 * (T - 2./3.*C) * ((eIdelta_0 + eIdelta_0_p) - cos2phi * (eIdelta_0 - eIdelta_0_p)) - 1./12./sqrt_10 * (7.*T - 8.*C + 2.*D)*eIdelta_1 - 1./2./sqrt_10 * (2.*T - 3.*C + D - Kp)*eIdelta_12s).abs();
    if (type == 2) return (9./8./sqrt_10 * (T - 2./3.*C) * ((eIdelta_0 + eIdelta_0_p) - cos2phi * (eIdelta_0 - eIdelta_0_p)) - 1./12./sqrt_10 * (7.*T - 8.*C + 2.*D)*eIdelta_1 - 1./2./sqrt_10 * (2.*T - 3.*C + D - Kp)*eIdelta_12s).arg();
    else return 0.;

}

RD_8_15_32::RD_8_15_32(const StandardModel& SM_i, int i)
: myObservables(SM_i)
{ type = i; }

double RD_8_15_32::computeThValue()
{
    updateParameters();

    if (type == 1) return ((-1./3./sqrt(5.) * (T + C - D)*eIdelta_1)).abs();
    if (type == 2) return ((-1./3./sqrt(5.) * (T + C - D)*eIdelta_1)).arg();
    else return 0.;

}

RD_27_15_12::RD_27_15_12(const StandardModel& SM_i)
: myObservables(SM_i)
{}

double RD_27_15_12::computeThValue()
{
    updateParameters();

    return -2.*sqrt(3./5.) * (T + C);

}

RD_27_15_32::RD_27_15_32(const StandardModel& SM_i)
: myObservables(SM_i)
{}

double RD_27_15_32::computeThValue()
{
    updateParameters();

    return sqrt(6./5.) * (T + C);

}
