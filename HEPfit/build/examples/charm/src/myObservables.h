/*
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MYOBSERVABLES_H
#define	MYOBSERVABLES_H

#pragma warning disable 2196
#include <HEPfit.h>
#include "myModel.h"


/**
 * @class myObservables
 * @brief A class for the gg -> 4l.
 */
class myObservables : public ThObservable {
public:
    myObservables(const StandardModel& SM_i);
    virtual ~myObservables();
    void updateParameters();

protected:

    /* Define the couplings here. */
    double T;
    double C;
    // double eps;
    double Tp;
    double Cp;
    // double eps_p;
    double TD;
    double CD;
    double kappa;
    double kappa_p;
    double D;
    double K;
    double Kp;
    double delta_0;
    double delta_0_p;
    double delta_12;
    double eps_d;
    double delta_12s;
    double delta_1;
    double delta_1s;
    double phi;

    double cos2phi, sin2phi;
    gslpp::complex eIdelta_0, eIdelta_0_p, eIdelta_1, eIdelta_12, eIdelta_1s, eIdelta_12s;

    /* CPV Parameters */

    double Delta_3;
    double Delta_4;
    double P;

    double GF;

    // CKM Parameters
    gslpp::complex V_ud, V_us;
    gslpp::complex V_cd, V_cs;
    gslpp::complex V_ub, V_cb;
    gslpp::complex ckm_scs, ckm_scs_bar;
    gslpp::complex ckm_b, ckm_b_bar;
    gslpp::complex ckm_cf, ckm_cf_bar;
    gslpp::complex ckm_dcs, ckm_dcs_bar;
    gslpp::complex ckm_CPV, ckm_CPV_bar;

    double PS(const double M, const double m1, const double m2) const;
    double BR_bar(double amp, double amp_bar);
    double aCP(double amp, double amp_bar);

    double MD0, MDP, MDS;
    double MPIP, MPIM, MKP, MKM;
    double MPI0, MKS, MKL;
    double MET, METP;
    double GD0, GDP, GDS;
    gslpp::complex ii;
    gslpp::matrix<gslpp::complex> VCKM;

    double sqrt_10;
    double sqrt_30;
    double sqrt_3_5;

private:
    const myModel * my_model;

};

class PENG_T : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    PENG_T(const StandardModel& SM_i);

    /**
     * @return BRpppm
     */
    double computeThValue ();

private:
};

class BRpppm : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BRpppm(const StandardModel& SM_i, int obstype_i);

    /**
     * @return BRpppm
     */
    double computeThValue ();

private:
    int obstype;

};

class BRp0p0 : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BRp0p0(const StandardModel& SM_i, int obstype_i);

    /**
     * @return BRp0p0
     */
    double computeThValue ();

private:
    int obstype;

};

class BRkpkm : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BRkpkm(const StandardModel& SM_i, int obstype_i);

    /**
     * @return BRkpkm
     */
    double computeThValue ();

private:
    int obstype;

};

class BRkSkS : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BRkSkS(const StandardModel& SM_i, int obstype_i);

    /**
     * @return BRksks
     */
    double computeThValue ();

private:
    int obstype;

};

class BRetet : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BRetet(const StandardModel& SM_i, int obstype_i);

    /**
     * @return BRetet
     */
    double computeThValue ();

private:
    int obstype;

};

class BRp0et : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BRp0et(const StandardModel& SM_i, int obstype_i);

    /**
     * @return BRp0et
     */
    double computeThValue ();

private:
    int obstype;

};

class BRkpkS : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BRkpkS(const StandardModel& SM_i, int obstype_i);

    /**
     * @return BRkpkS
     */
    double computeThValue ();

private:
    int obstype;

};

class BRppp0 : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BRppp0(const StandardModel& SM_i, int obstype_i);

    /**
     * @return BRppp0
     */
    double computeThValue ();

private:
    int obstype;

};

class BRppet : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BRppet(const StandardModel& SM_i, int obstype_i);

    /**
     * @return BRppet
     */
    double computeThValue ();

private:
    int obstype;

};

class BRDsppkS : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BRDsppkS(const StandardModel& SM_i, int obstype_i);

    /**
     * @return BRDsppkS
     */
    double computeThValue ();

private:
    int obstype;

};

class BRDsp0kp : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BRDsp0kp(const StandardModel& SM_i, int obstype_i);

    /**
     * @return BRDsp0kp
     */
    double computeThValue ();

private:
    int obstype;

};


class BRDskpet : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BRDskpet(const StandardModel& SM_i, int obstype_i);

    /**
     * @return BRDskpet
     */
    double computeThValue ();

private:
    int obstype;

};

/*                                                    */
/* amplitudes and CP-conjugate amplitudes A (CF Modes)*/
/*                                                    */

class BRppkm : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BRppkm(const StandardModel& SM_i, int obstype_i);

    /**
     * @return BRppkm
     */
    double computeThValue ();

private:
    int obstype;

};

class BRp0kS : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BRp0kS(const StandardModel& SM_i, int obstype_i);

    /**
     * @return BRp0kS
     */
    double computeThValue ();

private:
    int obstype;

};

class BRp0kL : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BRp0kL(const StandardModel& SM_i, int obstype_i);

    /**
     * @return BRp0kL
     */
    double computeThValue ();

private:
    int obstype;

};

class BRkSet : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BRkSet(const StandardModel& SM_i, int obstype_i);

    /**
     * @return BRkSet
     */
    double computeThValue ();

private:
    int obstype;

};

class BRppkS : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BRppkS(const StandardModel& SM_i, int obstype_i);

    /**
     * @return BRppkS
     */
    double computeThValue ();

private:
    int obstype;

};

class BRppkL : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BRppkL(const StandardModel& SM_i, int obstype_i);

    /**
     * @return BRppkL
     */
    double computeThValue ();

private:
    int obstype;

};

class BRDskpkS : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BRDskpkS(const StandardModel& SM_i, int obstype_i);

    /**
     * @return BRDskpkS
     */
    double computeThValue ();

private:
    int obstype;

};

class BRDskpkL : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BRDskpkL(const StandardModel& SM_i, int obstype_i);

    /**
     * @return BRDskpkS
     */
    double computeThValue ();

private:
    int obstype;

};

class BRDskpkS_kpkL : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BRDskpkS_kpkL(const StandardModel& SM_i);

    /**
     * @return BRDskpkS
     */
    double computeThValue ();

private:
    BRDskpkS myBRDskpkS;
    BRDskpkL myBRDskpkL;
};

class BRDsppet : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BRDsppet(const StandardModel& SM_i, int obstype_i);

    /**
     * @return BRDsppet
     */
    double computeThValue ();

private:
    int obstype;

};


/*                                                    */
/* amplitudes and CP-conjugate amplitudes A (DCS Modes)*/
/*                                                    */

class BRpmkp : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BRpmkp(const StandardModel& SM_i, int obstype_i);

    /**
     * @return BRpmkp
     */
    double computeThValue ();

private:
    int obstype;

};

class BRp0kp : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BRp0kp(const StandardModel& SM_i, int obstype_i);

    /**
     * @return BRp0kp
     */
    double computeThValue ();

private:
    int obstype;

};

class BRkpet : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BRkpet(const StandardModel& SM_i, int obstype_i);

    /**
     * @return BRkpet
     */
    double computeThValue ();

private:
    int obstype;

};

class BRppkm_pmkp : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    BRppkm_pmkp(const StandardModel& SM_i);

    /**
     * @return BRDskpkS
     */
    double computeThValue ();

private:
    BRppkm myBRppkm;
    BRpmkp myBRpmkp;
};

class delta_kpi : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    delta_kpi(const StandardModel& SM_i);

    /**
     * @return delta_kpi
     */
    double computeThValue ();

private:

};

class GR_R1 : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    GR_R1(const StandardModel& SM_i);

    /**
     * @return GR_R1
     */
    double computeThValue ();

private:

};

class GR_R2 : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    GR_R2(const StandardModel& SM_i);

    /**
     * @return GR_R2
     */
    double computeThValue ();

private:

};

class GR_R3 : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    GR_R3(const StandardModel& SM_i);

    /**
     * @return GR_R3
     */
    double computeThValue ();

private:

};

class GR_R4 : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    GR_R4(const StandardModel& SM_i);

    /**
     * @return GR_R4
     */
    double computeThValue ();

private:

};

class GR_DR : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    GR_DR(const StandardModel& SM_i);

    /**
     * @return GR_DR
     */
    double computeThValue ();

private:

};

class GR_ReEps1 : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    GR_ReEps1(const StandardModel& SM_i);

    /**
     * @return GR_ReEps1
     */
    double computeThValue ();

private:

};

class GR_ReEps2 : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    GR_ReEps2(const StandardModel& SM_i);

    /**
     * @return GR_ReEps2
     */
    double computeThValue ();

private:

};

class GR_ReEps0p : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    GR_ReEps0p(const StandardModel& SM_i);

    /**
     * @return GR_ReEps0p
     */
    double computeThValue ();

private:

    BRp0kS myBRp0kS;
    BRp0kL myBRp0kL;

};

class GR_RD0 : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    GR_RD0(const StandardModel& SM_i);

    /**
     * @return GR_RD0
     */
    double computeThValue ();

private:

    BRp0kS myBRp0kS;
    BRp0kL myBRp0kL;

};

class GR_RDp : public myObservables {
    public:

    /**
     * @brief Constructor.
     */
    GR_RDp(const StandardModel& SM_i);

    /**
     * @return GR_RDp
     */
    double computeThValue ();

    private:

    BRppkS myBRppkS;
    BRppkL myBRppkL;

};

class GR_RDsp : public myObservables {
    public:

    /**
     * @brief Constructor.
     */
    GR_RDsp(const StandardModel& SM_i);

    /**
     * @return GR_RDsp
     */
    double computeThValue ();

    private:

    BRDskpkS myBRDskpkS;
    BRDskpkL myBRDskpkL;

};

class DACP : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    DACP(const StandardModel& SM_i);

    /**
     * @return DACP
     */
    double computeThValue ();

private:

    BRpppm myBRpppm;
    BRkpkm myBRkpkm;

};

class delta_1_s : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    delta_1_s(const StandardModel& SM_i);

    /**
     * @return delta_1_s
     */
    double computeThValue ();

private:

};

class delta_12_s : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    delta_12_s(const StandardModel& SM_i);

    /**
     * @return delta_12_s
     */
    double computeThValue ();

private:

};

class A1_T : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    A1_T(const StandardModel& SM_i);

    /**
     * @return A1_T
     */
    double computeThValue ();

private:

};

class A2_T : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    A2_T(const StandardModel& SM_i);

    /**
     * @return A2_T
     */
    double computeThValue ();

private:

};

/* The Reduced Matrix Elements in the SU(3) formalism */
/* CF Sector */

class RD_8_6_1 : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    RD_8_6_1(const StandardModel& SM_i, int i);

    /**
     * @return RD_8_6_1
     */
    double computeThValue ();

private:

    int type;

};

class RD_8_15_1_CF : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    RD_8_15_1_CF(const StandardModel& SM_i, int i);

    /**
     * @return RD_8_15_1_CF
     */
    double computeThValue ();

private:

    int type;
};

class RD_27_15_1_CF : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    RD_27_15_1_CF(const StandardModel& SM_i);

    /**
     * @return RD_27_15_1_CF
     */
    double computeThValue ();

private:

};

/* DCS Sector */

class RD_8_6_0 : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    RD_8_6_0(const StandardModel& SM_i, int i);

    /**
     * @return RD_8_6_0
     */
    double computeThValue ();

private:

    int type;

};

class RD_8_15_1_DCS : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    RD_8_15_1_DCS(const StandardModel& SM_i, int i);

    /**
     * @return RD_8_15_1_DCS
     */
    double computeThValue ();

private:

    int type;
};

class RD_27_15_1_DCS : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    RD_27_15_1_DCS(const StandardModel& SM_i);

    /**
     * @return RD_27_15_1_DCS
     */
    double computeThValue ();

private:

};

/* SCS Sector */

class RD_1_3_12 : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    RD_1_3_12(const StandardModel& SM_i, int i);

    /**
     * @return RD_1_3_12
     */
    double computeThValue ();

private:

    int type;

};

class RD_8_3_12 : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    RD_8_3_12(const StandardModel& SM_i, int i);

    /**
     * @return RD_8_3_12
     */
    double computeThValue ();

private:

    int type;
};

class RD_8_6_12 : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    RD_8_6_12(const StandardModel& SM_i, int i);

    /**
     * @return RD_8_6_12
     */
    double computeThValue ();

private:

    int type;

};

class RD_8_15_12 : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    RD_8_15_12(const StandardModel& SM_i, int i);

    /**
     * @return RD_8_15_12
     */
    double computeThValue ();

private:

    int type;
};

class RD_8_15_32 : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    RD_8_15_32(const StandardModel& SM_i, int i);

    /**
     * @return RD_8_15_32
     */
    double computeThValue ();

private:

    int type;
};

class RD_27_15_12 : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    RD_27_15_12(const StandardModel& SM_i);

    /**
     * @return RD_27_15_12
     */
    double computeThValue ();

private:

};

class RD_27_15_32 : public myObservables {
public:

    /**
     * @brief Constructor.
     */
    RD_27_15_32(const StandardModel& SM_i);

    /**
     * @return RD_27_15_32
     */
    double computeThValue ();

private:

};


#endif	/* MYOBSERVABLES_H */
