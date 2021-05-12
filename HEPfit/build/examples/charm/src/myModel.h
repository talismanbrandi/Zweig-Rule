/*
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MYMODEL_H
#define	MYMODEL_H

#pragma warning disable 2196
#include <HEPfit.h>

/**
 * @class myModel
 * @brief My own Model.
 */
class myModel: public StandardModel {
public:

    static const int NmyModelvars = 16; /* Define number of mandatory parameters in the model. */
    static const std::string myModelvars[NmyModelvars]; /* Vector of model variable names. */

    /**
     * @brief myModel constructor
     */
    myModel();

    /**
     * @brief myModel destructor
     */
    ~myModel();

    virtual bool InitializeModel();

    virtual bool Init(const std::map<std::string, double>& DPars);

    virtual bool PreUpdate();

    virtual bool Update(const std::map<std::string, double>& DPars);

    virtual bool PostUpdate();

    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    //virtual bool setFlag(const std::string name, const bool value);


    /**
     *
     * @return the parameter T
     */
    double getT() const
    {
        return T;
    };

    /**
     *
     * @return the parameter C
     */
    double getC() const
    {
        return C;
    };

    /**
     *
     * @return the parameter kappa
     */
    double getkappa() const
    {
        return kappa;
    };

    /**
     *
     * @return the parameter kappa_p
     */
    double getkappa_p() const
    {
        return kappa_p;
    };

    /**
     *
     * @return the parameter D
     */
    double getD() const
    {
        return D;
    };

    /**
     *
     * @return the parameter K
     */
    double getK() const
    {
        return K;
    };

    /**
     *
     * @return the parameter Kp
     */
    double getKp() const
    {
        return Kp;
    };

    /**
     *
     * @return the parameter delta_0
     */
    double getdelta_0() const
    {
        return delta_0;
    };

    /**
     *
     * @return the parameter delta_0_p
     */
    double getdelta_0_p() const
    {
        return delta_0_p;
    };

    /**
     *
     * @return the parameter delta_12
     */
    double getdelta_12() const
    {
        return delta_12;
    };

    /**
     *
     * @return the parameter eps_d
     */
    double geteps_d() const
    {
        return eps_d;
    };

    /**
     *
     * @return the parameter delta_1
     */
    double getdelta_1() const
    {
        return delta_1;
    };

    /**
     *
     * @return the parameter phi
     */
    double getphi() const
    {
        return phi;
    };

    /**
     *
     * @return the cparameter P_T
     */
    double getP_T() const
    {
        return P_T;
    };

    /**
     *
     * @return the parameter Delta_3_T
     */
    double getDelta_3_T() const
    {
        return Delta_3_T;
    };

    /**
     *
     * @return the parameter Delta_4_T
     */
    double getDelta_4_T() const
    {
        return Delta_4_T;
    };

protected:

    virtual void setParameter(const std::string, const double&);

private:

    double T, C, kappa, kappa_p, D, K, Kp, delta_0, delta_0_p, delta_12, eps_d, delta_1, phi; /* Model Parameters */
    double P_T, Delta_3_T, Delta_4_T; /* Model Parameters CPV*/

};

#endif	/* MYMODEL_H */
