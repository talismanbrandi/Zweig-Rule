/*
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "myModel.h"
#include <boost/ref.hpp>

/* Define mandatory model parameters here. */
const std::string myModel::myModelvars[NmyModelvars] = {"T", "C", "kappa", "kappa_p", "D", "K", "Kp", "delta_0", "delta_0_p", "delta_12", "eps_d", "delta_1", "phi", "P_T", "Delta_3_T", "Delta_4_T"};

myModel::myModel()
:   StandardModel()
{
    /* Define all the parameters here and port them as observables too */
    ModelParamMap.insert(std::make_pair("T", std::cref(T)));
    ModelParamMap.insert(std::make_pair("C", std::cref(C)));
    ModelParamMap.insert(std::make_pair("kappa", std::cref(kappa)));
    ModelParamMap.insert(std::make_pair("kappa_p", std::cref(kappa_p)));
    ModelParamMap.insert(std::make_pair("D", std::cref(D)));
    ModelParamMap.insert(std::make_pair("K", std::cref(K)));
    ModelParamMap.insert(std::make_pair("Kp", std::cref(Kp)));
    ModelParamMap.insert(std::make_pair("delta_0", std::cref(delta_0)));
    ModelParamMap.insert(std::make_pair("delta_0_p", std::cref(delta_0_p)));
    ModelParamMap.insert(std::make_pair("delta_12", std::cref(delta_12)));
    ModelParamMap.insert(std::make_pair("eps_d", std::cref(eps_d)));
    ModelParamMap.insert(std::make_pair("delta_1", std::cref(delta_1)));
    ModelParamMap.insert(std::make_pair("phi", std::cref(phi)));
    ModelParamMap.insert(std::make_pair("P_T", std::cref(P_T)));
    ModelParamMap.insert(std::make_pair("Delta_3_T", std::cref(Delta_3_T)));
    ModelParamMap.insert(std::make_pair("Delta_4_T", std::cref(Delta_4_T)));
}

myModel::~myModel()
{
    if (IsModelInitialized()) {
        /* Destroy whatever you want, e.g. potentially dangerous pointers. */
    }
}

/* Initialize model here */
bool myModel::InitializeModel()
{
    //condition = false;
    setModelInitialized(StandardModel::InitializeModel());
    return(true);
}

bool myModel::Init(const std::map<std::string, double>& DPars)
{
    return(StandardModel::Init(DPars));
}

/* Do whatever is necessary before parameters are updated by the MCMC. */
bool myModel::PreUpdate()
{
    if(!StandardModel::PreUpdate()) return (false);
    return (true);
}

/* Model update method used be the MCMC to update the model parameters. */
bool myModel::Update(const std::map<std::string, double>& DPars)
{
    if(!PreUpdate()) return (false);

    UpdateError = false;

    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);

    if (UpdateError) return (false);

    if(!PostUpdate()) return (false);

    return (true);

    return (true);
}

/* Postupdate method to update whatever is needed after the model parameters are updated */
bool myModel::PostUpdate()
{
    if(!StandardModel::PostUpdate()) return (false);

    return (true);
}

/* Model parameters and their derived quantities can be set here. */
void myModel::setParameter(const std::string name, const double& value)
{
    if(name.compare("T") == 0)
        T = value;
    else if(name.compare("C") == 0)
        C = value;
    else if(name.compare("kappa") == 0)
        kappa = value;
    else if(name.compare("kappa_p") == 0)
        kappa_p = value;
    else if(name.compare("D") == 0)
        D = value;
    else if(name.compare("K") == 0)
        K = value;
    else if(name.compare("Kp") == 0)
        Kp = value;
    else if(name.compare("delta_0") == 0)
        delta_0 = value;
    else if(name.compare("delta_0_p") == 0)
        delta_0_p = value;
    else if(name.compare("delta_12") == 0)
        delta_12 = value;
    else if(name.compare("eps_d") == 0)
        eps_d = value;
    else if(name.compare("delta_1") == 0)
        delta_1 = value;
    else if(name.compare("phi") == 0)
        phi = value;
    else if(name.compare("P_T") == 0)
        P_T = value;
    else if(name.compare("Delta_3_T") == 0)
        Delta_3_T = value;
    else if(name.compare("Delta_4_T") == 0)
        Delta_4_T = value;
    else
        StandardModel::setParameter(name,value);
}

bool myModel::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NmyModelvars; i++) {
        if (DPars.find(myModelvars[i]) == DPars.end()) {
            std::cout << "missing mandatory myModel parameter " << myModelvars[i] << std::endl;
            return false;
        }
    }
    return(StandardModel::CheckParameters(DPars));
}
