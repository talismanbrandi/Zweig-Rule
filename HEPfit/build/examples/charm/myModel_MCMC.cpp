/*
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

/**
 * @example myModel_MCMC.cpp
 * @example myModel.h
 * @example myModel.cpp
 * @example myObservables.h
 * @example myObservables.cpp
 * This is an example of how to add user-defined model and observables
 * and to perform a Bayesian analysis with the Markov Chain Monte Carlo.
 *
 */

#pragma warning disable 2196
#include <iostream>
#include <HEPfit.h>
#include <boost/bind.hpp>
#include <TSystem.h>
#include "myModel.h"
#include "myObservables.h"

/* Necessary if MPI support is enabled during compilation. */
#ifdef _MPI
#include <mpi.h>
#endif

int main(int argc, char** argv)
{

/* Necessary if MPI support is enabled during compilation. */
#ifdef _MPI
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    /* In our MPI implementation the process with rank 0 is the master. */
    int rank = 0;
#endif

    try {

        if(argc < 3){
            /* Print usage and exit. */
            if (rank == 0) std::cout << "\nusage: " << argv[0] << " ModelConf.conf MonteCarlo.conf[--noMC]\n" << std::endl;
            return EXIT_SUCCESS;
        }

        /* Define the model configuration file.                                */
        /* Here it is passed as the first argument to the executable. The      */
        /* model configuration file provides the values with errors for the    */
        /* mandatory model parameters, as well as the list of observables,     */
        /* observables2D, correlated Gaussian observables.                     */
        /* See documentation for details.                                      */
        std::string ModelConf = argv[1];

        /* Define the Monte Carlo configuration file.                         */
        /* Here it is passed as the second argument to the executable. The    */
        /* Monte Carlo configuration file provides the parameters used in the */
        /* Monte Carlo run. See documentation for details.                    */
        std::string MCMCConf = argv[2];

        /* Define the ROOT output file (w/o extension, empty string will set it to MCout) */
        std::string FileOut = "";

        /* Define the optional job tag. */
        // std::string directory = "";
        std::string JobTag = "";
        // if(argc == 4) {
        //     directory = argv[3];
        //     JobTag = ("_" + directory).c_str();
        // }

        /* Create objects of the classes ModelFactory and ThObsFactory */
        ThObsFactory ThObsF;
        ModelFactory ModelF;
        myModel my_model;

        /* register user-defined model named ModelName defined in class ModelClass using the following syntax: */
        ModelF.addModelToFactory("myModel", boost::factory<myModel*>() );

        /* register user-defined ThObservable named ThObsName defined in class ThObsClass using the following syntax: */

        ThObsF.addObsToFactory("PENG_T", boost::factory<PENG_T*>() );

        ThObsF.addObsToFactory("BRpppm", boost::bind(boost::factory<BRpppm*>(), _1, 1) );
        ThObsF.addObsToFactory("BRp0p0", boost::bind(boost::factory<BRp0p0*>(), _1, 1) );
        ThObsF.addObsToFactory("BRkpkm", boost::bind(boost::factory<BRkpkm*>(), _1, 1) );
        ThObsF.addObsToFactory("BRkSkS", boost::bind(boost::factory<BRkSkS*>(), _1, 1) );
        ThObsF.addObsToFactory("BRetet", boost::bind(boost::factory<BRetet*>(), _1, 1) );
        ThObsF.addObsToFactory("BRp0et", boost::bind(boost::factory<BRp0et*>(), _1, 1) );
        ThObsF.addObsToFactory("BRkpkS", boost::bind(boost::factory<BRkpkS*>(), _1, 1) );
        ThObsF.addObsToFactory("BRkpkL", boost::bind(boost::factory<BRkpkS*>(), _1, 1) ); /* BRkpkS = BRkpkL by definition */
        ThObsF.addObsToFactory("BRppp0", boost::bind(boost::factory<BRppp0*>(), _1, 1) );
        ThObsF.addObsToFactory("BRppet", boost::bind(boost::factory<BRppet*>(), _1, 1) );
        ThObsF.addObsToFactory("BRDsppkS", boost::bind(boost::factory<BRDsppkS*>(), _1, 1) );
        ThObsF.addObsToFactory("BRDsp0kp", boost::bind(boost::factory<BRDsp0kp*>(), _1, 1) );
        ThObsF.addObsToFactory("BRkpet", boost::bind(boost::factory<BRkpet*>(), _1, 1) );
        ThObsF.addObsToFactory("BRppkm", boost::bind(boost::factory<BRppkm*>(), _1, 1) );
        ThObsF.addObsToFactory("BRp0kS", boost::bind(boost::factory<BRp0kS*>(), _1, 1) );
        ThObsF.addObsToFactory("BRp0kL", boost::bind(boost::factory<BRp0kL*>(), _1, 1) );
        ThObsF.addObsToFactory("BRkSet", boost::bind(boost::factory<BRkSet*>(), _1, 1) );
        ThObsF.addObsToFactory("BRppkS", boost::bind(boost::factory<BRppkS*>(), _1, 1) );
        ThObsF.addObsToFactory("BRppkL", boost::bind(boost::factory<BRppkL*>(), _1, 1) );
        ThObsF.addObsToFactory("BRDskpkS", boost::bind(boost::factory<BRDskpkS*>(), _1, 1) );
        ThObsF.addObsToFactory("BRDskpkL", boost::bind(boost::factory<BRDskpkL*>(), _1, 1) );
        ThObsF.addObsToFactory("BRDsppet", boost::bind(boost::factory<BRDsppet*>(), _1, 1) );
        ThObsF.addObsToFactory("BRpmkp", boost::bind(boost::factory<BRpmkp*>(), _1, 1) );
        ThObsF.addObsToFactory("BRp0kp", boost::bind(boost::factory<BRp0kp*>(), _1, 1) );
        ThObsF.addObsToFactory("BRDskpet", boost::bind(boost::factory<BRDskpet*>(), _1, 1) );
        ThObsF.addObsToFactory("BRDskpkS_kpkL", boost::factory<BRDskpkS_kpkL*>() );
        ThObsF.addObsToFactory("BRppkm_pmkp", boost::factory<BRppkm_pmkp*>() );

        ThObsF.addObsToFactory("delta_kpi", boost::factory<delta_kpi*>() );

        ThObsF.addObsToFactory("aCPpppm", boost::bind(boost::factory<BRpppm*>(), _1, 2) );
        ThObsF.addObsToFactory("aCPp0p0", boost::bind(boost::factory<BRp0p0*>(), _1, 2) );
        ThObsF.addObsToFactory("aCPkpkm", boost::bind(boost::factory<BRkpkm*>(), _1, 2) );
        ThObsF.addObsToFactory("aCPkSkS", boost::bind(boost::factory<BRkSkS*>(), _1, 2) );
        ThObsF.addObsToFactory("aCPp0et", boost::bind(boost::factory<BRp0et*>(), _1, 2) );
        ThObsF.addObsToFactory("aCPkpkS", boost::bind(boost::factory<BRkpkS*>(), _1, 2) );
        ThObsF.addObsToFactory("aCPppet", boost::bind(boost::factory<BRppet*>(), _1, 2) );

        ThObsF.addObsToFactory("aCPDsppkS", boost::bind(boost::factory<BRDsppkS*>(), _1, 2) );
        ThObsF.addObsToFactory("aCPDsp0kp", boost::bind(boost::factory<BRDsp0kp*>(), _1, 2) );

        ThObsF.addObsToFactory("DACP", boost::factory<DACP*>() );

        ThObsF.addObsToFactory("GR_R1", boost::factory<GR_R1*>() );
        ThObsF.addObsToFactory("GR_R2", boost::factory<GR_R2*>() );
        ThObsF.addObsToFactory("GR_R3", boost::factory<GR_R3*>() );
        ThObsF.addObsToFactory("GR_R4", boost::factory<GR_R4*>() );

        ThObsF.addObsToFactory("GR_DR", boost::factory<GR_DR*>() );
        ThObsF.addObsToFactory("GR_ReEps1", boost::factory<GR_ReEps1*>() );
        ThObsF.addObsToFactory("GR_ReEps2", boost::factory<GR_ReEps2*>() );
        ThObsF.addObsToFactory("GR_ReEps0p", boost::factory<GR_ReEps0p*>() );
        ThObsF.addObsToFactory("GR_RD0", boost::factory<GR_RD0*>() );
        ThObsF.addObsToFactory("GR_RDp", boost::factory<GR_RDp*>() );
        ThObsF.addObsToFactory("GR_RDsp", boost::factory<GR_RDsp*>() );

        ThObsF.addObsToFactory("A1_T", boost::factory<A1_T*>() );
        ThObsF.addObsToFactory("A2_T", boost::factory<A2_T*>() );

        ThObsF.addObsToFactory("delta_1_s", boost::factory<delta_1_s*>() );
        ThObsF.addObsToFactory("delta_12_s", boost::factory<delta_12_s*>() );

        ThObsF.addObsToFactory("AbsRD_8_6_1",      boost::bind(boost::factory<RD_8_6_1*>(),      _1, 1) );
        ThObsF.addObsToFactory("AbsRD_8_15_1_CF",  boost::bind(boost::factory<RD_8_15_1_CF*>(),  _1, 1) );
        ThObsF.addObsToFactory("AbsRD_27_15_1_CF",    boost::factory<RD_27_15_1_CF*>() );
        ThObsF.addObsToFactory("AbsRD_8_6_0",      boost::bind(boost::factory<RD_8_6_0*>(),      _1, 1) );
        ThObsF.addObsToFactory("AbsRD_8_15_1_DCS", boost::bind(boost::factory<RD_8_15_1_DCS*>(), _1, 1) );
        ThObsF.addObsToFactory("AbsRD_27_15_1_DCS",    boost::factory<RD_27_15_1_DCS*>() );
        ThObsF.addObsToFactory("AbsRD_1_3_12",     boost::bind(boost::factory<RD_1_3_12*>(),     _1, 1) );
        ThObsF.addObsToFactory("AbsRD_8_3_12",     boost::bind(boost::factory<RD_8_3_12*>(),     _1, 1) );
        ThObsF.addObsToFactory("AbsRD_8_6_12",     boost::bind(boost::factory<RD_8_6_12*>(),     _1, 1) );
        ThObsF.addObsToFactory("AbsRD_8_15_12",    boost::bind(boost::factory<RD_8_15_12*>(),    _1, 1) );
        ThObsF.addObsToFactory("AbsRD_8_15_32",    boost::bind(boost::factory<RD_8_15_32*>(),    _1, 1) );
        ThObsF.addObsToFactory("AbsRD_27_15_12",    boost::factory<RD_27_15_12*>() );
        ThObsF.addObsToFactory("AbsRD_27_15_32",    boost::factory<RD_27_15_32*>() );

        ThObsF.addObsToFactory("ArgRD_8_6_1",      boost::bind(boost::factory<RD_8_6_1*>(),      _1, 2) );
        ThObsF.addObsToFactory("ArgRD_8_15_1_CF",  boost::bind(boost::factory<RD_8_15_1_CF*>(),  _1, 2) );
        ThObsF.addObsToFactory("ArgRD_8_6_0",      boost::bind(boost::factory<RD_8_6_0*>(),      _1, 2) );
        ThObsF.addObsToFactory("ArgRD_8_15_1_DCS", boost::bind(boost::factory<RD_8_15_1_DCS*>(), _1, 2) );
        ThObsF.addObsToFactory("ArgRD_1_3_12",     boost::bind(boost::factory<RD_1_3_12*>(),     _1, 2) );
        ThObsF.addObsToFactory("ArgRD_8_3_12",     boost::bind(boost::factory<RD_8_3_12*>(),     _1, 2) );
        ThObsF.addObsToFactory("ArgRD_8_6_12",     boost::bind(boost::factory<RD_8_6_12*>(),     _1, 2) );
        ThObsF.addObsToFactory("ArgRD_8_15_12",    boost::bind(boost::factory<RD_8_15_12*>(),    _1, 2) );
        ThObsF.addObsToFactory("ArgRD_8_15_32",    boost::bind(boost::factory<RD_8_15_32*>(),    _1, 2) );


        /* Create an object of the class MonteCarlo. */
        MonteCarlo MC(ModelF, ThObsF, ModelConf, MCMCConf, FileOut, JobTag);

        /* Do a test run if you wish to see the values of the observables      */
        /* and the correlated Gaussian observables defined in the model        */
        /* configuration file computed with the central value of the mandatory */
        /* parameters defined in the same file.                                */
        if (MCMCConf.compare("--noMC") == 0) MC.TestRun(rank);


        /* Initiate the Mote Carlo run. */
        else MC.Run(rank);

        // if (MCMCConf.compare("--noMC") != 0 && rank == 0){

        //     FileStat_t info;
        //     if (gSystem->GetPathInfo(directory.c_str(), info) == 0)
        //         gSystem->Exec(("rm -rf " + directory).c_str());
        //     gSystem->MakeDirectory(directory.c_str());
        //     gSystem->MakeDirectory((directory + "/code").c_str());
        //     gSystem->MakeDirectory((directory + "/code/src").c_str());
        //     gSystem->Exec(("mv *" + JobTag + "* " + directory).c_str());
        //     gSystem->Exec(("cp " + ModelConf + " " + directory + "/code/.").c_str());
        //     gSystem->Exec(("cp " + MCMCConf + " " + directory + "/code/.").c_str());
        //     gSystem->Exec(("cp -r ../src/*.cpp ../src/*.h " + directory + "/code/src/.").c_str());
        //     gSystem->Exec(("cp ../myModel_MCMC.cpp ../Makefile " + directory + "/code/.").c_str());
        // }

/* Necessary if MPI support is enabled during compilation. */
#ifdef _MPI
        MPI_Finalize();
#endif

        return EXIT_SUCCESS;
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }
}
