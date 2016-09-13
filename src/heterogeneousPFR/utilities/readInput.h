/*##############################################################################################
#                                                                                              #
#     #############       #############       #############       ####                ####     #
#    #             #     #             #     #             #     #    #              #    #    #
#    #    #####    #     #    #########      #    #####    #     #    #              #    #    #
#    #    #   #    #     #    #              #    #   #    #     #    #              #    #    #
#    #    #####    #     #    #              #    #####    #     #    #              #    #    #
#    #             #     #    #########      #             #     #    #              #    #    #
#    #             #     #             #     #             #     #    #              #    #    #
#    #    #####    #      #########    #     #    #####    #     #    #              #    #    #
#    #    #   #    #              #    #     #    #   #    #     #    #              #    #    #
#    #    #   #    #      #########    #     #    #   #    #     #    #########      #    #    #
#    #    #   #    #     #             #     #    #   #    #     #             #     #    #    #
#     ####     ####       #############       ####     ####       #############       ####     #
#                                                                                              #
#   Department of Energy                                                                       #
#   Politecnico di Milano                                                                      #
#   Author: Stefano Rebughini <stefano.rebughini@polimi.it>                                    #
#                                                                                              #
################################################################################################
#                                                                                              #
#   License                                                                                    #
#                                                                                              #
#   This file is part of ASALI.                                                                #
#                                                                                              #
#   ASALI is free software: you can redistribute it and/or modify                              #
#   it under the terms of the GNU General Public License as published by                       #
#   the Free Software Foundation, either version 3 of the License, or                          #
#   (at your option) any later version.                                                        #
#                                                                                              #
#   ASALI is distributed in the hope that it will be useful,                                   #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of                             #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                              #
#   GNU General Public License for more details.                                               #
#                                                                                              #
#   You should have received a copy of the GNU General Public License                          #
#   along with ASALI. If not, see <http://www.gnu.org/licenses/>.                              #
#                                                                                              #
##############################################################################################*/

namespace ASALI
{
    class READinput 
    {
        public:

            READinput(std::string& file, std::string& options);

            #include "convert.h"

            inline double getAbsTol()              const  { return absTol_;};
            inline double getRelTol()              const  { return relTol_;};
            inline double getIntegrationTime()     const  { return tMax_;};
            inline double getGridError()           const  { return errorOnGrid_;};
            inline double getEnergyError()         const  { return errorOnEnergy_;};
            inline double getInertLength()         const  { return Linert_;};
            inline double getCatalyticLength()     const  { return Lcat_;};
            inline double getReactorLength()       const  { return Lcat_ + Linert_;};
            inline double getChannelDiameter()     const  { return Dchannel_;};
            inline double getParticleDiameter()    const  { return Dp_;};
            inline double getHoneycombDiameter()   const  { return Dmatrix_;};
            inline double getTubeDiameter()        const  { return Dmatrix_;};
            inline double getVoidFraction()        const  { return epsi_;};
            inline double getAlfa()                const  { return alfa_;};
            inline double getPressure()            const  { return p_;};
            inline double getGasTemperature()      const  { return Tgas_;};
            inline double getSolidTemperature()    const  { return Tsolid_;};
            inline double getVelocity()            const  { return v_;};
            inline double getSolidDensity()        const  { return rhoSolid_;};
            inline double getSolidConductivity()   const  { return condSolid_;};
            inline double getSolidSpecificHeat()   const  { return cpSolid_;};
            inline double getSpecificArea()        const  { return av_;};
            inline double getExternalArea()        const  { return aex_;};
            inline double getExternalTemperature() const  { return Tex_;};

            inline unsigned int getRestartPoints()         const  { return int(restartN_);};
            inline unsigned int getMaxPointsNumber()       const  { return int(maxN_);};
            inline unsigned int getMinPointsNumber()       const  { return int(minN_);};
            inline unsigned int getAddPointsNumber()       const  { return int(addN_);};
            inline unsigned int getSamplingGridDimension() const  { return int(NS_);};

            inline std::string getResults()              const  { return results_;};
            inline std::string getFeed()                 const  { return feed_;};
            inline std::string getStart()                const  { return start_;};
            inline std::string getChannelShape()         const  { return shape_;};
            inline std::string getLimitations()          const  { return limitations_;};
            inline std::string getInert()                const  { return inert_;};
            inline std::string getRestartResults()       const  { return restartResults_;};
            inline std::string getSolver()               const  { return solver_;};
            inline std::string getKineticsPath()         const  { return kineticsPath_;};
            inline std::string getDiscretizationScheme() const  { return discretizationScheme_;};
            inline std::string getReactorType()          const  { return reactorType_;};
            inline std::string getCorrelation()          const  { return limitations_;};
            inline std::string getCatalyst()             const  { return catalyst_;};

            inline bool getEnergy()                const  { return energy_;};
            inline bool getExternalHeatExchange()  const  { return ex_;};
            inline bool getHomogeneousReactions()  const  { return homo_;};
            inline bool getHeterogenousReactions() const  { return het_;};
            inline bool getGridType()              const  { return grow_;};
            inline bool getConstraints()           const  { return constraints_;};
            inline bool getDiffusion()             const  { return gasDiffusion_;};

            inline std::vector<double> getFraction()                const { return inletValue_;};
            inline std::vector<double> getRestartBulk()             const { return restartBulk_;};
            inline std::vector<double> getRestartWall()             const { return restartWall_;};
            inline std::vector<double> getRestartSite()             const { return restartSite_;};
            inline std::vector<double> getRestartGasTemperature()   const { return restartTgas_;};
            inline std::vector<double> getRestartSolidTemperature() const { return restartTsolid_;};
            inline std::vector<double> getRestartGrid()             const { return restartZ_;};
            inline std::vector<double> getSamplingGrid()            const { return samplingGrid_;};

            inline std::vector<std::string> getFractionName() const { return inletName_;};

            double getSherwood();
            void recapOnScreen();

        private:

            double minN_;
            double maxN_;
            double addN_;
            double restartN_;
            double absTol_;
            double relTol_;
            double tMax_;
            double errorOnGrid_;
            double errorOnEnergy_;
            double Lcat_;
            double Linert_;
            double Dchannel_;
            double Dp_;
            double Dmatrix_;
            double epsi_;
            double alfa_;
            double Sh_;
            double p_;
            double Tgas_;
            double Tsolid_;
            double Tex_;
            double v_;
            double av_;
            double aex_;
            double rhoSolid_;
            double cpSolid_;
            double condSolid_;

            std::string results_;
            std::string feed_;
            std::string start_;
            std::string shape_;
            std::string limitations_;
            std::string inert_;
            std::string latestGrid_;
            std::string restartResults_;
            std::string solver_;
            std::string kineticsPath_;
            std::string discretizationScheme_;
            std::string reactorType_;
            std::string catalyst_;

            bool energy_;
            bool constraints_;
            bool homo_;
            bool het_;
            bool grow_;
            bool gasDiffusion_;
            bool ex_;

            const std::string& file_;
            const std::string& options_;

            std::vector<double>      inletValue_;
            std::vector<double>      restartZ_;
            std::vector<double>      restartTgas_;
            std::vector<double>      restartTsolid_;
            std::vector<double>      restartBulk_;
            std::vector<double>      restartWall_;
            std::vector<double>      restartSite_;
            std::vector<double>      samplingGrid_;

            std::vector<std::string> inputVector_;
            std::vector<std::string> inletName_;

            unsigned int reactorIndex_;
            unsigned int solverIndex_;
            unsigned int catalystIndex_;
            unsigned int temperatureIndex_;
            unsigned int fractionIndex_;
            unsigned int solidIndex_;
            unsigned int numericalIndex_;
            unsigned int kineticsIndex_;
            unsigned int NS_;

            void error() { std::cout << "\nASALI::READinput::ERROR\n" << std::endl;};

            void save();
            void check();
            void solver();
            void kinetics();
            void grid();
            void reactor();
            void catalyst();
            void fractions();
            void pressure();
            void temperature();
            void velocity();
            void solid();
            void restarting();
            void sampling();

            std::string boolOnScreen(bool value);

            template < typename T > std::string to_string( const T& v )
            {
                std::ostringstream stm ;
                return ( stm << v ) ? stm.str() : "{*** error ***}" ;
            }
    };
    
    READinput::READinput(std::string& file, std::string& options):
    file_(file),
    options_(options)
    {
        NS_                = 0;
        av_                = 0;
        aex_               = 0;
        minN_              = 0;
        maxN_              = 0;
        addN_              = 0;
        restartN_          = 0;
        Dp_                = 0;
        absTol_            = 0;
        relTol_            = 0;
        tMax_              = 0;
        errorOnGrid_       = 0;
        errorOnEnergy_     = 0;
        Lcat_              = 0;
        Linert_            = 0;
        Dchannel_          = 0;
        Dmatrix_           = 0;
        epsi_              = 0;
        alfa_              = 0;
        Sh_                = 0;
        p_                 = 0;
        Tgas_              = 0;
        Tsolid_            = 0;
        Tex_               = 0;
        v_                 = 0;
        rhoSolid_          = 0;
        cpSolid_           = 0;
        condSolid_         = 0;

        energy_            = false;
        constraints_       = false;
        homo_              = false;
        het_               = false;
        grow_              = false;
        gasDiffusion_      = false;
        ex_                = false;

        reactorIndex_      = 0;
        solverIndex_       = 0;
        catalystIndex_     = 0;
        temperatureIndex_  = 0;
        fractionIndex_     = 0;
        solidIndex_        = 0;
        numericalIndex_    = 0;
        kineticsIndex_     = 0;

        //- Check input
        {
            if ( boost::filesystem::exists(file_) && !boost::filesystem::is_directory(file_))
            {}
            else
            {
                error();
                std::cout << file_ << " MUST BE a file" << std::endl;
                exit(EXIT_FAILURE);
            }

            if ( options_ != "none" )
            {
                if ( boost::filesystem::exists(options_) && !boost::filesystem::is_directory(options_))
                {
                    if ( options_ != "samplingGrid.txt")
                    {
                        error();
                        std::cout << options_ << " file MUST BE renamed as || samplingGrid.txt || " << std::endl;
                        exit(EXIT_FAILURE);
                    }
                }
                else if ( boost::filesystem::exists(options_) && boost::filesystem::is_directory(options_))
                {
                    error();
                    std::cout << "1/ 'sampling':     " << options_ << " MUST BE a file named || samplingGrid.txt || " << std::endl;
                    std::cout << "2/ 'RPA':          " << options_ << " MUST BE a specie in the database " << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
        
        //- Save all 
        save();
        check();
        kinetics();
        solver();
        grid();
        reactor();
        catalyst();
        fractions();
        pressure();
        temperature();
        velocity();
        solid();
        if ( start_ == "latest"      || 
             start_ == "converter"   ||
             start_ == "kinetic"     ||
             start_ == "RPA"         ||
             start_ == "conversion"  ||
             start_ == "heat")
        {
            restarting();
        }
        else if ( start_ == "sampling" )
        {
            sampling();
            restarting();
        }
    }

    void READinput::save()
    {
        std::string dummyString;
        const char *path = file_.c_str();
        std::ifstream is(path);
        unsigned int k = 0;
        while (getline(is,dummyString))
        {
            while (getline(is,dummyString))
            {
                std::istringstream iss(dummyString);
                while (iss >> dummyString)
                {
                    k++;
                    inputVector_.resize(k);
                    inputVector_[k-1] = dummyString;
                }
            }
        }
    }

    void READinput::check()
    {
        std::vector<bool>         checkWord(9);
        std::vector<std::string>  words(9);

        for (unsigned int i=0;i<checkWord.size();i++)
            checkWord[i] = false;

        words[0] = "Temperature";
        words[1] = "Pressure";
        words[2] = "Reactor";
        words[3] = "Catalyst";
        words[4] = "Mole fractions || Mass fractions";
        words[5] = "Solver options";
        words[6] = "Volumetric flow rate || Velocity";
        words[7] = "Solid";
        words[8] = "Kinetics path";

        for (unsigned int i=0;i<inputVector_.size();i++)
        {
            if      (inputVector_[i]   == "Temperature")         {checkWord[0]  = true; temperatureIndex_  = i;}
            else if (inputVector_[i]   == "Pressure")            {checkWord[1]  = true;}
            else if (inputVector_[i]   == "Reactor")             {checkWord[2]  = true; reactorIndex_      = i;}
            else if (inputVector_[i]   == "Catalyst")            {checkWord[3]  = true; catalystIndex_     = i;}
            else if (inputVector_[i]   == "Mole" &&
                     inputVector_[i+1] == "fractions" )          {checkWord[4]  = true; fractionIndex_     = i;}
            else if (inputVector_[i]   == "Mass" &&
                     inputVector_[i+1] == "fractions" )          {checkWord[4]  = true; fractionIndex_     = i;}
            else if (inputVector_[i]   == "Solver" &&
                     inputVector_[i+1] == "options" )            {checkWord[5]  = true; solverIndex_       = i;}
            else if (inputVector_[i]   == "Volumetric" &&
                     inputVector_[i+1] == "flow" &&
                     inputVector_[i+2] == "rate" )               {checkWord[6]  = true;}
            else if (inputVector_[i]   == "Velocity")            {checkWord[6]  = true;}
            else if (inputVector_[i]   == "Solid")               {checkWord[7]  = true; solidIndex_        = i;}
            else if (inputVector_[i]   == "Kinetics" &&
                     inputVector_[i+1] == "path" )               {checkWord[8]  = true; kineticsIndex_     = i;}
        }
        
        for (unsigned int i=0;i<checkWord.size();i++)
        {
            if ( checkWord[i] == false )
            {
                error();
                std::cout << "key word || " << words[i] << " || " << "is MISSING!\n" << std::endl;
                exit (EXIT_FAILURE);
            }
        }
    }

    void READinput::solver()
    {
        if ( inputVector_[solverIndex_+1+1] != "{" )
        {
            error();
            std::cout << "The Solver options sub-dictionary must start with {\n" << std::endl;
            exit (EXIT_FAILURE);
        }
        else
        {
            unsigned int finalCount = 1e05;
            for (unsigned int i=solverIndex_;i<inputVector_.size();i++)
            {
                if (inputVector_[i] == "}")
                {
                    finalCount = i;
                    break;
                }
            }
            
            if ( finalCount == 1e05 )
            {
                error();
                std::cout << "The Solver options sub-dictionary must finish with }\n" << std::endl;
                exit (EXIT_FAILURE);
            }

            std::vector<std::string> dummyVector;
            unsigned int k=0;
            for (unsigned int i=solverIndex_+1+1+1;i<=finalCount;i++)
            {
                if (inputVector_[i] == "{")
                {
                    error();
                    std::cout << "The Solver options sub-dictionary must finish with }\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
                else
                {
                    dummyVector.resize(k+1);
                    dummyVector[k] = inputVector_[i];
                    k++;
                }
            }

            std::vector<bool>        checkWord(12);
            std::vector<std::string> words(12);

            unsigned int absIndex;
            unsigned int relIndex;
            unsigned int tIndex;
            unsigned int energyIndex;
            unsigned int resultsIndex;
            unsigned int reactionIndex;
            unsigned int constraintIndex;
            unsigned int errorIndex;
            unsigned int diffIndex;
            unsigned int solverIndex;
            unsigned int exIndex;

            for (unsigned int i=0;i<checkWord.size();i++)
                checkWord[i] = false;

            words[0]  = "Absolute tollerance";
            words[1]  = "Relative tollerance";
            words[2]  = "Energy equation";
            words[3]  = "Reactions";
            words[4]  = "Grid";
            words[5]  = "Integration time";
            words[6]  = "Results";
            words[7]  = "Constraints";
            words[8]  = "Accepted errors";
            words[9]  = "Diffusion in gas";
            words[10] = "Numerical solver";
            words[11] = "External heat exchange";

            for (unsigned int i=0;i<dummyVector.size();i++)
            {
                if      (dummyVector[i] == "Absolute" &&
                         dummyVector[i+1] == "tollerance")         {checkWord[0]  = true; absIndex         = i;}
                else if (dummyVector[i] == "Relative" &&
                         dummyVector[i+1] == "tollerance")         {checkWord[1]  = true; relIndex         = i;}
                else if (dummyVector[i] == "Energy" &&
                         dummyVector[i+1] == "equation")           {checkWord[2]  = true; energyIndex      = i;}
                else if (dummyVector[i] == "Reactions")            {checkWord[3]  = true; reactionIndex    = i;}
                else if (dummyVector[i] == "Grid")                 {checkWord[4]  = true;}
                else if (dummyVector[i] == "Integration" &&
                         dummyVector[i+1] == "time")               {checkWord[5]  = true; tIndex           = i;}
                else if (dummyVector[i] == "Results")              {checkWord[6]  = true; resultsIndex     = i;}
                else if (dummyVector[i] == "Constraints")          {checkWord[7]  = true; constraintIndex  = i;}
                else if (dummyVector[i] == "Accepted" &&
                         dummyVector[i+1] == "errors")             {checkWord[8]  = true; errorIndex        = i;}
                else if (dummyVector[i] == "Diffusion" &&
                         dummyVector[i+1] == "in" &&
                         dummyVector[i+2] == "gas")                {checkWord[9]  = true; diffIndex         = i;}
                else if (dummyVector[i] == "Numerical" &&
                         dummyVector[i+1] == "solver")             {checkWord[10] = true; solverIndex       = i;}
                else if (dummyVector[i] == "External" &&
                         dummyVector[i+1] == "heat" &&
                         dummyVector[i+2] == "exchange")           {checkWord[11]  = true; exIndex          = i;}
            }

            for (unsigned int i=0;i<checkWord.size();i++)
            {
                if ( checkWord[i] == false)
                {
                    error();
                    std::cout << "key word || " << words[i] << " || is MISSING in Solver options sub-dictionary!\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
            }
            
            solver_ = dummyVector[solverIndex + 1 + 1];

            if (solver_ == "BzzMath" )
            {
                #if ASALI_USE_BZZ == 0
                error();
                std::cout << "key word || Numerical solver || cannot be || BzzMath || \n" << std::endl;
                exit(EXIT_FAILURE);
                #endif
            }
            else if ( solver_ == "OpenSMOKE" )
            {}
            else
            {
                error();
                std::cout << "key word || Numerical solver || could be || OpenSMOKE ||";
                #if ASALI_USE_BZZ == 1
                std::cout << " BzzMath ||";
                #endif
                std::cout<< std::endl;
                std::cout << "\n" << std::endl;
                exit(EXIT_FAILURE);
            }

            if ( dummyVector[exIndex+3] == "on" )
                ex_ = true;
            else if ( dummyVector[exIndex+3] == "true" )
                ex_ = true;
            else if ( dummyVector[exIndex+3] == "yes" )
                ex_ = true;
            else if ( dummyVector[exIndex+3] == "1" )
                ex_ = true;
            else
                ex_ = false;


            if ( dummyVector[diffIndex+3] == "on" )
                gasDiffusion_ = true;
            else if ( dummyVector[diffIndex+3] == "true" )
                gasDiffusion_ = true;
            else if ( dummyVector[diffIndex+3] == "yes" )
                gasDiffusion_ = true;
            else if ( dummyVector[diffIndex+3] == "1" )
                gasDiffusion_ = true;
            else
                gasDiffusion_ = false;

            if ( dummyVector[energyIndex+2] == "on" )
                energy_ = true;
            else if ( dummyVector[energyIndex+2] == "true" )
                energy_ = true;
            else if ( dummyVector[energyIndex+2] == "yes" )
                energy_ = true;
            else if ( dummyVector[energyIndex+2] == "1" )
                energy_ = true;
            else
                energy_ = false;

            if ( dummyVector[constraintIndex+1] == "on" )
                constraints_ = true;
            else if ( dummyVector[constraintIndex+1] == "true" )
                constraints_ = true;
            else if ( dummyVector[constraintIndex+1] == "yes" )
                constraints_ = true;
            else if ( dummyVector[constraintIndex+1] == "1" )
                constraints_ = true;
            else
                constraints_ = false;

            results_ = dummyVector[resultsIndex+1];
            if ( results_ != "mass" && results_ != "mole" )
            {
                error();
                std::cout << "key word || " << "Results" << " || MUST be || mass || mole ||\n" << std::endl;
                exit (EXIT_FAILURE);
            }
            
            tMax_ = boost::lexical_cast<double>(dummyVector[tIndex+2]);
            std::string tdim = dummyVector[tIndex+3];
            ConvertsToSecond(tMax_, tdim);

            if ( dummyVector[reactionIndex+1] != "(" )
            {
                error();
                std::cout << "The Solver options/Reactions sub-dictionary must start with (\n" << std::endl;
                exit (EXIT_FAILURE);
            }
            else
            {
                unsigned int finalCount = 1e05;
                for (unsigned int i=reactionIndex;i<dummyVector.size();i++)
                {
                    if (dummyVector[i] == ")")
                    {
                        finalCount = i;
                        break;
                    }
                }
                
                if ( finalCount == 1e05 )
                {
                    error();
                    std::cout << "The Solver options/Reactions sub-dictionary must finish with )\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
                else
                {
                    unsigned int k=0;
                    for (unsigned int i=reactionIndex+1+1;i<finalCount;i++)
                    {
                        if (dummyVector[i] == "(")
                        {
                            error();
                            std::cout << "The Solver options/Reactions sub-dictionary must finish with )\n" << std::endl;
                            exit (EXIT_FAILURE);
                        }
                        else
                        {
                            if ( dummyVector[i] == "Homogeneous" )
                            {
                                if ( dummyVector[i+1] == "on" )
                                    homo_ = true;
                                else if ( dummyVector[i+1] == "true" )
                                    homo_ = true;
                                else if ( dummyVector[i+1] == "yes" )
                                    homo_ = true;
                                else if ( dummyVector[i+1] == "1" )
                                    homo_ = true;
                                else
                                    homo_ = false;
                                
                                i++;
                            }
                            else if ( dummyVector[i] == "Heterogeneous" )
                            {
                                if ( dummyVector[i+1] == "on" )
                                    het_ = true;
                                else if ( dummyVector[i+1] == "true" )
                                    het_ = true;
                                else if ( dummyVector[i+1] == "yes" )
                                    het_ = true;
                                else if ( dummyVector[i+1] == "1" )
                                    het_ = true;
                                else
                                    het_ = false;
                                
                                i++;
                            }
                            else
                            {
                                error();
                                std::cout << "key word || Homogeneous || Heterogeneous || is MISSING in Solver options/Reactions sub-dictionary!\n" << std::endl;
                                exit (EXIT_FAILURE);
                            }
                        }
                    }
                }
            }
            
            
            absTol_ = boost::lexical_cast<double>(dummyVector[absIndex+1+1]);
            relTol_ = boost::lexical_cast<double>(dummyVector[relIndex+1+1]);

            if ( dummyVector[errorIndex+1+1] != "(" )
            {
                error();
                std::cout << "The Solver options/Accepted errors sub-dictionary must start with (\n" << std::endl;
                exit (EXIT_FAILURE);
            }
            else
            {
                unsigned int finalCount = 1e05;
                for (unsigned int i=errorIndex;i<dummyVector.size();i++)
                {
                    if (dummyVector[i] == ")")
                    {
                        finalCount = i;
                        break;
                    }
                }
                
                if ( finalCount == 1e05 )
                {
                    error();
                    std::cout << "The Solver options/Accepted errors sub-dictionary must finish with )\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
                else
                {
                    unsigned int k=0;
                    for (unsigned int i=errorIndex+1+1+1;i<finalCount;i++)
                    {
                        if (dummyVector[i] == "(")
                        {
                            error();
                            std::cout << "The Solver options/Accepted errors sub-dictionary must finish with )\n" << std::endl;
                            exit (EXIT_FAILURE);
                        }
                        else
                        {
                            if ( dummyVector[i] == "grid" )
                            {
                                errorOnGrid_ = boost::lexical_cast<double>(dummyVector[i+1]);
                                i++;
                            }
                            else if ( dummyVector[i] == "energy" )
                            {
                                errorOnEnergy_ = boost::lexical_cast<double>(dummyVector[i+1]);
                                i++;
                            }
                            else
                            {
                                error();
                                std::cout << "key word || grid || energy || is MISSING in Solver options/Accepted errors sub-dictionary!\n" << std::endl;
                                exit (EXIT_FAILURE);
                            }
                        }
                    }
                }
            }
        }
    }

    void READinput::grid()
    {
        unsigned int finalCount = 1e05;
        for (unsigned int i=solverIndex_;i<inputVector_.size();i++)
        {
            if (inputVector_[i] == "}")
            {
                finalCount = i;
                break;
            }
        }
        std::vector<std::string> dummyVector;
        unsigned int k = 0;
        for (unsigned int i=solverIndex_+1+1+1;i<=finalCount;i++)
        {
            dummyVector.resize(k+1);
            dummyVector[k] = inputVector_[i];
            k++;
        }

        double gridIndex;
        for (unsigned int i=0;i<dummyVector.size();i++)
        {
            if (dummyVector[i] == "Grid")                    {gridIndex     = i;}
        }

        if ( dummyVector[gridIndex+1] != "(" )
        {
            error();
            std::cout << "The Solver options/Grid sub-dictionary must start with (\n" << std::endl;
            exit (EXIT_FAILURE);
        }
        else
        {
            finalCount = 1e05;
            for (unsigned int i=gridIndex;i<dummyVector.size();i++)
            {
                if (dummyVector[i] == ")")
                {
                    finalCount = i;
                    break;
                }
            }

            if ( finalCount == 1e05 )
            {
                error();
                std::cout << "The Solver options/Grid sub-dictionary must finish with )\n" << std::endl;
                exit (EXIT_FAILURE);
            }
            else
            {
                std::vector<std::string> vectorDummy;
                k = 0;
                for (unsigned int i=gridIndex+1+1;i<=finalCount;i++)
                {
                    if (dummyVector[i] == "(")
                    {
                        error();
                        std::cout << "The Solver options/Grid sub-dictionary must finish with )\n" << std::endl;
                        exit (EXIT_FAILURE);
                    }
                    else
                    {
                        vectorDummy.resize(k+1);
                        vectorDummy[k] = dummyVector[i];
                        k++;
                    }
                }

                std::vector<bool>        checkWord(6);
                std::vector<std::string> words(6);

                double growIndex;
                double minIndex;
                double maxIndex;
                double addIndex;
                double startIndex;
                double discreIndex;

                for (unsigned int i=0;i<checkWord.size();i++)
                    checkWord[i] = false;

                words[0] = "Growing";
                words[1] = "Minimum number of points";
                words[2] = "Maximum number of points";
                words[3] = "Add number of points";
                words[4] = "Resolution type";
                words[5] = "Discretization scheme";

                for (unsigned int i=0;i<vectorDummy.size();i++)
                {
                    if      (vectorDummy[i] == "Growing")             {checkWord[0] = true; growIndex     = i;}
                    else if (vectorDummy[i] == "Minimum" &&
                             vectorDummy[i+1] == "number" &&
                             vectorDummy[i+2] == "of" &&
                             vectorDummy[i+3] == "points")            {checkWord[1] = true; minIndex      = i;}
                    else if (vectorDummy[i] == "Maximum" &&
                             vectorDummy[i+1] == "number" &&
                             vectorDummy[i+2] == "of" &&
                             vectorDummy[i+3] == "points")            {checkWord[2] = true; maxIndex      = i;}
                    else if (vectorDummy[i] == "Add" &&
                             vectorDummy[i+1] == "number" &&
                             vectorDummy[i+2] == "of" &&
                             vectorDummy[i+3] == "points")            {checkWord[3] = true; addIndex      = i;}
                    else if (vectorDummy[i] == "Resolution" &&
                             vectorDummy[i+1] == "type")              {checkWord[4] = true; startIndex    = i;}
                    else if (vectorDummy[i] == "Discretization" &&
                             vectorDummy[i+1] == "scheme")            {checkWord[5] = true; discreIndex    = i;}
                }
            
                for (unsigned int i=0;i<checkWord.size();i++)
                {
                    if ( checkWord[i] == false)
                    {
                        error();
                        std::cout << "key word || " << words[i] << " || is MISSING in Solver options/Grid sub-dictionary!\n" << std::endl;
                        exit (EXIT_FAILURE);
                    }
                }

                minN_ = boost::lexical_cast<double>(vectorDummy[minIndex+1+1+1+1]);
                if ( minN_ < 10 )
                {
                    error();
                    std::cout << "Minimum number of points MUST BE > 10" << std::endl;
                    exit (EXIT_FAILURE);
                }
                maxN_ = boost::lexical_cast<double>(vectorDummy[maxIndex+1+1+1+1]);
                addN_ = boost::lexical_cast<double>(vectorDummy[addIndex+1+1+1+1]);

                if ( vectorDummy[growIndex+1] == "on" )
                    grow_ = true;
                else if ( vectorDummy[growIndex+1] == "true" )
                    grow_ = true;
                else if ( vectorDummy[growIndex+1] == "yes" )
                    grow_ = true;
                else if ( vectorDummy[growIndex+1] == "1" )
                    grow_ = true;
                else
                    grow_ = false;

                start_ = vectorDummy[startIndex+1+1];
                if ( start_ != "new"        && 
                     start_ != "latest"     && 
                     start_ != "converter"  && 
                     start_ != "sampling"   && 
                     start_ != "help"       &&
                     start_ != "kinetic"    &&
                     start_ != "RPA"        &&
                     start_ != "conversion" &&
                     start_ != "heat")
                {
                    error();
                    std::cout << "key word || " << "Resolution type" << " || MUST be || new || latest || converter || sampling || help || kinetic || RPA || conversion || heat ||\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
                else if ( start_ == "help" )
                {
                    error();
                    std::cout << "'Resolution type' options:" << std::endl;
                    std::cout << "1/ 'new'         : starting from a new uniform grid" << std::endl;
                    std::cout << "2/ 'latest'      : starting from the latest solution in the results folder" << std::endl;
                    std::cout << "3/ 'converter'   : last results are converted from mass to mole fraction and viceversa" << std::endl;
                    std::cout << "4/ 'sampling'    : last results are sampled on a grid provided by the user in the additional file" << std::endl;
                    std::cout << "5/ 'kinetic'     : kinetic analysis of the results" << std::endl;
                    std::cout << "6/ 'RPA'         : reaction path analysis of the results" << std::endl;
                    std::cout << "7/ 'conversion'  : conversion of the specie is obtained from the last results" << std::endl;
                    std::cout << "8/ 'heat'        : heat of reaction is obtained from the last results\n" << std::endl;
                    exit (EXIT_FAILURE);
                }

                discretizationScheme_ = vectorDummy[discreIndex+1+1];
                if ( discretizationScheme_ != "BDS" && discretizationScheme_ != "CDS")
                {
                    error();
                    std::cout << "key word || " << "Discretization scheme" << " || MUST be || CDS || BDS ||\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
            }
        }
    }
    
    void READinput::reactor()
    {
        if ( inputVector_[reactorIndex_+1] != "{" )
        {
            error();
            std::cout << "The Reactor sub-dictionary must start with {\n" << std::endl;
            exit (EXIT_FAILURE);
        }
        else
        {
            unsigned int finalCount = 1e05;
            for (unsigned int i=reactorIndex_;i<inputVector_.size();i++)
            {
                if (inputVector_[i] == "}")
                {
                    finalCount = i;
                    break;
                }
            }
            
            if ( finalCount == 1e05 )
            {
                error();
                std::cout << "The Reactor sub-dictionary must finish with }\n" << std::endl;
                exit (EXIT_FAILURE);
            }

            std::vector<std::string> dummyVector;
            unsigned int k=0;
            for (unsigned int i=reactorIndex_+1+1;i<=finalCount;i++)
            {
                if (inputVector_[i] == "{")
                {
                    error();
                    std::cout << "The Reactor sub-dictionary must finish with }\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
                else
                {
                    dummyVector.resize(k+1);
                    dummyVector[k] = inputVector_[i];
                    k++;
                }
            }
            
            std::vector<bool>        checkWord(11);
            std::vector<std::string> words(11);

            unsigned int inertIndex;
            unsigned int catIndex;
            unsigned int hydraulicIndex;
            unsigned int diameterIndex;
            unsigned int shapeIndex;
            unsigned int epsiIndex;
            unsigned int limitationIndex;
            unsigned int typeIndex;
            unsigned int particleIndex;
            unsigned int tubeIndex;
            unsigned int washcoatIndex;

            for (unsigned int i=0;i<checkWord.size();i++)
                checkWord[i] = false;

            words[0]  = "inert length";
            words[1]  = "catalytic length";
            words[2]  = "hydraulic diameter";
            words[3]  = "channel diameter";
            words[4]  = "channel shape";
            words[5]  = "void fraction";
            words[6]  = "gas-to-solid correlation";
            words[7]  = "type";
            words[8]  = "particle diameter";
            words[9]  = "washcoat thickness";
            words[10] = "tube shape";

            for (unsigned int i=0;i<dummyVector.size();i++)
            {
                if      (dummyVector[i] == "inert" &&
                         dummyVector[i+1] == "length")             {checkWord[0]  = true; inertIndex       = i;}
                else if (dummyVector[i] == "catalytic" &&
                         dummyVector[i+1] == "length")             {checkWord[1]  = true; catIndex         = i;}
                else if (dummyVector[i] == "hydraulic" &&
                         dummyVector[i+1] == "diameter")           {checkWord[2]  = true; hydraulicIndex   = i;}
                else if (dummyVector[i] == "channel" &&
                         dummyVector[i+1] == "diameter")           {checkWord[3]  = true; diameterIndex    = i;}
                else if (dummyVector[i] == "channel" &&
                         dummyVector[i+1] == "shape")              {checkWord[4]  = true; shapeIndex       = i;}
                else if (dummyVector[i] == "void" &&
                         dummyVector[i+1] == "fraction")           {checkWord[5]  = true; epsiIndex        = i;}
                else if (dummyVector[i] == "gas-to-solid" &&
                         dummyVector[i+1] == "correlation")        {checkWord[6]  = true; limitationIndex  = i;}
                else if (dummyVector[i] == "type")                 {checkWord[7]  = true; typeIndex        = i;}
                else if (dummyVector[i] == "particle" &&
                         dummyVector[i+1] == "diameter")           {checkWord[8]  = true; particleIndex    = i;}
                else if (dummyVector[i] == "washcoat" &&
                         dummyVector[i+1] == "thickness")          {checkWord[9]  = true; washcoatIndex    = i;}
                else if (dummyVector[i] == "tube" &&
                         dummyVector[i+1] == "shape")              {checkWord[10] = true; tubeIndex        = i;}
            }

            {
                reactorType_ = dummyVector[typeIndex + 1];

                if ( reactorType_ == "honeyComb" )
                {
                    checkWord[8]  = true;
                    checkWord[9]  = true;
                    checkWord[10] = true;
                }
                else if ( reactorType_ == "packedBed" )
                {
                    checkWord[3]  = true;
                    checkWord[4]  = true;
                    checkWord[5]  = true;
                    checkWord[9]  = true;
                    checkWord[10] = true;
                }
                else if ( reactorType_ == "tubular" )
                {
                    checkWord[3] = true;
                    checkWord[4] = true;
                    checkWord[5] = true;
                    checkWord[8] = true;
                }
                else
                {
                    error();
                    std::cout << "key word ||" << " type " << "|| MUST be || honeyComb || packedBed || tubular ||\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
            }


            for (unsigned int i=0;i<checkWord.size();i++)
            {
                if ( checkWord[i] == false)
                {
                    error();
                    std::cout << "key word || " << words[i] << " || is MISSING in Reactor sub-dictionary!\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
            }
            
            Lcat_ = boost::lexical_cast<double>(dummyVector[catIndex+2]);
            std::string dimCat = dummyVector[catIndex+3];
            ConvertsToMeter(Lcat_, dimCat);
            
            Linert_ = boost::lexical_cast<double>(dummyVector[inertIndex+2]);
            std::string dimInert = dummyVector[inertIndex+3];
            ConvertsToMeter(Linert_, dimInert);

            Dmatrix_ = boost::lexical_cast<double>(dummyVector[hydraulicIndex+2]);
            std::string dimMatrix = dummyVector[hydraulicIndex+3];
            ConvertsToMeter(Dmatrix_, dimMatrix);

            if ( reactorType_ == "honeyComb" )
            {
                Dchannel_ = boost::lexical_cast<double>(dummyVector[diameterIndex+2]);
                std::string dimChannel = dummyVector[diameterIndex+3];
                ConvertsToMeter(Dchannel_, dimChannel);

                epsi_ = boost::lexical_cast<double>(dummyVector[epsiIndex+2]);
                
                av_   = 4.*epsi_/Dchannel_;
                aex_  = 4./Dmatrix_;

                shape_ = dummyVector[shapeIndex+1+1];
                if ( shape_ != "circle" && shape_ != "triangle" && shape_ != "square")
                {
                    error();
                    std::cout << "key word ||" << " channel shape " << "|| MUST be || circle || triangle || square ||\n" << std::endl;
                    exit (EXIT_FAILURE);
                }

                limitations_ = dummyVector[limitationIndex+1+1];
                if ( limitations_ != "massTransfer" && limitations_ != "chemicalRegime")
                {
                    error();
                    std::cout << "key word ||" << " gas-to-particle solid " << "|| MUST be || massTransfer || chemicalRegime ||\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
            }
            else if ( reactorType_ == "packedBed" )
            {
                Dp_ = boost::lexical_cast<double>(dummyVector[particleIndex+2]);
                std::string dimParticle = dummyVector[particleIndex+3];
                ConvertsToMeter(Dp_, dimParticle);

                epsi_ = 0.4 + 0.05*Dp_/Dmatrix_ + 0.4*std::pow((Dp_/Dmatrix_),2.);
                
                av_   = 6.*(1. - epsi_)/Dp_;
                aex_  = 4./Dmatrix_;

                limitations_ = dummyVector[limitationIndex+1+1];
                if ( limitations_ != "Yoshida" && limitations_ != "Wakao" && limitations_ != "Petrovic" )
                {
                    error();
                    std::cout << "key word ||" << " gas-to-particle solid " << "|| MUST be || Yoshida || Wakao || Petrovic || \n" << std::endl;
                    exit (EXIT_FAILURE);
                }

            }
            else if ( reactorType_ == "tubular" )
            {
                double Sw   = boost::lexical_cast<double>(dummyVector[washcoatIndex+2]);
                std::string dimSw = dummyVector[washcoatIndex+3];
                ConvertsToMeter(Sw, dimSw);

                double Dint = Dmatrix_ - 2.*Sw;

                epsi_ = std::pow(Dint,2.)/std::pow(Dmatrix_,2.);
                
                av_   = 4./Dmatrix_;
                aex_  = 4./Dmatrix_;

                shape_ = dummyVector[tubeIndex+1+1];
                if ( shape_ != "circle" && shape_ != "triangle" && shape_ != "square")
                {
                    error();
                    std::cout << "key word ||" << " channel shape " << "|| MUST be || circle || triangle || square ||\n" << std::endl;
                    exit (EXIT_FAILURE);
                }

                limitations_ = dummyVector[limitationIndex+1+1];
                if ( limitations_ != "massTransfer" && limitations_ != "chemicalRegime")
                {
                    error();
                    std::cout << "key word ||" << " gas-to-particle solid " << "|| MUST be || massTransfer || chemicalRegime ||\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
            }
        }
    }
    
    void READinput::catalyst()
    {
        if ( inputVector_[catalystIndex_+1] != "{" )
        {
            error();
            std::cout << "The Catalyst sub-dictionary must start with {\n" << std::endl;
            exit (EXIT_FAILURE);
        }
        else
        {
            unsigned int finalCount = 1e05;
            for (unsigned int i=catalystIndex_;i<inputVector_.size();i++)
            {
                if (inputVector_[i] == "}")
                {
                    finalCount = i;
                    break;
                }
            }
            
            if ( finalCount == 1e05 )
            {
                error();
                std::cout << "The Catalyst sub-dictionary must finish with }\n" << std::endl;
                exit (EXIT_FAILURE);
            }

            std::vector<std::string> dummyVector;
            unsigned int k=0;
            for (unsigned int i=catalystIndex_+1+1;i<=finalCount;i++)
            {
                if (inputVector_[i] == "{")
                {
                    error();
                    std::cout << "The Catalyst sub-dictionary must finish with }\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
                else
                {
                    dummyVector.resize(k+1);
                    dummyVector[k] = inputVector_[i];
                    k++;
                }
            }

            std::vector<bool>        checkWord(6);
            std::vector<std::string> words(6);

            unsigned int massIndex;
            unsigned int dispersionIndex;
            unsigned int RhIndex;
            unsigned int deactivationIndex;
            unsigned int alfaIndex;
            unsigned int activeIndex;

            for (unsigned int i=0;i<checkWord.size();i++)
                checkWord[i] = false;

            words[0] = "mass";
            words[1] = "dispersion";
            words[2] = "Rh fraction";
            words[3] = "deactivation factor";
            words[4] = "alfa";
            words[5] = "active site";


            for (unsigned int i=0;i<dummyVector.size();i++)
            {
                if         (dummyVector[i] == "mass")                 {checkWord[0] = true; massIndex          = i;}
                else if (dummyVector[i] == "dispersion")              {checkWord[1] = true; dispersionIndex    = i;}
                else if (dummyVector[i] == "Rh" &&
                         dummyVector[i+1] == "fraction")              {checkWord[2] = true; RhIndex            = i;}
                else if (dummyVector[i] == "deactivation" &&
                         dummyVector[i+1] == "factor")                {checkWord[3] = true; deactivationIndex  = i;}
                else if (dummyVector[i] == "active" &&
                         dummyVector[i+1] == "site")                  {checkWord[5] = true; activeIndex        = i;}
            }

            for (unsigned int i=0;i<dummyVector.size();i++)
            {
                if(dummyVector[i] == "alfa")
                {
                    for (unsigned int k=0;k<checkWord.size();k++)
                        checkWord[k] = true;

                    alfaIndex = i;
                    break;
                }
                else if (dummyVector[i] == "mass")                {checkWord[4] = true; massIndex           = i;}
                else if (dummyVector[i] == "dispersion")          {checkWord[4] = true; dispersionIndex     = i;}
                else if (dummyVector[i] == "Rh" &&
                         dummyVector[i+1] == "fraction")          {checkWord[4] = true; RhIndex             = i;}
                else if (dummyVector[i] == "deactivation" &&
                         dummyVector[i+1] == "factor")            {checkWord[4] = true; deactivationIndex    = i;}
            }

            for (unsigned int i=0;i<checkWord.size();i++)
            {
                if ( checkWord[i] == false)
                {
                    error();
                    std::cout << "key word || " << words[i] << " || is MISSING in Catalyst sub-dictionary!\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
            }
            
            catalyst_ = dummyVector[activeIndex+2];
            
            for (unsigned int i=0;i<dummyVector.size();i++)
            {
                if(dummyVector[i] == "alfa")
                {
                    alfa_ = boost::lexical_cast<double>(dummyVector[alfaIndex+1]);
                    std::string dim = dummyVector[alfaIndex+2];
                    ConvertsToOneOnMeter(alfa_,dim);
                    break;
                }
                else
                {
                    double Wcat = boost::lexical_cast<double>(dummyVector[massIndex+1]);
                    std::string dim = dummyVector[massIndex+2];
                    ConvertsToKg(Wcat, dim);
            
                    double RhDispersion      = boost::lexical_cast<double>(dummyVector[dispersionIndex+1]);
                    double RhMassFraction    = boost::lexical_cast<double>(dummyVector[RhIndex+2]);
                    double dectivationFactor = boost::lexical_cast<double>(dummyVector[deactivationIndex+2]);

                    double RhPM = 102.9;
                    double SiteD = 2.49e-08;
                    double ActiveArea = RhMassFraction*RhDispersion*Wcat/(RhPM*SiteD);

                    double Ain = 0.25*3.14*Dmatrix_*Dmatrix_;
                    double ReactorVolume = Ain*Lcat_;
                    alfa_ = dectivationFactor*ActiveArea/ReactorVolume;
                }
            }
        }
    }

    double READinput::getSherwood()
    {
        if (limitations_ == "chemicalRegime")
        {
            if ( shape_ == "square")                 {Sh_ = 3.087;}
            else if ( shape_ == "circle")            {Sh_ = 4.362;}
            else if ( shape_ == "triangle")          {Sh_ = 1.891;}
        }
        else if (limitations_ == "massTransfer")
        {
            if ( shape_ == "square")                {Sh_ = 2.977;}
            else if ( shape_ == "circle")           {Sh_ = 3.659;}
            else if ( shape_ == "triangle")         {Sh_ = 2.494;}
        }
        
        return Sh_;
    }

    void READinput::fractions()
    {
        if ( inputVector_[fractionIndex_+1+1] != "{" )
        {
            error();
            std::cout << "The Mole/Mass fractions sub-dictionary must start with {\n" << std::endl;
            exit (EXIT_FAILURE);
        }
        else
        {
            unsigned int finalCount = 1e05;
            for (unsigned int i=fractionIndex_;i<inputVector_.size();i++)
            {
                if (inputVector_[i] == "}")
                {
                    finalCount = i;
                    break;
                }
            }
            
            if ( finalCount == 1e05 )
            {
                error();
                std::cout << "The Mole/Mass fractions sub-dictionary must finish with }\n" << std::endl;
                exit (EXIT_FAILURE);
            }

            unsigned int k=0;
            for (unsigned int i=fractionIndex_+1+1+1;i<finalCount;i++)
            {
                if (inputVector_[i] == "{")
                {
                    error();
                    std::cout << "The Mole/Mass fractions sub-dictionary must finish with }\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
                else if ( inputVector_[i] == "inert" )
                {
                    inert_ = inputVector_[i+1];
                    i++;
                }
                else
                {
                    inletValue_.resize(k+1);
                    inletName_.resize(k+1);
                    inletValue_[k] = boost::lexical_cast<double>(inputVector_[i+1]);
                    inletName_[k]  = inputVector_[i];
                    k++;
                    i++;
                }
            }

            bool test = false;
            for (unsigned int i=0;i<inletName_.size();i++)
            {
                if ( inletName_[i] == inert_)
                {
                    test = true;
                }
            }
            
            if ( test == false )
            {
                error();
                std::cout << "key word || inert || MUST BE one of the feed species!\n" << std::endl;
                exit (EXIT_FAILURE);
            }

            feed_ = inputVector_[fractionIndex_];
        }
    }

    void READinput::pressure()
    {
        for (unsigned int i=0;i<inputVector_.size();i++)
        {
            if ( inputVector_[i] == "Pressure")
            {
                p_ = boost::lexical_cast<double>(inputVector_[i+1]);
                std::string dim = inputVector_[i+2];
                ConvertsToPascal(p_, dim);
            }
        }
    }

    void READinput::kinetics()
    {
        kineticsPath_ = inputVector_[kineticsIndex_+2];
    }

    void READinput::temperature()
    {
        if ( inputVector_[temperatureIndex_+1] != "{" )
        {
            error();
            std::cout << "The Temperature sub-dictionary must start with {\n" << std::endl;
            exit (EXIT_FAILURE);
        }
        else
        {
            unsigned int finalCount = 1e05;
            for (unsigned int i=temperatureIndex_;i<inputVector_.size();i++)
            {
                if (inputVector_[i] == "}")
                {
                    finalCount = i;
                    break;
                }
            }
            
            if ( finalCount == 1e05 )
            {
                error();
                std::cout << "The Temperature sub-dictionary must finish with }\n" << std::endl;
                exit (EXIT_FAILURE);
            }

            std::vector<std::string> dummyVector;
            unsigned int k=0;
            for (unsigned int i=temperatureIndex_+1+1;i<finalCount;i++)
            {
                if (inputVector_[i] == "{")
                {
                    error();
                    std::cout << "The Temperature sub-dictionary must finish with }\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
                else
                {
                    dummyVector.resize(k+1);
                    dummyVector[k] = inputVector_[i];
                    k++;
                }
            }

            std::vector<bool>        checkWord(3);
            std::vector<std::string> words(3);

            unsigned int gasIndex;
            unsigned int solidIndex;
            unsigned int exIndex;

            for (unsigned int i=0;i<checkWord.size();i++)
                checkWord[i] = false;

            words[0] = "gas";
            words[1] = "solid";
            words[2] = "external";


            for (unsigned int i=0;i<dummyVector.size();i++)
            {
                if      (dummyVector[i] == "gas")               {checkWord[0] = true; gasIndex    = i;}
                else if (dummyVector[i] == "solid")             {checkWord[1] = true; solidIndex  = i;}
                else if (dummyVector[i] == "external")          {checkWord[2] = true; exIndex     = i;}
            }

            if ( ex_ == false )
                checkWord[2] = true;


            for (unsigned int i=0;i<checkWord.size();i++)
            {
                if ( checkWord[i] == false)
                {
                    error();
                    std::cout << "key word || " << words[i] << " || alfa || is MISSING in Reactor sub-dictionary!\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
            }

            Tgas_ = boost::lexical_cast<double>(dummyVector[gasIndex+1]);
            std::string dimGas = dummyVector[gasIndex+2];
            if ( dimGas == "°C")
                FromCelsiusToKelvin(Tgas_,dimGas);

            Tsolid_ = boost::lexical_cast<double>(dummyVector[solidIndex+1]);
            std::string dimSolid = dummyVector[solidIndex+2];
            if ( dimSolid == "°C")
                FromCelsiusToKelvin(Tsolid_,dimSolid);

            if ( ex_ == true )
            {
                Tex_ = boost::lexical_cast<double>(dummyVector[exIndex+1]);
                std::string dimEx = dummyVector[exIndex+2];
                if ( dimEx == "°C")
                    FromCelsiusToKelvin(Tex_,dimEx);
            }
            
        }
    }

    void READinput::velocity()
    {
        for (unsigned int i=0;i<inputVector_.size();i++)
        {
            if (inputVector_[i] == "Velocity")
            {
                v_ = boost::lexical_cast<double>(inputVector_[i+1]);
                std::string dim = inputVector_[i+2];
                ConvertsToMeterPerSecond(v_,dim);
            }
            else if (inputVector_[i]   == "Volumetric" &&
                     inputVector_[i+1] == "flow" &&
                     inputVector_[i+2] == "rate" )
            {
                v_ = boost::lexical_cast<double>(inputVector_[i+3]);
                std::string dim = inputVector_[i+4];
                ConvertsToNm3perSecond(v_, dim);
                double Ain = 0.25*3.14*pow(Dmatrix_,2.);
                double T0 = 273.15;
                double P0 = 101325.;
                double nR = v_*P0/T0;
                double vnew = nR*Tgas_/p_;
                v_ = vnew/Ain;
            }
        }
    }
    
    void READinput::solid()
    {
        if ( inputVector_[solidIndex_+1] != "{" )
        {
            error();
            std::cout << "The Solid sub-dictionary must start with {\n" << std::endl;
            exit (EXIT_FAILURE);
        }
        else
        {
            unsigned int finalCount = 1e05;
            for (unsigned int i=solidIndex_;i<inputVector_.size();i++)
            {
                if (inputVector_[i] == "}")
                {
                    finalCount = i;
                    break;
                }
            }
            
            if ( finalCount == 1e05 )
            {
                error();
                std::cout << "The Solid sub-dictionary must finish with }\n" << std::endl;
                exit (EXIT_FAILURE);
            }

            std::vector<std::string> dummyVector;
            unsigned int k=0;
            for (unsigned int i=solidIndex_+1+1;i<=finalCount;i++)
            {
                if (inputVector_[i] == "{")
                {
                    error();
                    std::cout << "The Solid sub-dictionary must finish with }\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
                else
                {
                    dummyVector.resize(k+1);
                    dummyVector[k] = inputVector_[i];
                    k++;
                }
            }

            std::vector<bool>        checkWord(3);
            std::vector<std::string> words(3);

            double rhoIndex;
            double condIndex;
            double cpIndex;

            for (unsigned int i=0;i<checkWord.size();i++)
                checkWord[i] = false;

            words[0] = "density";
            words[1] = "conductivity";
            words[2] = "specific heat";


            for (unsigned int i=0;i<dummyVector.size();i++)
            {
                if         (dummyVector[i] == "density")                 {checkWord[0] = true; rhoIndex      = i;}
                else if (dummyVector[i] == "conductivity")             {checkWord[1] = true; condIndex     = i;}
                else if (dummyVector[i] == "specific" &&
                         dummyVector[i+1] == "heat")                {checkWord[2] = true; cpIndex       = i;}
            }

            for (unsigned int i=0;i<checkWord.size();i++)
            {
                if ( checkWord[i] == false)
                {
                    error();
                    std::cout << "key word || " << words[i] << " || alfa || is MISSING in Reactor sub-dictionary!\n" << std::endl;
                    exit (EXIT_FAILURE);
                }
            }
            
            if ( dummyVector[rhoIndex + 1 + 1] == "Kg/m3" )
            {
                rhoSolid_ = boost::lexical_cast<double>(dummyVector[rhoIndex + 1]);
            }
            else
            {
                error();
                std::cout << "key word || " << "density" << " || accepted unit is Kg/m3 in Solid sub-dictionary!\n" << std::endl;
                exit (EXIT_FAILURE);
            }

            if ( dummyVector[condIndex + 1 + 1] == "W/m/K" )
            {
                condSolid_ = boost::lexical_cast<double>(dummyVector[condIndex + 1]);
            }
            else
            {
                error();
                std::cout << "key word || " << "conductivity" << " || accepted unit is W/m/K in Solid sub-dictionary!\n" << std::endl;
                exit (EXIT_FAILURE);
            }

            if ( dummyVector[cpIndex + 1 + 1 + 1] == "J/Kg/K" )
            {
                cpSolid_ = boost::lexical_cast<double>(dummyVector[cpIndex + 1 + 1]);
            }
            else
            {
                error();
                std::cout << "key word || " << "specific heat" << " || accepted unit is J/Kg/K in Solid sub-dictionary!\n" << std::endl;
                exit (EXIT_FAILURE);
            }

        }
    }

    void READinput::restarting()
    {
        std::string latestGridpath;
        std::string bulkSpeciespath;
        std::string adsorbedSpeciespath;
        std::string gridpath;
        std::string informationpath;
        std::string temperaturepath;
        std::string wallSpeciespath;

        for (unsigned int i=1;i<=1e05;i++)
        {
            std::string path = "results/" + to_string(i);
            if (boost::filesystem::exists(path))
            {
                if (boost::filesystem::is_directory(path))  
                {
                    latestGrid_         = to_string(i);
                    bulkSpeciespath     = path + "/bulkSpecies.txt";
                    adsorbedSpeciespath = path + "/adsorbedSpecies.txt";
                    gridpath            = path + "/grid.txt";
                    informationpath     = path + "/information.txt";
                    temperaturepath     = path + "/temperature.txt";
                    wallSpeciespath     = path + "/wallSpecies.txt";
                }
            }
        }

        if (boost::lexical_cast<double>(latestGrid_) >= maxN_ && start_ == "latest")
        {
            error();
            std::cout << "Max number of grid points already reached" << std::endl;
            exit (EXIT_FAILURE);
        }

        std::string dummyString;
        std::vector<double> dummyVector;
        double dummyDouble;

        if (boost::filesystem::exists(gridpath))
        {
            const char *pathGrid = gridpath.c_str();
            std::ifstream is(pathGrid);
            unsigned int k=0;
            while (getline(is,dummyString))
            {
                dummyVector.resize(k+1);
                dummyVector[k] = boost::lexical_cast<double>(dummyString);
                k++;
            }

            restartN_ = dummyVector.size();
            restartZ_.resize(restartN_);

            for (unsigned int k=0;k<restartZ_.size();k++)
                restartZ_[k] = dummyVector[k];
        }
        else
        {
            error();
            std::cout << "grid.txt file NOT exist" << std::endl;
            exit (EXIT_FAILURE);
        }

        if (boost::filesystem::exists(informationpath))
        {
            const char *pathInfo = informationpath.c_str();
            std::ifstream is(pathInfo);
            while (getline(is,dummyString))
            {
                std::istringstream iss(dummyString);
                while (iss >> dummyString)
                {
                    restartResults_ = dummyString;
                    if ( restartResults_ == "mole" || restartResults_ == "mass")
                        break;
                }
                if ( restartResults_ == "mole" || restartResults_ == "mass")
                    break;
            }
        }
        else
        {
            error();
            std::cout << "information.txt file NOT exist" << std::endl;
            exit (EXIT_FAILURE);
        }

        if (boost::filesystem::exists(temperaturepath))
        {
            const char *pathTemperature = temperaturepath.c_str();
            std::ifstream is(pathTemperature);

            dummyVector.resize(restartN_*2);
            for (unsigned int k=0;k<dummyVector.size();k++)
                dummyVector[k] = 0.;

            unsigned int k=0;
            while (getline(is, dummyString))
            {
                std::istringstream iss(dummyString);
                while (iss >> dummyDouble)
                {
                    dummyVector[k] = dummyDouble;
                    k++;
                }
            }
            restartTgas_.resize(restartN_);
            restartTsolid_.resize(restartN_);
            unsigned int counter=0;
            for (unsigned int k=0;k<restartN_;k++)
            {
                restartTgas_[k]   = std::max(0.,dummyVector[counter]);
                restartTsolid_[k] = std::max(0.,dummyVector[counter + 1]);
                counter = counter + 2;
            }
        }
        else
        {
            error();
            std::cout << "temperature.txt file NOT exist" << std::endl;
            exit (EXIT_FAILURE);
        }

        if (boost::filesystem::exists(bulkSpeciespath))
        {
            const char *pathBulk = bulkSpeciespath.c_str();
            std::ifstream is(pathBulk);

            unsigned int k=0;
            while (getline(is, dummyString))
            {
                std::istringstream iss(dummyString);
                while (iss >> dummyDouble)
                {
                    restartBulk_.resize(k+1);
                    restartBulk_[k] = std::max(0.,dummyDouble);
                    k++;
                }
            }
        }
        else
        {
            error();
            std::cout << "bulkSpecies.txt file NOT exist" << std::endl;
            exit (EXIT_FAILURE);
        }

        if (boost::filesystem::exists(wallSpeciespath))
        {
            const char *pathWall = wallSpeciespath.c_str();
            std::ifstream is(pathWall);

            unsigned int k=0;
            while (getline(is, dummyString))
            {
                std::istringstream iss(dummyString);
                while (iss >> dummyDouble)
                {
                    restartWall_.resize(k+1);
                    restartWall_[k] = std::max(0.,dummyDouble);
                    k++;
                }
            }
        }
        else
        {
            error();
            std::cout << "wallSpecies.txt file NOT exist" << std::endl;
            exit (EXIT_FAILURE);
        }

        if (boost::filesystem::exists(adsorbedSpeciespath))
        {
            const char *pathSite = adsorbedSpeciespath.c_str();
            std::ifstream is(pathSite);

            unsigned int k=0;
            while (getline(is, dummyString))
            {
                std::istringstream iss(dummyString);
                while (iss >> dummyDouble)
                {
                    restartSite_.resize(k+1);
                    restartSite_[k] = std::max(0.,dummyDouble);
                    k++;
                }
            }
        }
        else
        {
            error();
            std::cout << "adsorbedSpecies.txt file NOT exist" << std::endl;
            exit (EXIT_FAILURE);
        }
    }

    void READinput::recapOnScreen()
    {
        std::cout.precision(6);
        std::cout << "\n################################################################################################" << std::endl;
        std::cout << "                                          GENERAL INPUT                                         \n" << std::endl;
        std::cout << "Reactor lenght                           = " << Lcat_ + Linert_ << "\t[m]" << std::endl;
        std::cout << "Inert lenght                             = " << Linert_ << "\t[m]" << std::endl;
        std::cout << "Catalytic lenght                         = " << Lcat_ << "\t[m]" << std::endl;
        if ( reactorType_ == "honeyComb" )
        {
            std::cout << "Channel diameter                         = " << Dchannel_ << "\t[m]" << std::endl;
            std::cout << "Honeycomb diameter                       = " << Dmatrix_ << "\t[m]" << std::endl;
        }
        else if ( reactorType_ == "packedBed" )
        {
            std::cout << "Particle diameter                        = " << Dp_ << "\t[m]" << std::endl;
            std::cout << "Tube diameter                            = " << Dmatrix_ << "\t[m]" << std::endl;
        }
        std::cout << "\n################################################################################################" << std::endl;
        std::cout << "                                       REACTOR PROPERTIES                                       \n" << std::endl;
        std::cout << "Specific area (av)                       = " << av_ << "\t[1/m]" << std::endl;
        std::cout << "Specific catalytic area (alfa)           = " << alfa_ << "\t[1/m]" << std::endl;
        std::cout << "Void fraction                            = " << epsi_ << "\t[-]" << std::endl;
        if ( reactorType_ == "honeyComb" )
        {
            std::cout << "Transport regime:                          " << limitations_ << std::endl;
            std::cout << "Channel shape:                             " << shape_ << std::endl;
        }
        std::cout << "\n################################################################################################" << std::endl;
        std::cout << "                                        SOLID PROPERTIES                                        \n" << std::endl;
        std::cout << "Density                                  = " << rhoSolid_ << "\t[Kg/m3]" << std::endl;
        std::cout << "Specific heat                            = " << cpSolid_ << "\t[J/Kg/K]" << std::endl;
        std::cout << "Conductivity                             = " << condSolid_ << "\t[W/m/K]" << std::endl;
        std::cout << "\n################################################################################################" << std::endl;
        std::cout << "                                      OPERATING CONDITIONS                                      \n" << std::endl;
        std::cout << "Feed velocity                            = " << v_ << "\t[m/s]" << std::endl;
        std::cout << "Feed pressure                            = " << p_ << "\t[Pa]" << std::endl;
        std::cout << "Initial solid temperature                = " << Tsolid_ << "\t[K]" << std::endl;
        std::cout << "Feed temperature                         = " << Tgas_ << "\t[K]" << std::endl;
        std::cout << "\n################################################################################################" << std::endl;
        std::cout << "                                         SOLVER OPTIONS                                         \n" << std::endl;
        std::cout << "Homogeneous  reactions:                    " << boolOnScreen(homo_) << std::endl;
        std::cout << "Heterogeneus reactions:                    " << boolOnScreen(het_) << std::endl;
        std::cout << "Energy balance:                            " << boolOnScreen(energy_) << std::endl;
        std::cout << "Feed in:                                   " << feed_ << " fractions" << std::endl;
        std::cout << "Results in:                                " << results_ << " fractions" << std::endl;
        std::cout << "\n################################################################################################" << std::endl;
        std::cout << "                                       AXIAL GRID OPTIONS                                       \n" << std::endl;
        std::cout << "Growing grid:                               " << boolOnScreen(grow_) << std::endl;
        if (grow_ == true)
        {

            if ( start_ == "new" )
                std::cout << "Min number of points                      = " << minN_ << std::endl;
            else
                std::cout << "Starting folder:                            " << latestGrid_ << std::endl;

            std::cout << "Max number of points                      = " << maxN_ << std::endl;
            std::cout << "Add points per iteration                  = " << addN_ << std::endl;
        }
        else
        {
            std::cout << "Number of points                          = " << minN_ << std::endl;
        }
        std::cout << "Integration time                          = " << tMax_ << "\t[s]" << std::endl;
        std::cout << "\n################################################################################################" << std::endl;
        std::cout << "                                        NUMERICAL SOLVER                                        \n" << std::endl;
        std::cout << "Solver compiled with                        ";
        #if ASALI_USE_BZZ == 1
        std::cout << "|| BzzMath ";
        #endif
        std::cout << "|| OpenSMOKE ";
        std::cout << "||" << std::endl;
        std::cout << "Chosen solver:                              || " << solver_ << " ||" << std::endl;
        std::cout << "\n################################################################################################" << std::endl;
    }
    
    std::string READinput::boolOnScreen(bool value)
    {
        std::string value_;
        if ( value == true)
            value_ = "on";
        else
            value_ = "off";
        
        return value_;
    }

    void READinput::sampling()
    {
        if ( options_ == "samplingGrid.txt" )
        {
            if ( boost::filesystem::is_empty(options_) )
            {
                
                error();
                std::cout << "No sampling points in || samplingGrid.txt || " << std::endl;
                exit (EXIT_FAILURE);
            }
            else
            {
                std::string dummyString;
                const char *path = options_.c_str();
                std::ifstream is(path);
                unsigned int k=0;
                while (getline(is,dummyString))
                {
                    k++;
                    samplingGrid_.resize(k);
                    samplingGrid_[k-1] = boost::lexical_cast<double>(dummyString);
                }

                NS_ = samplingGrid_.size();
            }


            if ( samplingGrid_[NS_- 1] > (Lcat_ + Linert_) )
            {
                error();
                std::cout << "Max value in || samplingGrid.txt || file is bigger than reactor length ( " << Lcat_ + Linert_ << " m )" << std::endl;
                exit (EXIT_FAILURE);
            }
            else if ( samplingGrid_[0] < 0. )
            {
                error();
                std::cout << "Min value in || samplingGrid.txt || file is negative ( reactor length > 0) " << std::endl;
                exit (EXIT_FAILURE);
            }
            
        }
    }
}
