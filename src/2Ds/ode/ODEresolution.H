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

if ( input.getODEsolver() == "BzzMath")
{
    #if ASALI_USE_BZZ == 1
    ode->resize(input.getModelName()[j]);

    BzzVector yFbzz;
    BzzVector x0bzz;
    ChangeDimensions( ode->NumberOfEquations(), &yFbzz, true);
    ChangeDimensions( ode->NumberOfEquations(), &x0bzz, true);
    ChangeDimensions( ode->NumberOfEquations(), &yF,    true);
    ChangeDimensions( ode->NumberOfEquations(), &x0,    true);

    if (input.getModelName()[j] == "honeyComb")
    {
        unsigned int counter = 1;
        for (unsigned int a=0;a<input.axialPoints();a++)
        {
            for (unsigned int r=0;r<input.radialPoints();r++)
            {
                for(unsigned int s=1;s<=input.numberOfSpecies();s++)
                    x0bzz[counter++] = input.getMassFraction()[s-1];
                for(unsigned int s=1;s<=input.numberOfSpecies();s++)
                    x0bzz[counter++] = input.getMassFraction()[s-1];
                x0bzz[counter++] = input.getCoolantTemperature();
                x0bzz[counter++] = input.getCoolantTemperature();
            }
        }
    }
    else if (input.getModelName()[j] == "packedBed")
    {
        unsigned int counter = 1;
        for (unsigned int a=0;a<input.axialPoints();a++)
        {
            for (unsigned int r=0;r<input.radialPoints();r++)
            {
                for(unsigned int s=1;s<=input.numberOfSpecies();s++)
                    x0bzz[counter++] = input.getMassFraction()[s-1];
                x0bzz[counter++] = input.getCoolantTemperature();
            }
        }
    }
    else if (input.getModelName()[j] == "microBed")
    {
        unsigned int counter = 1;
        for (unsigned int a=0;a<input.axialPoints();a++)
        {
            for (unsigned int r=0;r<input.radialPoints();r++)
            {
                for(unsigned int s=1;s<=input.numberOfSpecies();s++)
                    x0bzz[counter++] = input.getMassFraction()[s-1];
                x0bzz[counter++] = input.getCoolantTemperature();
                x0bzz[counter++] = input.getCoolantTemperature();
            }
        }
    }

    BzzOdeSparseStiffObject o;
    o(x0bzz,0.,ode,ode->BlockDimension());

    {
        BzzVector yMax(ode->NumberOfEquations());
        if (input.getModelName()[j] == "honeyComb")
        {
            unsigned int counter = 1;
            for (unsigned int a=0;a<input.axialPoints();a++)
            {
                for (unsigned int r=0;r<input.radialPoints();r++)
                {
                    for(unsigned int s=1;s<=input.numberOfSpecies();s++)
                        yMax[counter++] = 1.;
                    for(unsigned int s=1;s<=input.numberOfSpecies();s++)
                        yMax[counter++] = 1.;
                    yMax[counter++] = 5000.;
                    yMax[counter++] = 5000.;
                }
            }
        }
        else if (input.getModelName()[j] == "packedBed")
        {
            unsigned int counter = 1;
            for (unsigned int a=0;a<input.axialPoints();a++)
            {
                for (unsigned int r=0;r<input.radialPoints();r++)
                {
                    for(unsigned int s=1;s<=input.numberOfSpecies();s++)
                        yMax[counter++] = 1.;
                    yMax[counter++] = 5000.;
                }
            }
        }
        else if (input.getModelName()[j] == "microBed")
        {
            unsigned int counter = 1;
            for (unsigned int a=0;a<input.axialPoints();a++)
            {
                for (unsigned int r=0;r<input.radialPoints();r++)
                {
                    for(unsigned int s=1;s<=input.numberOfSpecies();s++)
                        yMax[counter++] = 1.;
                    yMax[counter++] = 5000.;
                    yMax[counter++] = 5000.;
                }
            }
        }
        BzzVector yMin(ode->NumberOfEquations());
        if (input.getModelName()[j] == "honeyComb")
        {
            unsigned int counter = 1;
            for (unsigned int a=0;a<input.axialPoints();a++)
            {
                for (unsigned int r=0;r<input.radialPoints();r++)
                {
                    for(unsigned int s=1;s<=input.numberOfSpecies();s++)
                        yMin[counter++] = 0.;
                    for(unsigned int s=1;s<=input.numberOfSpecies();s++)
                        yMin[counter++] = 0.;
                    yMin[counter++] = 0.;
                    yMin[counter++] = 0.;
                }
            }
        }
        else if (input.getModelName()[j] == "packedBed")
        {
            unsigned int counter = 1;
            for (unsigned int a=0;a<input.axialPoints();a++)
            {
                for (unsigned int r=0;r<input.radialPoints();r++)
                {
                    for(unsigned int s=1;s<=input.numberOfSpecies();s++)
                        yMin[counter++] = 0.;
                    yMin[counter++] = 0.;
                }
            }
        }
        else if (input.getModelName()[j] == "microBed")
        {
            unsigned int counter = 1;
            for (unsigned int a=0;a<input.axialPoints();a++)
            {
                for (unsigned int r=0;r<input.radialPoints();r++)
                {
                    for(unsigned int s=1;s<=input.numberOfSpecies();s++)
                        yMin[counter++] = 0.;
                    yMin[counter++] = 0.;
                    yMin[counter++] = 0.;
                }
            }
        }
        o.SetMinimumConstraints(&yMin);
        o.SetMaximumConstraints(&yMax);
    }

    o.SetTolAbs(1.e-12);
    o.SetTolRel(1.e-07);

    ode->start(input.getModelName()[j]);
    yFbzz = o(1e10);

    FromBzzToOS(yFbzz,yF);

    for (unsigned int i=1;i<=ode->NumberOfEquations();i++)
    {
        x0[i] = yF[i];
        yF[i] = 0.;
    }
    
    ode->end(input.getModelName()[j]);
    #endif
}
else if ( input.getODEsolver() == "Sundials" )
{
    #if ASALI_USE_SUNDIALS == 1
    ode->resize(input.getModelName()[j]);

    ChangeDimensions( ode->NumberOfEquations(), &yF, true);
    ChangeDimensions( ode->NumberOfEquations(), &x0, true);
    if (input.getModelName()[j] == "honeyComb")
    {
        unsigned int counter = 1;
        for (unsigned int a=0;a<input.axialPoints();a++)
        {
            for (unsigned int r=0;r<input.radialPoints();r++)
            {
                for(unsigned int s=1;s<=input.numberOfSpecies();s++)
                    x0[counter++] = input.getMassFraction()[s-1];
                for(unsigned int s=1;s<=input.numberOfSpecies();s++)
                    x0[counter++] = input.getMassFraction()[s-1];
                x0[counter++] = input.getCoolantTemperature();
                x0[counter++] = input.getCoolantTemperature();
            }
        }
    }
    else if (input.getModelName()[j] == "packedBed")
    {
        unsigned int counter = 1;
        for (unsigned int a=0;a<input.axialPoints();a++)
        {
            for (unsigned int r=0;r<input.radialPoints();r++)
            {
                for(unsigned int s=1;s<=input.numberOfSpecies();s++)
                    x0bzz[counter++] = input.getMassFraction()[s-1];
                x0bzz[counter++] = input.getCoolantTemperature();
            }
        }
    }
    else if (input.getModelName()[j] == "microBed")
    {
        unsigned int counter = 1;
        for (unsigned int a=0;a<input.axialPoints();a++)
        {
            for (unsigned int r=0;r<input.radialPoints();r++)
            {
                for(unsigned int s=1;s<=input.numberOfSpecies();s++)
                    x0bzz[counter++] = input.getMassFraction()[s-1];
                x0bzz[counter++] = input.getCoolantTemperature();
                x0bzz[counter++] = input.getCoolantTemperature();
            }
        }
    }

    OpenSMOKE::ODESystem_CVODE_Template *odeSystemObject;
    odeSystemObject = OpenSMOKE::ODESystem_CVODE_Template::GetInstance();
    odeSystemObject->SetOdeSystem(ode);

    OpenSMOKE::OpenSMOKE_CVODE_Sundials<OpenSMOKE::ODESystem_CVODE_Template>   o(odeSystemObject);
    o.SetDimensions(ode->NumberOfEquations());
    o.SetAbsoluteTolerance(1.e-12);
    o.SetRelativeTolerance(1.e-07);
    o.SetMaximumNumberOfSteps(10000);
    o.SetAnalyticalJacobian(false);
    o.SetInitialValues(0., x0.GetHandle());

    ode->start(input.getModelName()[j]);
    double tf = 1e10;
    o.Solve(tf);
    o.Solution(yF.GetHandle());
    
    for (unsigned int i=1;i<=ode->NumberOfEquations();i++)
    {
        x0[i] = yF[i];
        yF[i] = 0.;
    }
    
    ode->end(input.getModelName()[j]);
    #endif
}
