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

if ( input.getDaeSolver() == "BzzMath")
{
	#if ASALI_USE_BZZ == 1

	BzzVector yFbzz;
	BzzVector x0bzz;
	ChangeDimensions( dae->NumberOfEquations(), &yFbzz, true);
	ChangeDimensions( dae->NumberOfEquations(), &x0bzz, true);
	ChangeDimensions( dae->NumberOfEquations(), &yF, true);

	FromOSToBzz(x0,x0bzz);

	BzzVectorInt iAO;
	dae->AlgebraicEquations(iAO);

	BzzDaeObject o;
	o(x0bzz,0.,iAO,dae);

	if ( input.getConstraints() == true)
	{
		BzzVector yMax(dae->NumberOfEquations());
		{
			unsigned counter = 1;
			{
				for (unsigned int j=1;j<=NC;j++)
					yMax[counter++] = 1.;
				for (unsigned int j=1;j<=NC;j++)
					yMax[counter++] = 1.;
				for (unsigned int j=1;j<=SURF_NC;j++)
					yMax[counter++] = 1.;
			}
		}
		BzzVector yMin(dae->NumberOfEquations());
		{
			unsigned counter = 1;
			{
				for (unsigned int j=1;j<=NC;j++)
					yMin[counter++] = 0.;
				for (unsigned int j=1;j<=NC;j++)
					yMin[counter++] = 0.;
				for (unsigned int j=1;j<=SURF_NC;j++)
					yMin[counter++] = 0.;
			}
		}
		o.SetMinimumConstraints(&yMin);
		o.SetMaximumConstraints(&yMax);
	}

	{
		BzzVector tolAbs(dae->NumberOfEquations());
		{
			unsigned counter = 1;
			{
				for (unsigned int j=1;j<=NC;j++)
					tolAbs[counter++] = input.getSpecieAbsTol();
				for (unsigned int j=1;j<=NC;j++)
					tolAbs[counter++] = input.getSpecieAbsTol();
				for (unsigned int j=1;j<=SURF_NC;j++)
					tolAbs[counter++] = input.getSpecieAbsTol();

			}
		}
		BzzVector tolRel(dae->NumberOfEquations());
		{
			unsigned counter = 1;
			{
				for (unsigned int j=1;j<=NC;j++)
					tolRel[counter++] = input.getSpecieRelTol();
				for (unsigned int j=1;j<=NC;j++)
					tolRel[counter++] = input.getSpecieRelTol();
				for (unsigned int j=1;j<=SURF_NC;j++)
					tolRel[counter++] = input.getSpecieRelTol();
			}
		}

		o.SetTolAbs(tolAbs);
		o.SetTolRel(tolRel);
	}

	dae->start();
	yFbzz = o(input.getReactorLength());

	FromBzzToOS(yFbzz,yF);

	dae->end();
	#endif
}
else if ( input.getDaeSolver() == "Sundials" )
{
	#if ASALI_USE_SUNDIALS == 1

	OpenSMOKE::OpenSMOKEVectorDouble residuals(dae->NumberOfEquations());
	OpenSMOKE::OpenSMOKEVectorBool   algebraic(dae->NumberOfEquations());

	dae->AlgebraicEquations(algebraic);
	ChangeDimensions( dae->NumberOfEquations(), &yF, true);

	OpenSMOKE::DAESystem_IDA_Template *daeSystemObject;
	daeSystemObject = OpenSMOKE::DAESystem_IDA_Template::GetInstance(dae->NumberOfEquations(),algebraic.GetHandle());
	daeSystemObject->SetDaeSystem(dae);

	OpenSMOKE::OpenSMOKE_IDA_Sundials<OpenSMOKE::DAESystem_IDA_Template>   o(daeSystemObject);
	o.SetDimensions(dae->NumberOfEquations());
	o.SetAbsoluteTolerance(input.getSpecieAbsTol());
	o.SetRelativeTolerance(input.getSpecieRelTol());
	o.SetMaximumNumberOfSteps(100000);
	o.SetAnalyticalJacobian(false);
	o.SetAlgebraicEquations(algebraic.GetHandle());
	o.SetInitialValues(0., x0.GetHandle(), residuals.GetHandle());

	dae->start();
	o.Solve(input.getIntegrationTime());
	o.Solution(yF.GetHandle(), residuals.GetHandle());

	dae->end();
	#endif
}
