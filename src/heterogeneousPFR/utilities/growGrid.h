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
	class GROWgrid 
	{
		public:

			GROWgrid(unsigned int& NC, unsigned int& SURF_NC);

			#include "vector.h"

			void resize(const unsigned int NP);

			void setVariables  (const OpenSMOKE::OpenSMOKEVectorDouble  y);
			void setGrid       (const OpenSMOKE::OpenSMOKEVectorDouble z);
			void setAddPoints  (const unsigned int addN);
			void setMaxPoints  (const unsigned int maxN);
			void setTol        (const double tol);

			void calculateErrors();
			void addition();
			void interpolation();
			void recapOnScreen();

			unsigned int					 getInterpolatedPoints()      { return interpolatedNP_;};
			bool							 getExitCondition();
			OpenSMOKE::OpenSMOKEVectorDouble getInterpolatedGrid()        { return interpolatedZ_;};
			OpenSMOKE::OpenSMOKEVectorDouble getInterpolatedVariables()   { return interpolatedVariables_;};

			~GROWgrid(){ };

		private:

			double a_;
			double b_;
			double tol_;

			unsigned int addN_;
			unsigned int addNfixed_;
			unsigned int maxN_;
			unsigned int NE_;
			unsigned int NP_;
			unsigned int interpolatedNP_;

			const unsigned int& NC_;
			const unsigned int& SURF_NC_;

			bool ended;

			OpenSMOKE::OpenSMOKEVectorDouble* xBulk_;
			OpenSMOKE::OpenSMOKEVectorDouble* xWall_;
			OpenSMOKE::OpenSMOKEVectorDouble* teta_;
			OpenSMOKE::OpenSMOKEVectorDouble* errorMatrix_;
			OpenSMOKE::OpenSMOKEVectorDouble* addxBulk_;
			OpenSMOKE::OpenSMOKEVectorDouble* addxWall_;
			OpenSMOKE::OpenSMOKEVectorDouble* addteta_;
			OpenSMOKE::OpenSMOKEVectorDouble* interpolatedxBulk_;
			OpenSMOKE::OpenSMOKEVectorDouble* interpolatedxWall_;
			OpenSMOKE::OpenSMOKEVectorDouble* interpolatedteta_;

			OpenSMOKE::OpenSMOKEVectorDouble Tbulk_;
			OpenSMOKE::OpenSMOKEVectorDouble Twall_;
			OpenSMOKE::OpenSMOKEVectorDouble addTbulk_;
			OpenSMOKE::OpenSMOKEVectorDouble addTwall_;
			OpenSMOKE::OpenSMOKEVectorDouble interpolatedTbulk_;
			OpenSMOKE::OpenSMOKEVectorDouble interpolatedTwall_;
			OpenSMOKE::OpenSMOKEVectorDouble z_;
			OpenSMOKE::OpenSMOKEVectorDouble addZ_;
			OpenSMOKE::OpenSMOKEVectorDouble interpolatedZ_;
			OpenSMOKE::OpenSMOKEVectorDouble interpolatedVariables_;
			OpenSMOKE::OpenSMOKEVectorDouble errors_;
			OpenSMOKE::OpenSMOKEVectorDouble Dz_;

			OpenSMOKE::OpenSMOKEVectorBool errorBool_;

			void weightedError(const OpenSMOKE::OpenSMOKEVectorDouble variable, OpenSMOKE::OpenSMOKEVectorDouble &error);
			void insertion(const OpenSMOKE::OpenSMOKEVectorDouble Old, const OpenSMOKE::OpenSMOKEVectorDouble Add, OpenSMOKE::OpenSMOKEVectorDouble &New);
	};
	
	GROWgrid::GROWgrid(unsigned int& NC, unsigned int& SURF_NC):
	NC_(NC),
	SURF_NC_(SURF_NC)
	{
		a_   = 0.;
		b_   = 0.;
		tol_ = 0.;

		addN_           = 0;
		addNfixed_      = 0;
		maxN_           = 0;
		NE_             = 0;
		NP_             = 0;
		interpolatedNP_ = 0;

		ended = false;
	}

	void GROWgrid::resize(const unsigned int NP)
	{
		NP_ = NP;

		NE_ = NC_ + NC_ + SURF_NC_ + 1 + 1;

		xBulk_         = new OpenSMOKE::OpenSMOKEVectorDouble[NP_];
		xWall_         = new OpenSMOKE::OpenSMOKEVectorDouble[NP_];
		teta_          = new OpenSMOKE::OpenSMOKEVectorDouble[NP_];
		errorMatrix_   = new OpenSMOKE::OpenSMOKEVectorDouble[NP_];

		for (unsigned int i=0;i<NP_;i++)
		{
			ChangeDimensions(NC_,      &xBulk_[i],       true);
			ChangeDimensions(NC_,      &xWall_[i],       true);
			ChangeDimensions(SURF_NC_, &teta_[i],        true);
			ChangeDimensions(NE_,      &errorMatrix_[i], true);
		}
		
		ChangeDimensions(NP_, &Tbulk_,         true);
		ChangeDimensions(NP_, &Twall_,         true);
		ChangeDimensions(NP_, &errors_,        true);
		ChangeDimensions(NP_, &errorBool_,     true);
		ChangeDimensions(NP_, &interpolatedZ_, true);

		a_ = 0.6;
		b_ = 1.;
	}

	void GROWgrid::setVariables(const OpenSMOKE::OpenSMOKEVectorDouble  y)
	{
		unsigned int counter = 1;
		for (unsigned int i=0;i<NP_;i++)
		{
			for (unsigned int j=1;j<=NC_;j++)
				xBulk_[i][j] = y[counter++];
			for (unsigned int j=1;j<=NC_;j++)
				xWall_[i][j] = y[counter++];
			for (unsigned int j=1;j<=SURF_NC_;j++)
				teta_[i][j] = y[counter++];
			Tbulk_[i+1] = y[counter++];
			Twall_[i+1] = y[counter++];
		}
	}
	
	void GROWgrid::setGrid(const OpenSMOKE::OpenSMOKEVectorDouble z)
	{
		ChangeDimensions(z.Size(), &z_, true);
		for (unsigned int k=1;k<=z_.Size();k++)
			z_[k] = z[k];
		
		ChangeDimensions((NP_- 1), &Dz_, true);
		for (unsigned int k=1;k<=Dz_.Size();k++)
			Dz_[k] = z_[k+1] - z_[k];
	}

	void GROWgrid::setAddPoints(const unsigned int addN)
	{
		addNfixed_ = addN;
	}

	void GROWgrid::setMaxPoints(const unsigned int maxN)
	{
		maxN_ = maxN;
	}

	void GROWgrid::setTol(const double tol)
	{
		tol_ = tol;
	}

	void GROWgrid::insertion(const OpenSMOKE::OpenSMOKEVectorDouble Old, const OpenSMOKE::OpenSMOKEVectorDouble Add, OpenSMOKE::OpenSMOKEVectorDouble &New)
	{
		ChangeDimensions(interpolatedNP_, &New, true);
		unsigned int counter = 0;
		for(unsigned int k=1;k<=NP_;k++)
		{
			New[k+counter] = Old[k];
			if ( errorBool_[k] == true)
			{
				counter++;
				New[k+counter] = Add[counter];
			}
		}
	}

	void GROWgrid::weightedError(const OpenSMOKE::OpenSMOKEVectorDouble variable, OpenSMOKE::OpenSMOKEVectorDouble &error)
	{
		OpenSMOKE::OpenSMOKEVectorDouble errorFirst(NP_);
		for (unsigned int k=1;k<NP_;k++)
			errorFirst[k]=fabs((variable[k]-variable[k-1])/Dz_[k-1]);
		errorFirst[0] = errorFirst[1];
		
		OpenSMOKE::OpenSMOKEVectorDouble errorSecond(NP_ - 1);
		errorSecond[0] = 0.;
		for (unsigned int k=1;k<(NP_-1);k++)
			errorSecond[k]=fabs((variable[k+1]*Dz_[k-1] + variable[k-1]*Dz_[k]-variable[k]*(Dz_[k]+Dz_[k-1]))/(Dz_[k]*Dz_[k]*Dz_[k-1]));
		
		for (unsigned int k=0;k<(NP_-1);k++)
			error[k] = a_*errorFirst[k] + b_*errorSecond[k];
		error[NP_-1]=0.;
	}

	void GROWgrid::calculateErrors()
	{
		for (unsigned int j=1;j<=NC_;j++)
		{
			OpenSMOKE::OpenSMOKEVectorDouble variable(NP_ );
			for ( unsigned int k=0;k<NP_;k++)
				variable[k+1] = std::max(0.,xBulk_[k][j]);

			weightedError(variable, errors_);

			for ( unsigned int k=0;k<NP_;k++)
				errorMatrix_[k][j] = errors_[k+1];
		}

		for (unsigned int j=1;j<=NC_;j++)
		{
			OpenSMOKE::OpenSMOKEVectorDouble variable(NP_ );
			for ( unsigned int k=0;k<NP_;k++)
				variable[k+1] = std::max(0.,xWall_[k][j]);

			weightedError(variable, errors_);

			for ( unsigned int k=0;k<NP_;k++)
				errorMatrix_[k][j+NC_] = errors_[k+1];
		}

		for (unsigned int j=1;j<=SURF_NC_;j++)
		{
			OpenSMOKE::OpenSMOKEVectorDouble variable(NP_ );
			for ( unsigned int k=0;k<NP_;k++)
				variable[k+1] = std::max(0.,teta_[k][j]);

			weightedError(variable, errors_);

			for ( unsigned int k=0;k<NP_;k++)
				errorMatrix_[k][j+NC_+NC_] = errors_[k+1];
		}

		weightedError(Tbulk_, errors_);
		for ( unsigned int k=0;k<NP_;k++)
			errorMatrix_[k][NC_+NC_+SURF_NC_+1] = errors_[k+1];

		weightedError(Twall_, errors_);
		for ( unsigned int k=0;k<NP_;k++)
			errorMatrix_[k][NC_+NC_+SURF_NC_+1+1] = errors_[k+1];

		for ( unsigned int k=0;k<NP_;k++)
			errors_[k+1] = errorMatrix_[k].Max();

		for (unsigned int k=1;k<=errorBool_.Size();k++)
			errorBool_[k] = false;

		if (addNfixed_ >= NP_)
		{
			for (unsigned int k=1;k<=NP_;k++)
				errorBool_[k] = true;
			
			errorBool_[NP_] = false;
		}
		else
		{
			OpenSMOKE::OpenSMOKEVectorDouble fakeError(NP_);

			for (unsigned int k=1;k<=NP_;k++)
				fakeError[k] = fabs(errors_[k]);
			
			unsigned int counter = 1;
			for (unsigned int j=1;j<=addNfixed_;j++)
			{
				double max = fakeError.Max();
				for (unsigned int k=1;k<=NP_;k++)
				{
					if (fakeError[k] == max)
					{
						fakeError[k] = -1.;
						errorBool_[k] = true;
					}
				}
			}

			for (unsigned int k=1;k<=NP_;k++)
			{
				if (errors_[k] <= tol_)
					errorBool_[k] = false;
			}

			errorBool_[NP_] = false;
		}

		addN_ = 0;
		for (unsigned int k=1;k<=NP_;k++)
		{
			if (errorBool_[k] == true)
			{
				addN_++;
			}
		}

		interpolatedNP_ = addN_ + NP_;

	}

	void GROWgrid::addition()
	{
		addxBulk_  = new OpenSMOKE::OpenSMOKEVectorDouble[addN_];
		addxWall_  = new OpenSMOKE::OpenSMOKEVectorDouble[addN_];
		addteta_   = new OpenSMOKE::OpenSMOKEVectorDouble[addN_];

		for (unsigned int i=0;i<addN_;i++)
		{
			ChangeDimensions(NC_,      &addxBulk_[i], true);
			ChangeDimensions(NC_,      &addxWall_[i], true);
			ChangeDimensions(SURF_NC_, &addteta_[i],  true);
		}
		ChangeDimensions(addN_, &addTbulk_, true);
		ChangeDimensions(addN_, &addTwall_, true);
		ChangeDimensions(addN_, &addZ_, true);

		unsigned int counter = 0;
		for (unsigned int k=0;k<NP_;k++)
		{
			if (errorBool_[k+1] == true)
			{
				for (unsigned int j=1;j<=NC_;j++)
				{
					addxBulk_[counter][j] = (xBulk_[k][j] + xBulk_[k+1][j])*0.5;
					addxWall_[counter][j] = (xWall_[k][j] + xWall_[k+1][j])*0.5;
				}
				for (unsigned int j=1;j<=SURF_NC_;j++)
				{
					addteta_[counter][j] = (teta_[k][j] + teta_[k+1][j])*0.5;
				}
				addTbulk_[counter+1] = (Tbulk_[k+1] + Tbulk_[k+1+1])*0.5;
				addTwall_[counter+1] = (Twall_[k+1] + Twall_[k+1+1])*0.5;
				addZ_[counter+1]     = (z_[k+1]     + z_[k+1+1])*0.5;
				counter++;
			}
		}
	}

	void GROWgrid::interpolation()
	{
		interpolatedxBulk_ = new OpenSMOKE::OpenSMOKEVectorDouble[interpolatedNP_];
		interpolatedxWall_ = new OpenSMOKE::OpenSMOKEVectorDouble[interpolatedNP_];
		interpolatedteta_  = new OpenSMOKE::OpenSMOKEVectorDouble[interpolatedNP_];

		for (unsigned int i=0;i<interpolatedNP_;i++)
		{
			ChangeDimensions(NC_,      &interpolatedxBulk_[i], true);
			ChangeDimensions(NC_,      &interpolatedxWall_[i], true);
			ChangeDimensions(SURF_NC_, &interpolatedteta_[i],  true);
		}
		
		ChangeDimensions(interpolatedNP_, &interpolatedTbulk_, true);
		ChangeDimensions(interpolatedNP_, &interpolatedTwall_, true);
		ChangeDimensions(interpolatedNP_, &interpolatedZ_, true);

		for (unsigned int j=1;j<=NC_;j++)
		{
			OpenSMOKE::OpenSMOKEVectorDouble old_(NP_);
			OpenSMOKE::OpenSMOKEVectorDouble add_(addN_);
			OpenSMOKE::OpenSMOKEVectorDouble new_(interpolatedNP_);

			for ( unsigned int k=0;k<NP_;k++)
				old_[k+1] = xBulk_[k][j];
			
			for ( unsigned int k=0;k<addN_;k++)
				add_[k+1] = addxBulk_[k][j];

			insertion(old_,add_,new_);

			for ( unsigned int k=0;k<interpolatedNP_;k++)
				interpolatedxBulk_[k][j] = new_[k+1];
		}

		for (unsigned int j=1;j<=NC_;j++)
		{
			OpenSMOKE::OpenSMOKEVectorDouble old_(NP_);
			OpenSMOKE::OpenSMOKEVectorDouble add_(addN_);
			OpenSMOKE::OpenSMOKEVectorDouble new_(interpolatedNP_);

			for ( unsigned int k=0;k<NP_;k++)
				old_[k+1] = xWall_[k][j];
			
			for ( unsigned int k=0;k<addN_;k++)
				add_[k+1] = addxWall_[k][j];

			insertion(old_,add_,new_);

			for ( unsigned int k=0;k<interpolatedNP_;k++)
				interpolatedxWall_[k][j] = new_[k+1];
		}

		for (unsigned int j=1;j<=SURF_NC_;j++)
		{
			OpenSMOKE::OpenSMOKEVectorDouble old_(NP_);
			OpenSMOKE::OpenSMOKEVectorDouble add_(addN_);
			OpenSMOKE::OpenSMOKEVectorDouble new_(interpolatedNP_);

			for ( unsigned int k=0;k<NP_;k++)
				old_[k+1] = teta_[k][j];
			
			for ( unsigned int k=0;k<addN_;k++)
				add_[k+1] = addteta_[k][j];

			insertion(old_,add_,new_);

			for ( unsigned int k=0;k<interpolatedNP_;k++)
				interpolatedteta_[k][j] = new_[k+1];
		}

		insertion(Tbulk_,addTbulk_,interpolatedTbulk_);
		insertion(Twall_,addTwall_,interpolatedTwall_);

		insertion(z_,addZ_,interpolatedZ_);


		ChangeDimensions(NE_*interpolatedNP_, &interpolatedVariables_, true);
		unsigned int counter = 1;
		for (unsigned int i=0;i<interpolatedNP_;i++)
		{
			for (unsigned int j=1;j<=NC_;j++)
				interpolatedVariables_[counter++] = interpolatedxBulk_[i][j]/interpolatedxBulk_[i].SumElements();
			for (unsigned int j=1;j<=NC_;j++)
				interpolatedVariables_[counter++] = interpolatedxWall_[i][j]/interpolatedxWall_[i].SumElements();
			for (unsigned int j=1;j<=SURF_NC_;j++)
				interpolatedVariables_[counter++] = interpolatedteta_[i][j]/interpolatedteta_[i].SumElements();
			interpolatedVariables_[counter++] = interpolatedTbulk_[i+1];
			interpolatedVariables_[counter++] = interpolatedTwall_[i+1];
		}
	}

	bool GROWgrid::getExitCondition()
	{
		ended = false;
		if ( interpolatedNP_ >= maxN_)
		{
			ended = true;
		}
		else if ( interpolatedNP_ == NP_ )
		{
			ended = true;
		}
		return ended;
	}

	void GROWgrid::recapOnScreen()
	{
		delete [] xBulk_;
		delete [] xWall_;
		delete [] teta_;
		delete [] addxBulk_;
		delete [] addxWall_;
		delete [] addteta_;
		delete [] interpolatedxBulk_;
		delete [] interpolatedxWall_;
		delete [] interpolatedteta_;
		delete [] errorMatrix_;

		std::cout << "\n#################################" << std::endl;
		std::cout << "          INTERPOLATION        " << std::endl;
		std::cout << "START number of points:      " << NP_ << std::endl;
		std::cout << "FINAL number of points:      " << interpolatedNP_ << std::endl;
		std::cout << "ADD   number of points:      " << addN_ << std::endl;
		std::cout << "#################################" << std::endl;
	}
}
