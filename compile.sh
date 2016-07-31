reset
echo '/*##############################################################################################'
echo '#                                                                                              #'
echo '#     #############       #############       #############       ####                ####     #'
echo '#    #             #     #             #     #             #     #    #              #    #    #'
echo '#    #    #####    #     #    #########      #    #####    #     #    #              #    #    #'
echo '#    #    #   #    #     #    #              #    #   #    #     #    #              #    #    #'
echo '#    #    #####    #     #    #              #    #####    #     #    #              #    #    #'
echo '#    #             #     #    #########      #             #     #    #              #    #    #'
echo '#    #             #     #             #     #             #     #    #              #    #    #'
echo '#    #    #####    #      #########    #     #    #####    #     #    #              #    #    #'
echo '#    #    #   #    #              #    #     #    #   #    #     #    #              #    #    #'
echo '#    #    #   #    #      #########    #     #    #   #    #     #    #########      #    #    #'
echo '#    #    #   #    #     #             #     #    #   #    #     #             #     #    #    #'
echo '#     ####     ####       #############       ####     ####       #############       ####     #'
echo '#                                                                                              #'
echo '#   Department of Energy                                                                       #'
echo '#   Politecnico di Milano                                                                      #'
echo '#   Author: Stefano Rebughini <stefano.rebughini@polimi.it>                                    #'
echo '#                                                                                              #'
echo '################################################################################################'
echo '#                                                                                              #'
echo '#   License                                                                                    #'
echo '#                                                                                              #'
echo '#   This file is part of ASALI.                                                                #'
echo '#                                                                                              #'
echo '#   ASALI is free software: you can redistribute it and/or modify                              #'
echo '#   it under the terms of the GNU General Public License as published by                       #'
echo '#   the Free Software Foundation, either version 3 of the License, or                          #'
echo '#   (at your option) any later version.                                                        #'
echo '#                                                                                              #'
echo '#   ASALI is distributed in the hope that it will be useful,                                   #'
echo '#   but WITHOUT ANY WARRANTY; without even the implied warranty of                             #'
echo '#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                              #'
echo '#   GNU General Public License for more details.                                               #'
echo '#                                                                                              #'
echo '#   You should have received a copy of the GNU General Public License                          #'
echo '#   along with ASALI. If not, see <http://www.gnu.org/licenses/>.                              #'
echo '#                                                                                              #'
echo '##############################################################################################*/'
echo ' '
echo ' '
if [ $# != 1 ]
then
	echo "Available compiling options: || complete || mkl || bzz || basic ||"
	echo ' '
	echo "Example: ./compile.sh basic"
	echo ' '
	echo ' '	
	exit 1
fi

v=$1
name=mybashrc.$v

echo "Compiling $v version"
echo ' '
echo ' '

source $name

export BOOST_LIBS='-Wl,--start-group -Wl,-Bstatic -lboost_program_options -lboost_system -lboost_filesystem -lboost_regex -lboost_date_time -Wl,--end-group -Wl,-Bdynamic'

mkdir exe

#2D models
rm -f exe/2Ds-$v.sh
g++ -Ofast -m64 -Wno-write-strings -fpermissive src/2Ds/main.C $SUPPORT $ASALI_SUPPORT -I$MKL_INCLUDE -I2Ds -Isrc/2Ds/ode -Isrc/2Ds/bvp -Isrc/2Ds/utilities -I$BZZ_INCLUDE -I$EIGEN -I$BOOST/include/ -I$OPENSMOKE/ -I$RAPIDXML -L$BZZ_LIBS_PATH -L$BOOST/lib/ $BZZ_LIBS -lgfortran -L$MKL_LIBS_PATH $MKL_LIBS -Wl,--no-as-needed -ldl -lpthread -lm $BOOST_LIBS -o exe/2Ds-$v.sh

#heat transfer pfr
rm -f exe/heatTransferPFR-$v.sh
g++ -Ofast -m64 -Wno-write-strings -fpermissive src/heatTransferPFR/main.C $SUPPORT $ASALI_SUPPORT -I$MKL_INCLUDE -IheatTransferPFR -Isrc/heatTransferPFR/utilities -Isrc/heatTransferPFR/ode -I. -I$BZZ_INCLUDE -I$EIGEN -I$BOOST/include/ -I$OPENSMOKE/ -I$RAPIDXML -L$BZZ_LIBS_PATH -L$BOOST/lib/ $BZZ_LIBS -lgfortran -L$MKL_LIBS_PATH $MKL_LIBS -Wl,--no-as-needed -ldl -lpthread -lm $BOOST_LIBS -o exe/heatTransferPFR-$v.sh

#heterogeneous pfr
rm -f exe/heterogeneousPFR-$v.sh
g++ -Ofast -m64 -Wno-write-strings -fpermissive src/heterogeneousPFR/main.C $SUPPORT $ASALI_SUPPORT -I$MKL_INCLUDE -IheterogeneousPFR -Isrc/heterogeneousPFR/utilities -Isrc/heterogeneousPFR/ode -Isrc/heterogeneousPFR/bvp -Isrc/heterogeneousPFR/post -I. -I$BZZ_INCLUDE -I$EIGEN -I$BOOST/include/ -I$OPENSMOKE/ -I$RAPIDXML -L$BZZ_LIBS_PATH -L$BOOST/lib/ $BZZ_LIBS -lgfortran -L$MKL_LIBS_PATH $MKL_LIBS -Wl,--no-as-needed -ldl -lpthread -lm $BOOST_LIBS -o exe/heterogeneousPFR-$v.sh

#ideal batch
rm -f exe/idealBATCH-$v.sh
g++ -Ofast -m64 -Wno-write-strings -fpermissive src/idealBATCH/main.C $SUPPORT $ASALI_SUPPORT -I$MKL_INCLUDE -IidealBATCH -Isrc/idealBATCH/utilities -Isrc/idealBATCH/ode -Isrc/idealBATCH/ic -I. -I. -I$BZZ_INCLUDE -I$EIGEN -I$BOOST/include/ -I$OPENSMOKE/ -I$RAPIDXML -L$BZZ_LIBS_PATH -L$BOOST/lib/ $BZZ_LIBS -lgfortran -L$MKL_LIBS_PATH $MKL_LIBS -Wl,--no-as-needed -ldl -lpthread -lm $BOOST_LIBS -o exe/idealBATCH-$v.sh

#ideal pfr
rm -f exe/idealPFR-$v.sh
g++ -Ofast -m64 -Wno-write-strings -fpermissive src/idealPFR/main.C $SUPPORT $ASALI_SUPPORT -I$MKL_INCLUDE -IidealPFR -Isrc/idealPFR/utilities -Isrc/idealPFR/ode -Isrc/idealPFR/dae -I. -I$BZZ_INCLUDE -I$EIGEN -I$BOOST/include/ -I$OPENSMOKE/ -I$RAPIDXML -L$BZZ_LIBS_PATH -L$BOOST/lib/ $BZZ_LIBS -lgfortran -L$MKL_LIBS_PATH $MKL_LIBS -Wl,--no-as-needed -ldl -lpthread -lm $BOOST_LIBS -o exe/idealPFR-$v.sh

#mass transfer pfr
rm -f exe/massTransferPFR-$v.sh
g++ -Ofast -m64 -Wno-write-strings -fpermissive src/massTransferPFR/main.C $SUPPORT $ASALI_SUPPORT -I$MKL_INCLUDE -Isrc/massTransferPFR/. -Isrc/massTransferPFR/utilities -Isrc/massTransferPFR/dae -Isrc/massTransferPFR/ode  -I$BZZ_INCLUDE -I$EIGEN -I$BOOST/include/ -I$OPENSMOKE/ -I$RAPIDXML -L$BZZ_LIBS_PATH -L$BOOST/lib/ $BZZ_LIBS -lgfortran -L$MKL_LIBS_PATH $MKL_LIBS -Wl,--no-as-needed -ldl -lpthread -lm $BOOST_LIBS -o exe/massTransferPFR-$v.sh

#packed bed
rm -f exe/packedBed-$v.sh
g++ -Ofast -m64 -Wno-write-strings -fpermissive src/packedBed/main.C $SUPPORT $ASALI_SUPPORT -I$MKL_INCLUDE -IpackedBed -Isrc/packedBed/ode -Isrc/packedBed/dae -Isrc/packedBed/utilities -I$BZZ_INCLUDE -I$EIGEN -I$BOOST/include/ -I$OPENSMOKE/ -I$RAPIDXML -L$BZZ_LIBS_PATH -L$BOOST/lib/ $BZZ_LIBS -lgfortran -L$MKL_LIBS_PATH $MKL_LIBS -Wl,--no-as-needed -ldl -lpthread -lm $BOOST_LIBS -o exe/packedBed-$v.sh
