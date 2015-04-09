source mybashrc
export ASALI_SUPPORT='-DASALI_USE_SUNDIALS=1 -DASALI_USE_BZZ=1'
mkdir exe
rm -f exe/heterogeneousPFR_with_bzz.sh
rm -f exe/idealPFR_heatTransfer_with_bzz.sh

cd src
g++ -O3 -m64 -Wno-write-strings -fpermissive heterogeneousPFR/main.C $SUPPORT $ASALI_SUPPORT -I$MKL_INCLUDE -IheterogeneousPFR -IheterogeneousPFR/utilities -IheterogeneousPFR/ode -IheterogeneousPFR/bvp -I. -I$BZZ_INCLUDE -I$EIGEN -I$BOOST/include/ -I$OPENSMOKE/ -I$RAPIDXML -L$BZZ_LIBS_PATH -L$BOOST/lib/ $BZZ_LIBS -lgfortran -L$MKL_LIBS_PATH $MKL_LIBS -Wl,--no-as-needed -ldl -lpthread -lm -lboost_program_options -lboost_system -lboost_filesystem -lboost_regex -lboost_date_time -o ../exe/heterogeneousPFR.sh

g++ -O3 -m64 -Wno-write-strings -fpermissive idealPFR_heatTransfer/main.C $SUPPORT $ASALI_SUPPORT -I$MKL_INCLUDE -IidealPFR_heatTransfer -IidealPFR_heatTransfer/utilities -IidealPFR_heatTransfer/ode -I. -I$BZZ_INCLUDE -I$EIGEN -I$BOOST/include/ -I$OPENSMOKE/ -I$RAPIDXML -L$BZZ_LIBS_PATH -L$BOOST/lib/ $BZZ_LIBS -lgfortran -L$MKL_LIBS_PATH $MKL_LIBS -Wl,--no-as-needed -ldl -lpthread -lm -lboost_program_options -lboost_system -lboost_filesystem -lboost_regex -lboost_date_time -o ../exe/idealPFR_heatTransfer.sh
cd ..
