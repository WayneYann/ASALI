source mybashrc
export ASALI_SUPPORT='-DASALI_USE_SUNDIALS=1 -DASALI_USE_BZZ=1'
mkdir exe
rm -f exe/heterogeneousPFR_with.sh
rm -f exe/idealPFR_heatTransfer.sh
rm -f exe/idealPFR_massTransfer.sh
rm -f exe/idealPFR.sh

cd src
g++ -O3 -m64 -Wno-write-strings -fpermissive heterogeneousPFR/main.C $SUPPORT $ASALI_SUPPORT -I$MKL_INCLUDE -IheterogeneousPFR -IheterogeneousPFR/utilities -IheterogeneousPFR/ode -IheterogeneousPFR/bvp -I. -I$BZZ_INCLUDE -I$EIGEN -I$BOOST/include/ -I$OPENSMOKE/ -I$RAPIDXML -I$SUNDIALS_INCLUDE -L$SUNDIALS_LIBS_PATH $SUNDIALS_LIBS -L$BZZ_LIBS_PATH -L$BOOST/lib/ $BZZ_LIBS -lgfortran -L$MKL_LIBS_PATH $MKL_LIBS -Wl,--no-as-needed -ldl -lpthread -lm -lboost_program_options -lboost_system -lboost_filesystem -lboost_regex -lboost_date_time -o ../exe/heterogeneousPFR.sh

g++ -O3 -m64 -Wno-write-strings -fpermissive idealPFR_heatTransfer/main.C $SUPPORT $ASALI_SUPPORT -I$MKL_INCLUDE -IidealPFR_heatTransfer -IidealPFR_heatTransfer/utilities -IidealPFR_heatTransfer/ode -I. -I$BZZ_INCLUDE -I$EIGEN -I$BOOST/include/ -I$OPENSMOKE/ -I$RAPIDXML -I$SUNDIALS_INCLUDE -L$SUNDIALS_LIBS_PATH $SUNDIALS_LIBS -L$BZZ_LIBS_PATH -L$BOOST/lib/ $BZZ_LIBS -lgfortran -L$MKL_LIBS_PATH $MKL_LIBS -Wl,--no-as-needed -ldl -lpthread -lm -lboost_program_options -lboost_system -lboost_filesystem -lboost_regex -lboost_date_time -o ../exe/idealPFR_heatTransfer.sh

g++ -O3 -m64 -Wno-write-strings -fpermissive idealPFR/main.C $SUPPORT $ASALI_SUPPORT -I$MKL_INCLUDE -IidealPFR -IidealPFR/utilities -IidealPFR/ode -IidealPFR/dae -I. -I$BZZ_INCLUDE -I$EIGEN -I$BOOST/include/ -I$OPENSMOKE/ -I$RAPIDXML -I$SUNDIALS_INCLUDE -L$SUNDIALS_LIBS_PATH $SUNDIALS_LIBS -L$BZZ_LIBS_PATH -L$BOOST/lib/ $BZZ_LIBS -lgfortran -L$MKL_LIBS_PATH $MKL_LIBS -Wl,--no-as-needed -ldl -lpthread -lm -lboost_program_options -lboost_system -lboost_filesystem -lboost_regex -lboost_date_time -o ../exe/idealPFR.sh

g++ -O3 -m64 -Wno-write-strings -fpermissive idealPFR_massTransfer/main.C $SUPPORT $ASALI_SUPPORT -I$MKL_INCLUDE -IidealPFR_massTransfer -IidealPFR_massTransfer/utilities -IidealPFR_massTransfer/ode -I. -I$BZZ_INCLUDE -I$EIGEN -I$BOOST/include/ -I$OPENSMOKE/ -I$RAPIDXML -I$SUNDIALS_INCLUDE -L$SUNDIALS_LIBS_PATH $SUNDIALS_LIBS -L$BZZ_LIBS_PATH -L$BOOST/lib/ $BZZ_LIBS -lgfortran -L$MKL_LIBS_PATH $MKL_LIBS -Wl,--no-as-needed -ldl -lpthread -lm -lboost_program_options -lboost_system -lboost_filesystem -lboost_regex -lboost_date_time -o ../exe/idealPFR_massTransfer.sh

cd ..
