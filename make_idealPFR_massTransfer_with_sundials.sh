source mybashrc
export ASALI_SUPPORT='-DASALI_USE_SUNDIALS=1 -DASALI_USE_BZZ=0'
mkdir exe
rm -f exe/idealPFR_massTransfer_with_sundials.sh

cd src
g++ -O3 -m64 -Wno-write-strings -fpermissive idealPFR_massTransfer/main.C $SUPPORT $ASALI_SUPPORT -I$MKL_INCLUDE -IidealPFR_massTransfer -IidealPFR_massTransfer/utilities -IidealPFR_massTransfer/ode -IidealPFR_massTransfer/dae -I. -I$EIGEN -I$BOOST/include/ -I$OPENSMOKE/ -I$RAPIDXML -I$SUNDIALS_INCLUDE -L$BOOST/lib/ -lgfortran -L$SUNDIALS_LIBS_PATH $SUNDIALS_LIBS -L$MKL_LIBS_PATH $MKL_LIBS -Wl,--no-as-needed -ldl -lpthread -lm -lboost_program_options -lboost_system -lboost_filesystem -lboost_regex -lboost_date_time -o ../exe/idealPFR_massTransfer_with_sundials.sh
cd ..
