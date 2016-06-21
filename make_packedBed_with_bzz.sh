source mybashrc
export ASALI_SUPPORT='-DASALI_USE_SUNDIALS=0 -DASALI_USE_BZZ=1'
mkdir exe
rm -f exe/packedBed_with_bzz.sh

cd src
g++ -O3 -m64 -Wno-write-strings -fpermissive packedBed/main.C $SUPPORT $ASALI_SUPPORT -I$MKL_INCLUDE -IpackedBed -IpackedBed/ode -IpackedBed/dae -IpackedBed/utilities -I$BZZ_INCLUDE -I$EIGEN -I$BOOST/include/ -I$OPENSMOKE/ -I$RAPIDXML -L$BZZ_LIBS_PATH -L$BOOST/lib/ $BZZ_LIBS -lgfortran -L$MKL_LIBS_PATH $MKL_LIBS -Wl,--no-as-needed -ldl -lpthread -lm -lboost_program_options -lboost_system -lboost_filesystem -lboost_regex -lboost_date_time -o ../exe/packedBed_with_bzz.sh
cd ..
