#!/bin/bash

# Abort the script if any command failed to return 0
set -e

# Print the instructios
set -x

rm -f val*

function check_kats {
    # HQC 1 LEVEL 3
    diff hqc-128-1/PQCkemKAT_3165.rsp ../KATs/Constant_Time_HQC/hqc-128-1/hqc_128-1_kat.rsp > /dev/null
    sha1sum hqc-128-1/PQCkemKAT_3165.rsp
    sha1sum ../KATs/Constant_Time_HQC/hqc-128-1/hqc_128-1_kat.rsp
    # HQC 2 LEVEL 2
    diff hqc-192-1/PQCkemKAT_5539.rsp ../KATs/Constant_Time_HQC/hqc-192-1/hqc_192-1_kat.rsp > /dev/null
    sha1sum hqc-192-1/PQCkemKAT_5539.rsp
    sha1sum ../KATs/Constant_Time_HQC/hqc-192-1/hqc_192-1_kat.rsp
    # HQC 2 LEVEL 3
    diff hqc-192-2/PQCkemKAT_5924.rsp ../KATs/Constant_Time_HQC/hqc-192-2/hqc_192-2_kat.rsp > /dev/null
    sha1sum hqc-192-2/PQCkemKAT_5924.rsp
    sha1sum ../KATs/Constant_Time_HQC/hqc-192-2/hqc_192-2_kat.rsp
    # HQC 3 LEVEL 2
    diff hqc-256-1/PQCkemKAT_8029.rsp ../KATs/Constant_Time_HQC/hqc-256-1/hqc_256-1_kat.rsp > /dev/null
    sha1sum hqc-256-1/PQCkemKAT_8029.rsp
    sha1sum ../KATs/Constant_Time_HQC/hqc-256-1/hqc_256-1_kat.rsp
    # HQC 3 LEVEL 3
    diff hqc-256-2/PQCkemKAT_8543.rsp ../KATs/Constant_Time_HQC/hqc-256-2/hqc_256-2_kat.rsp > /dev/null
    sha1sum hqc-256-2/PQCkemKAT_8543.rsp
    sha1sum ../KATs/Constant_Time_HQC/hqc-256-2/hqc_256-2_kat.rsp
    # HQC 3 LEVEL 4
    diff hqc-256-3/PQCkemKAT_8937.rsp ../KATs/Constant_Time_HQC/hqc-256-3/hqc_256-3_kat.rsp > /dev/null
    sha1sum hqc-256-3/PQCkemKAT_8937.rsp
    sha1sum ../KATs/Constant_Time_HQC/hqc-256-3/hqc_256-3_kat.rsp
    # Clean
    rm hqc-128-1/PQCkemKAT*
    rm hqc-192-1/PQCkemKAT*
    rm hqc-192-2/PQCkemKAT*
    rm hqc-256-1/PQCkemKAT*
    rm hqc-256-2/PQCkemKAT*
    rm hqc-256-3/PQCkemKAT*
}

# Compile all
cd hqc-128-1
make clean; 
make hqc-128-1-kat;  > /dev/null
./bin/hqc-128-1-kat
rm -rf ./bin
rm -rf doc/html

make hqc-128-1 RUN=1;  > /dev/null
./bin/hqc-128-1
rm -rf ./bin
rm -rf doc/html

make hqc-128-1 VALGRIND=1;
valgrind --log-file=../val_Basic3 --leak-check=full ./bin/hqc-128-1
date >> ../val_Basic3
rm -rf ./bin

cd ../hqc-192-1
make clean; 
make hqc-192-1-kat;  > /dev/null
./bin/hqc-192-1-kat
rm -rf ./bin
rm -rf doc/html

make hqc-192-1 RUN=1;  > /dev/null
./bin/hqc-192-1
rm -rf ./bin
rm -rf doc/html

make hqc-192-1 VALGRIND=1;
valgrind --log-file=../val_Advanced2 --leak-check=full ./bin/hqc-192-1
date >> ../val_Advanced2
rm -rf ./bin

cd ../hqc-192-2
make clean; 
make hqc-192-2-kat;  > /dev/null
./bin/hqc-192-2-kat
rm -rf ./bin
rm -rf doc/html

make hqc-192-2 RUN=1;  > /dev/null
./bin/hqc-192-2
rm -rf ./bin
rm -rf doc/html

make hqc-192-2 VALGRIND=1;
valgrind --log-file=../val_Advanced3 --leak-check=full ./bin/hqc-192-2
date >> ../val_Advanced3
rm -rf ./bin

cd ../hqc-256-1
make clean; 
make hqc-256-1-kat;  > /dev/null
./bin/hqc-256-1-kat
rm -rf ./bin
rm -rf doc/html

make hqc-256-1 RUN=1;  > /dev/null
./bin/hqc-256-1
rm -rf ./bin
rm -rf doc/html

make hqc-256-1 VALGRIND=1;
valgrind --log-file=../val_Paranoiac2 --leak-check=full ./bin/hqc-256-1
date >> ../val_Paranoiac2
rm -rf ./bin

cd ../hqc-256-2
make clean; 
make hqc-256-2-kat;  > /dev/null
./bin/hqc-256-2-kat
rm -rf ./bin
rm -rf doc/html

make hqc-256-2 RUN=1;  > /dev/null
./bin/hqc-256-2
rm -rf ./bin
rm -rf doc/html

make hqc-256-2 VALGRIND=1;
valgrind --log-file=../val_Paranoiac3 --leak-check=full ./bin/hqc-256-2
date >> ../val_Paranoiac3
rm -rf ./bin

cd ../hqc-256-3
make clean; 
make hqc-256-3-kat;  > /dev/null
./bin/hqc-256-3-kat
rm -rf ./bin
rm -rf doc/html

make hqc-256-3 RUN=1;  > /dev/null
./bin/hqc-256-3
rm -rf ./bin
rm -rf doc/html

make hqc-256-3 VALGRIND=1;
valgrind --log-file=../val_Paranoiac4 --leak-check=full ./bin/hqc-256-3
date >> ../val_Paranoiac4
rm -rf ./bin

cd ..
check_kats
