#!/bin/bash

cd src/

touch patch.h
touch results.txt

cd ..
make hqc-256-3-kat
cd bin/
echo 'kat execution'
pwd
./hqc-256-3-kat

cp PQCkemKAT_8937.rsp PQCkemKAT_8937SsPatch.rsp
cd ..
make hqc-256-3
cd bin/

echo 'Perf execution'
pwd
echo 'Perf execution without PATCH'

echo > ../results.txt
echo 'Origin, without patch' >> ../results.txt
echo >> ../results.txt

./hqc-256-3 >> ../results.txt

cd ..
cd src/

echo
echo 'PATCH fast convolution'
echo

echo >> ../results.txt
echo 'PATCH fast convolution' >> ../results.txt
echo >> ../results.txt

echo '#define PATCH'  > patch.h

cat patch.h  >> ../results.txt

cd ..
pwd
make hqc-256-3-kat
cd bin/

echo 'kat execution'
pwd
./hqc-256-3-kat

cp PQCkemKAT_8937.rsp PQCkemKAT_8937P1.rsp
echo 'diff PQCs' >> ../results.txt
OUT=$(diff PQCkemKAT_8937P1.rsp PQCkemKAT_8937SsPatch.rsp)

#echo 'OUT' $OUT
if [ $OUT -eq '']; then
	echo 'diff OK'
	echo 'diff OK' >> ../results.txt
else
	echo 'diff non OK !!!'
	echo 'diff non OK !!!' >> ../results.txt
fi

cd ..
make hqc-256-3
cd bin/

echo 'Perf execution with PATCH'
pwd
./hqc-256-3 >> ../results.txt

echo >> ../results.txt
echo >> ../results.txt



echo
echo 'PATCH syndrome_gen2'
echo

echo >> ../results.txt
echo 'PATCH syndrome_gen2' >> ../results.txt
echo >> ../results.txt

cd ..
cd src/
echo '#define PATCH2'  > patch.h

cat patch.h  >> ../results.txt

cd ..
pwd
make hqc-256-3-kat
cd bin/

echo 'kat execution'
pwd
./hqc-256-3-kat

cp PQCkemKAT_8937.rsp PQCkemKAT_8937P2.rsp
echo 'diff PQCs' >> ../results.txt
OUT=$(diff PQCkemKAT_8937P2.rsp PQCkemKAT_8937SsPatch.rsp)

#echo 'OUT' $OUT
if [ $OUT -eq '']; then
	echo 'diff OK'
	echo 'diff OK' >> ../results.txt
else
	echo 'diff non OK !!!'
	echo 'diff non OK !!!' >> ../results.txt
fi

cd ..
make hqc-256-3
cd bin/

echo 'Perf execution with PATCH'
pwd
./hqc-256-3 >> ../results.txt

echo >> ../results.txt
echo >> ../results.txt


echo
echo 'PATCH cumulation'
echo

echo >> ../results.txt
echo 'PATCH cumulation' >> ../results.txt
echo >> ../results.txt

cd ..
cd src/
echo '#define PATCH'  > patch.h
echo '#define PATCH2'  >> patch.h

cat patch.h  >> ../results.txt


cd ..
pwd
make hqc-256-3-kat
cd bin/

echo 'kat execution'
pwd
./hqc-256-3-kat

cp PQCkemKAT_8937.rsp PQCkemKAT_8937P1P2.rsp
echo 'diff PQCs' >> ../results.txt
OUT=$(diff PQCkemKAT_8937P1P2.rsp PQCkemKAT_8937SsPatch.rsp)

#echo 'OUT' $OUT
if [ $OUT -eq '']; then
	echo 'diff OK'
	echo 'diff OK' >> ../results.txt
else
	echo 'diff non OK !!!'
	echo 'diff non OK !!!' >> ../results.txt
fi

cd ..
make hqc-256-3
cd bin/

echo 'Perf execution with PATCH'
pwd
./hqc-256-3 >> ../results.txt

echo >> ../results.txt
echo >> ../results.txt

