#! /bin/sh

tar -xvf sdsl-lite.tar.gz
cd sdsl-lite
./install.sh "$(pwd)"/libsdsl
mv libsdsl/ ..

cd ../

tar -xvf pSAscan-0.1.0.tar.gz
cd pSAscan-0.1.0
cd src
make
cd ../../

tar -xjvf LCPscan-0.1.0.tar.bz2
cd LCPscan-0.1.0
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE="Release" ..
make 
cd ../../

mkdir tmpdata
cd tmpdata
mkdir data01
mkdir data02
cd ../
cwd=$(pwd)
echo -e  "disk="$cwd"/tmpdata/data01/stxxl,200000,syscall_unlink\ndisk="$cwd"/tmpdata/data02/stxxl,200000,syscall_unlink\n">./.stxxl
