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

tar -xvf EM-SparsePhi-0.1.0.tar.gz
cd EM-SparsePhi-0.1.0
cd src
make
cd ../../

mkdir tmpdata
cd tmpdata
mkdir data01
mkdir data02
cd ../
cwd=$(pwd)
echo -e  "disk="$cwd"/tmpdata/data01/stxxl,200000,syscall_unlink\ndisk="$cwd"/tmpdata/data02/stxxl,200000,syscall_unlink\n">./.stxxl
