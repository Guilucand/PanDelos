#!/bin/bash
mkdir native/build
cd native/build
cmake .. -DCMAKE_BUILD_TYPE=RELEASE
make -j
cd ../../
javac -classpath ext/commons-io-2.6.jar  -sourcepath ./ infoasys/cli/pangenes/Pangenes.java & jar cvf ig.jar infoasys/
