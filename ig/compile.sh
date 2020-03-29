#!/bin/bash
mkdir -q native/build
cd native/build
cmake .. -DCMAKE_BUILD_TYPE=RELEASE
make -j
cd ../../

javac -classpath ext/commons-io-2.6.jar -classpath ext/commons-cli-2.6.jar  -sourcepath ./ infoasys/cli/pangenes/Pangenes.java && jar cvf ig.jar infoasys/
rm $(find -name '*.class')


