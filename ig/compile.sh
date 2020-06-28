#!/bin/bash
mkdir -p native/build
cd native/build
cmake .. -DCMAKE_BUILD_TYPE=RELEASE
make -j
cd ../../

javac -classpath ext/commons-io-2.6.jar -classpath ext/commons-cli-1.4.jar  -sourcepath ./ infoasys/cli/pangenes/Pangenes.java && jar cvfm ig.jar META-INF/MANIFEST.MF infoasys/ ext/ native/build/*native*
rm $(find -name '*.class')
