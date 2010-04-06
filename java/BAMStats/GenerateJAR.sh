#!/bin/bash

picardPath="/stornext/snfs1/next-gen/software/picard-tools/current"

samJarName=`ls $picardPath"/"sam-*.jar`
picardJarName=`ls $picardPath"/"picard-*.jar`

rm -rf BAMStats* *.jar *.class

echo "SAM Jar : "$samJarName
echo "Picard Jar : "$picardJarName

echo "Compiling project"
javac -classpath $samJarName *.java

echo "Generating Manifest file"
manifestFile=`pwd`"/BAMStatsManifest.txt"

echo -e "Class-Path: "$samJarName" "$picardJarName"\nMain-Class: BAMStats.Driver\n" > $manifestFile

echo "Building Jar file"

cd ..
jar cvfm BAMStats.jar $manifestFile BAMStats/*.class
mv BAMStats.jar ./BAMStats
