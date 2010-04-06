#!/bin/bash

picardPath="/stornext/snfs1/next-gen/software/picard-tools/current"

samJarName=`ls $picardPath"/"sam-*.jar`
picardJarName=`ls $picardPath"/"picard-*.jar`
outputJarName="raw.bam.reads.validator.jar"
rm -rf *.jar *.class

echo "SAM Jar : "$samJarName
echo "Picard Jar : "$picardJarName

echo "Compiling project"
javac -classpath $samJarName *.java

echo "Generating Manifest file"
manifestFile=`pwd`"/ValidateSOLiDBAMManifest.txt"

echo -e "Class-Path: "$samJarName" "$picardJarName"\nMain-Class: Driver\n" > $manifestFile

echo "Building Jar file"

jar cvfm $outputJarName $manifestFile *.class
