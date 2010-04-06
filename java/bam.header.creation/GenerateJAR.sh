#!/bin/bash

picardPath="/stornext/snfs1/next-gen/software/picard-tools/current"

samJarName=`ls $picardPath"/"sam-*.jar`
picardJarName=`ls $picardPath"/"picard-*.jar`

echo "SAM Jar : "$samJarName
echo "Picard Jar : "$picardJarName

echo "Compiling project"
javac -classpath $samJarName":"$picardJarName *.java

echo "Generating Manifest file"
manifestFile=`pwd`"/BAMHeaderCreationManifest.txt"

echo -e "Class-Path: "$samJarName" "$picardJarName"\nMain-Class: Driver\n" > $manifestFile

echo "Building Jar file"

#cd ..
jar cvfm bam.header.creation.jar $manifestFile *.class ReferenceList.xml
#mv RegenerateBAMHeader.jar ./RegenerateBAMHeader
