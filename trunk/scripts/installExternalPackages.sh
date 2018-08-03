#!/bin/bash

#Jenkins xml files for individual packages
EXTERNAL_TEST_FILE="$ISSM_DIR/nightlylog/results/external.xml"
mkdir -p $ISSM_DIR/nightlylog/results
echo "<testsuite tests=\"$NUMPACKAGES\">" > $EXTERNAL_TEST_FILE

source $ISSM_DIR/BuildConfig/externalpackages

#number of packages: 
NUMPACKAGES=$(($(echo $EXTERNALPACKAGES | wc -w )/2))
for ((i=1;i<=$NUMPACKAGES;i++))
do
	NUM1=$((2*$i-1))
	NUM2=$((2*$i))
	PACKAGENAME=$(echo $EXTERNALPACKAGES | cut -d " " -f $NUM1-$NUM1)
	PACKAGEINST=$(echo $EXTERNALPACKAGES | cut -d " " -f $NUM2-$NUM2)

	cd $ISSM_DIR/externalpackages/$PACKAGENAME

	#install if requested or if previous install has not been successful
	echo "======================================================";
	echo "       Installing $PACKAGENAME                        ";
	echo "======================================================";
	echo "<testcase classname=\"externalpackages\" name=\"$PACKAGENAME\">" >> $EXTERNAL_TEST_FILE
	./$PACKAGEINST |  tee compil.log
	if [ $? -ne 0 ]; then
		echo "======================================================";
		echo "    ERROR: installation of $PACKAGENAME failed        ";
		echo "======================================================";
		echo '<failure message="failure">' >> $EXTERNAL_TEST_FILE
		cat ./compil.log >> $EXTERNAL_TEST_FILE
		echo '</failure>' >> $EXTERNAL_TEST_FILE
		exit 1
	else
		echo '<success message="success">' >> $EXTERNAL_TEST_FILE
		cat ./compil.log >> $EXTERNAL_TEST_FILE
		echo '</success>' >> $EXTERNAL_TEST_FILE
		touch SUCCESS
	fi
	echo '</testcase>' >> $EXTERNAL_TEST_FILE
	source $ISSM_DIR/etc/environment.sh
done
echo '</testsuite>' >> $EXTERNAL_TEST_FILE
