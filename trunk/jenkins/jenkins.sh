#!/bin/bash
#This bash script calls the nightlyrun.m matlab file to run our nightly test decks. 
#It then processes the results and sends an email to the Ice developpers.

echo "Cleaning up execution directory"
rm -rf $ISSM_DIR/execution/*
rm -rf $ISSM_DIR/nightlylog
mkdir  $ISSM_DIR/nightlylog

#Server URI
SERVER='https://ross.ics.uci.edu:8080'
#SERVER='http://ross.ics.uci.edu:8080'

#Get configuration
#Source config file{{{
if [ $# -ne 1 ];
then
	#no config file specified: exit
	echo "no config file specified. Exiting..." >&2 # Error message to stderr.
	exit 1
fi
if [ ! -f "$1" ]
then
	echo "File $1 not found!" >&2   # Error message to stderr.
	exit 1
fi 

#Defaults (to avoid -eq: unary operator expected)
EXAMPLES_TEST=0
MATLAB_TEST=0
PYTHON_TEST=0
JAVASCRIPT_TEST=0

#source configuration script
source $1;
#}}}
#Get Operating system (OS) name{{{
OS=$(uname -s)
if [[ $OS == CYGWIN* ]]; then 
	OS="win";
fi
#}}}

#Install ISSM
#Determinig Installation type depending on svn changes{{{
echo "======================================================";
echo "             Determining Installation type            "
echo "======================================================";
if [ -a $ISSM_DIR/svn_revision_old ]
then
	SVN_PREVIOUS=$(cat $ISSM_DIR/svn_revision_old)
	SVN_CURRENT=$SVN_REVISION_1
	echo "Previous revision number: $SVN_PREVIOUS"
	echo "Current  revision number: $SVN_CURRENT"

	#Get changes from jenkins itself (svn requires credentials)
	rm -rf changes
	wget $SERVER/job/$JOB_NAME/$BUILD_NUMBER/changes > /dev/null 2>&1

	#Process html page and get the list of files that has changed (tricky...)
	#cat changes | grep '="The file was modified"' | sed -e 's/.*<\/td><td><a>\(.*\)<\/a><\/td><td>.*/\1/' > $ISSM_DIR/TEMP
	#cat changes | grep 'document_edit' |sed -e 's/document_edit.png/document_edit.png\
		#/g' | sed -e 's/.*<\/a><\/td><td>\(.*\)<\/td><\/tr>.*/\1/' | grep -v 'document_edit.png' > $ISSM_DIR/TEMP
	cat changes  | tr " " "\n" | grep trunk |  sed -e 's/.*<a>\(.*\)<\/a>.*/\1/' > $ISSM_DIR/TEMP

	#Print list of files
	echo "   "
	echo "List of updated files"
	cat $ISSM_DIR/TEMP
	echo "   "

	#Do we need to reinstall externalpackages?
	echo "Determining installation type"
	if [ ! -z "$(cat $ISSM_DIR/TEMP | grep externalpackages)" ] ; then
		echo "  -- checking for changed externalpackages... yes";
		ISSM_EXTERNALPACKAGES="yes"
	else
		echo "  -- checking for changed externalpackages... no";
		ISSM_EXTERNALPACKAGES="no"
	fi

	#Do we need to reconfigure
	if [ ! -z "$(cat $ISSM_DIR/TEMP | grep -e "Makefile.am" -e "m4" )" ] ||
		[ ! -f "$ISSM_DIR/bin/issm.exe" ] ||
		[ "$ISSM_EXTERNALPACKAGES" == "yes" ] ;
	then
		echo "  -- checking for reconfiguration... yes";
		ISSM_RECONFIGURE="yes"
	else
		echo "  -- checking for reconfiguration... no";
		ISSM_RECONFIGURE="no"
	fi

	#Do we need to recompile
	if [ ! -z "$(cat $ISSM_DIR/TEMP | grep -e "\.cpp" -e "\.h" )" ] ||
		[ "$ISSM_RECONFIGURE" == "yes" ] ;
	then
		echo "  -- checking for recompilation... yes";
		ISSM_COMPILATION="yes"
	else
		echo "  -- checking for recompilation... no";
		ISSM_COMPILATION="no"
	fi
	#ISSM_COMPILATION="yes"

else
	echo "Previous revision not found, this must be a fresh install"
	echo "  -- checking for changed externalpackages... yes";
	echo "  -- checking for reconfiguration... yes";
	echo "  -- checking for recompilation... yes";
	ISSM_EXTERNALPACKAGES="yes"
	ISSM_RECONFIGURE="yes"
	ISSM_COMPILATION="yes"
fi
echo " "
echo "Recording current svn version: $SVN_REVISION_1"
echo $SVN_REVISION_1 > $ISSM_DIR/svn_revision_old
#}}}
#Install external packages    (ISSM_EXTERNALPACKAGES){{{

#number of packages: 
NUMPACKAGES=$(($(echo $EXTERNALPACKAGES | wc -w )/2))

#Jenkins xml files for individual packages
EXTERNAL_TEST_FILE="$ISSM_DIR/nightlylog/results/external.xml"
mkdir -p $ISSM_DIR/nightlylog/results
echo "<testsuite tests=\"$NUMPACKAGES\">" > $EXTERNAL_TEST_FILE

# Need a source here for when builds start midway through installation of externalpackages.
source $ISSM_DIR/etc/environment.sh

if [ "$OS" == "win" ]; then
	echo " == WINDOWS ENVIRONMENT DETECTED =="
	source $ISSM_DIR/externalpackages/windows/windows_environment.sh
fi

EXTERNALPACKAGES_FAILED=0;

for ((i=1;i<=$NUMPACKAGES;i++))
do
	NUM1=$((2*$i-1))
	NUM2=$((2*$i))
	PACKAGENAME=$(echo $EXTERNALPACKAGES | cut -d " " -f $NUM1-$NUM1)
	PACKAGEINST=$(echo $EXTERNALPACKAGES | cut -d " " -f $NUM2-$NUM2)

	cd $ISSM_DIR/externalpackages/$PACKAGENAME

	#install if requested or if previous install has not been successful
	if [ "$ISSM_EXTERNALPACKAGES" == "yes" ] || [ ! -d ./install -a ! -d ./install-javascript ] || ([ -d ./install ] && [ ! "$(ls -A ./install)" ]) || ([ -d ./install-javascript ] && [ ! "$(ls -A ./install-javascript)" ]); then
		echo "======================================================";
		echo "       Installing $PACKAGENAME                        ";
		echo "======================================================";

		./$PACKAGEINST &> compil.log
		if [ $? -ne 0 ] && [ "$PACKAGENAME" != "boost" ]; then
			echo "======================================================";
			echo "    ERROR: installation of $PACKAGENAME failed        ";
			echo "======================================================";
			echo "<testcase classname=\"externalpackages\" name=\"$PACKAGENAME\">" >> $EXTERNAL_TEST_FILE
			echo '<failure message="failure">External packages did not install right. Check it.' >> $EXTERNAL_TEST_FILE
			cat compil.log >> $EXTERNAL_TEST_FILE
			echo '</failure>' >> $EXTERNAL_TEST_FILE
			echo '</testcase>' >> $EXTERNAL_TEST_FILE
			EXTERNALPACKAGES_FAILED=1;
		else
			echo "<testcase classname=\"externalpackages\" name=\"$PACKAGENAME\"/>" >> $EXTERNAL_TEST_FILE
		fi
		source $ISSM_DIR/etc/environment.sh

		#If external package is rebuilt, we also need to recompile
		ISSM_RECONFIGURE="yes"
		ISSM_COMPILATION="yes"
	else
		echo "======================================================";
		echo "       Skipping $PACKAGENAME                          ";
		echo "======================================================";
		echo "<testcase classname=\"externalpackages\" name=\"$PACKAGENAME\"/>" >> $EXTERNAL_TEST_FILE
	fi
done
echo '</testsuite>' >> $EXTERNAL_TEST_FILE

if [ $EXTERNALPACKAGES_FAILED -eq 1 ]; then
	echo "===================================================================================================";
	echo "    ERROR: One or more of the externalpackages has failed. Skipping everything remaining steps.    ";
	echo "===================================================================================================";
	exit 1;
fi

# Source here to include any newly installed externalpackages on the path. 
source $ISSM_DIR/etc/environment.sh

if [ "$OS" == "win" ]; then
	echo " == WINDOWS ENVIRONMENT DETECTED =="
	source $ISSM_DIR/externalpackages/windows/windows_environment.sh
fi

#Set CXX/CC flags for JS runs after exnteralpackages to avoid conflicts during their compilation
CXX_PREVIOUS=$CXX
CC_PREVIOUS=$CC
if [ $JAVASCRIPT_TEST -eq 1 ]; then
	export CXX=em++
	export CC=emcc
	source $ISSM_DIR/externalpackages/emscripten/install/emsdk_env.sh
fi

#}}}
#ISSM compilation yes/no                (ISSM_COMPILATION) {{{
if [ "$ISSM_COMPILATION" == "yes" ]
then
	cd $ISSM_DIR
	if [ "$ISSM_RECONFIGURE" == "yes" ]
	then
		echo "======================================================";
		echo "             Cleaning up and reconfiguring            "
		echo "======================================================";
		make uninstall
		make distclean
		./scripts/automakererun.sh
		if [ $? -ne 0 ]; then 
			echo "autoreconf failed!"
			exit 1
		fi
		eval "./configure $ISSM_CONFIG"
		if [ $? -ne 0 ]; then 
			echo "ISSM configuration failed (see options below)"
			echo $ISSM_CONFIG
			echo "ISSM configuration failed!"
			exit 1
		fi
	fi

	#4: compile and install ISSM
	echo "======================================================";
	echo "                    Compiling ISSM                    "
	echo "======================================================";
	if [ $NUMCPUS_INSTALL -gt 1 ]; then
		echo "Making with " $NUMCPUS_INSTALL " cpus"
		make -j $NUMCPUS_INSTALL
	else
		make
	fi
	if [ $? -ne 0 ]; then 
		echo "ISSM_COMPILATION failed!"
		exit 1
	fi
	make install
elif [ "$ISSM_COMPILATION" == "no" ]
then
	echo "Skipping ISSM compilation"
else
	echo "ISSM_COMPILATION supported values are: yes and no. Exiting..." >&2 # Error message to stderr.
	exit 1
fi
#}}}

#Restore CXX/CC to their previous values 
export CXX=$CXX_PREVIOUS
export CC=$CC_PREVIOUS

#matlab tests
# {{{
if [ $MATLAB_TEST -eq 1 ]; then
#Launch all tests on different cpus
for (( i=1;i<=$NUMCPUS_RUN;i++ ))
do
	#Launch matlab and the nightly run script
	cat > $ISSM_DIR/nightlylog/matlab_run$i.m << EOF
	warning off %necessary to avoid a log of several Go for parallel runs
	try,
	$(if [ "$MATLAB_NROPTIONS" = ""  ]
	then
		echo "runme('output','nightly','rank',$i,'numprocs',$NUMCPUS_RUN);"
	else
		echo "runme($MATLAB_NROPTIONS,'output','nightly','rank',$i,'numprocs',$NUMCPUS_RUN);"
	fi
	)
	catch me,
		%An error occured, get report and exit
		message=getReport(me)
		directory=strsplit(pwd,'/');
		fid=fopen([issmdir '/nightlylog/matlaberror.log'], 'at');
		fprintf(fid,'\nMatlab error occured in: %s\n\n',directory{end});
		fprintf(fid,'%s',message);
		fclose(fid);
	end
	disp('MATLABEXITEDCORRECTLY');
	exit
EOF
	cd $ISSM_DIR/test/NightlyRun
	if [ "$OS" = "win" ]; then
		$MATLAB_PATH/bin/matlab -nojvm -nosplash -r "addpath $ISSM_DIR_WIN/src/m/dev; devpath; addpath $ISSM_DIR_WIN/nightlylog/; matlab_run$i" -logfile $ISSM_DIR_WIN/nightlylog/matlab_log$i.log &
	else
		$MATLAB_PATH/bin/matlab -nojvm -nosplash -r "addpath $ISSM_DIR/src/m/dev; devpath; addpath $ISSM_DIR/nightlylog/; matlab_run$i" -logfile $ISSM_DIR/nightlylog/matlab_log$i.log &
	fi
done

#wait until matlab closes
if [ "$OS" = "win" ]; then
	sleep 5;
	echo "Waiting for matlab on windows"
	pid=$(ps aux -W | grep MATLAB | awk '{printf("%s\n","MATLAB");}')
	echo '-----------------------------'
	echo "pid: $pid"
	echo '-----------------------------'
	while [ -n "$pid" ]
	do
		pid=$(ps aux -W | grep MATLAB | awk '{printf("%s\n","MATLAB");}')
		sleep 1;
	done
	echo "DONE!"
else
	wait
fi

#concatenate reports
cd $ISSM_DIR/nightlylog/
echo 'CHECKING NIGHTLYLOG DIRECTORY'
echo '-----------------------------'
ls -la
echo '-----------------------------'

rm matlab_log.log

for job in `jobs -p`
do
echo "Waiting on: $job"
    wait $job
done

for (( i=1;i<=$NUMCPUS_RUN;i++ ))
do
	cat matlab_log$i.log >> matlab_log.log
done

#filter out windows characters: 
cat matlab_log.log | tr -cd '\11\12\40-\176' > matlab_log.log2 && mv matlab_log.log2 matlab_log.log
fi
# }}}

#python tests
# {{{
if [ $PYTHON_TEST -eq 1 ]; then
#Launch all tests on different cpus
PYTHON_START_TIME=$(timer);
export PYTHONSTARTUP=$ISSM_DIR/src/m/dev/devpath.py
export PYTHONUNBUFFERED=1 #we don't want python to buffer otherwise issm.exe output is not captured
for (( i=1;i<=$NUMCPUS_RUN;i++ ))
do
	cd $ISSM_DIR/test/NightlyRun
	echo "--------------Running Python test for Rank $i---------------------"
	./runme.py --output=nightly --rank=$i --numprocs=$NUMCPUS_RUN $PYTHON_NROPTIONS &> $ISSM_DIR/nightlylog/python_log$i.log &
	echo "--------------Running Python test for Rank $i---------------------"
done

# concatenate reports
cd $ISSM_DIR/nightlylog/
rm python_log.log

for job in `jobs -p`
do
echo "Waiting on: $job"
    wait $job
done

for (( i=1;i<=$NUMCPUS_RUN;i++ ))
do
	echo "This is the concatenation phase for rank: python_log$i.log"
	cat python_log$i.log >> python_log.log
done
fi
# }}}

# Examples Test
# {{{
# This test will allow us to check on the status of the examples.
if [ $EXAMPLES_TEST -eq 1 ];
then
	FILE='runme.m'
	cd $ISSM_DIR/examples

	for dir in ./* ;
	do
		if [ -d "${dir}" ];
		then
		# Some of the examples are incomplete (on purpose). As such, we will have to populate the
		# missing steps in order to make sure that everything is working.
			echo "Testing directory example: $(basename $dir)"

			# Greenland is missing step 8
			if [ -z "$SED" ];
			then
				SED='sed'
			fi

			cd ${dir}

			if [ "${dir}" == "./Greenland" ];
			then
				STEP_EIGHT="\n	disp('   Step 8: Plotting exercise');\n	md = loadmodel('.\/Models\/Greenland.HistoricTransient');\n	figure\n	surfmb=[]; for i=2:201; surfmb=[surfmb ...\n		md.results.TransientSolution(i).SmbMassBalance]; end\n	subplot(3,1,1); plot([1:200],mean(surfmb));\n	title('Mean Surface mass balance');\n	vel=[]; for i=2:201; vel=[vel md.results.TransientSolution(i).Vel]; end\n	subplot(3,1,2); plot([1:200],mean(vel));\n	title('Mean Velocity');\n	volume=[]; for i=2:201; volume=[volume md.results.TransientSolution(i).IceVolume]; end\n	subplot(3,1,3); plot([1:200],volume);\n	title('Ice Volume'); xlabel('years');"

				$SED -i.bak 's/steps=\[1\];/steps=\[1:8\];\n\ntry\n/' $FILE
				$SED -i.bak "s/if any(steps==8)/&${STEP_EIGHT}/" $FILE
			elif [ "${dir}" == "./IceBridge" ];	
			then
				$SED -i.bak 's/steps=\[1\];/steps=\[1:5\];\n\ntry\n/' $FILE
			elif [ "${dir}" == "./IceflowModels" ];	
			then
				# Almost nothing to this example
				$SED -i.bak '1 s/^.*$/try\n\n&/' $FILE
			elif [ "${dir}" == "./ISMIP" ];	
			then
				# Eight steps... none of which are implmented in the script...
				$SED -i.bak '1 s/^.*$/try\n\n&/' $FILE
			elif [ "${dir}" == "./Inversion" ];	
			then
				$SED -i.bak 's/steps=\[1\];/steps=\[1:4\];\n\ntry\n/' $FILE
			elif [ "${dir}" == "./Jakobshavn" ];	
			then
				$SED -i.bak 's/steps=\[1\];/steps=\[1:4\];\n\ntry\n/' $FILE
			elif [ "${dir}" == "./Jakobshavn" ];	
			then
				$SED -i.bak 's/steps=\[1\];/steps=\[1:4\];\n\ntry\n/' $FILE
			elif [ "${dir}" == "./Pig" ];	
			then
				# Step 6 is needed
				STEP_SIX="\n disp('Needs work!')"
				$SED -i.bak 's/steps=\[1\];/steps=\[1:7\];\n!mv DomainOutline.bkp DomainOutline.exp;\n\ntry\n/' $FILE
				$SED -i.bak "s/if any(steps==6)/&${STEP_SIX}/" $FILE
			elif [ "${dir}" == "./PigSensitivity" ];	
			then
				# Step 4 is needed
				STEP_FOUR="\n disp('Needs work!')"
				$SED -i.bak 's/steps=\[1\];/steps=\[1:4\];\n\ntry\n/' $FILE
				$SED -i.bak "s/if any(steps==6)/&${STEP_FOUR}/" $FILE
			elif [ "${dir}" == "./SquareIceShelf" ];	
			then
				# Almost nothing to this example
				$SED -i.bak '1 s/^.*$/try\n\n&/' $FILE
			elif [ "${dir}" == "./UncertaintyQuantification" ];	
			then
				$SED -i.bak 's/steps=\[1\];/steps=\[1:7\];\n\ntry\n/' $FILE
			elif [ "${dir}" == "./Data" ];	
			then
				echo "Data directory is used by examples. No modifications required."
			else
				echo "Not implemented yet!"
				$SED -i.bak '1 s/^.*$/try\n\n&/' $FILE
			fi

			if [ "${dir}" == "./Data" ];
			then
				./Download.sh
			else
				LOG_FILE="matlab_log_$(basename $dir)_examples.log"
				echo "disp('SUCCESS');" >> $FILE
				echo 'catch me' >> $FILE
				echo 'message=getReport(me);' >> $FILE
				echo "fprintf('%s',message);" >> $FILE
				echo "disp('FAILURE');" >> $FILE
				echo 'end' >> $FILE

				$MATLAB_PATH/bin/matlab -nosplash -nodisplay -r "addpath $ISSM_DIR/src/m/dev; devpath; addpath $ISSM_DIR/nightlylog/; runme" -logfile $ISSM_DIR/nightlylog/$LOG_FILE
				echo "starting: $(basename $dir)" >> $ISSM_DIR/nightlylog/matlab_log_examples.log
				cat $ISSM_DIR/nightlylog/$LOG_FILE >> $ISSM_DIR/nightlylog/matlab_log_examples.log
				echo "finished: $(basename $dir)" >> $ISSM_DIR/nightlylog/matlab_log_examples.log
			fi
			cd ..
		fi
	done
fi
# }}}

#process logs to be junit compatible
#{{{
cd $ISSM_DIR/nightlylog/
source $ISSM_DIR/externalpackages/shell2junit/install/sh2ju.sh
juLogClean

if [ $MATLAB_TEST -eq 1 ]; then
	#number tests:
	numtests=`cat matlab_log.log  | grep "\-\-\-\-\-\-\-\-starting" | wc -l`
	testlist=`cat matlab_log.log  | grep "\-\-\-\-\-\-\-\-starting" | sed 's/----------------starting://g'  | sed 's/-//g'`

	#look through numtests:
	for i in `echo $testlist`
	do
		juLog  -test=MATLAB-$i -name=Error -error=ERROR awk "/starting:$i/{flag=1;next}/finished/{flag=0} flag{print}" matlab_log.log
		juLog  -test=MATLAB-$i -name=Failure -error=FAILURE awk "/starting:$i/{flag=1;next}/finished/{flag=0} flag{print}" matlab_log.log
	done
fi
if [ $PYTHON_TEST -eq 1 ]; then
	#number tests:
	numtests=`cat python_log.log  | grep "\-\-\-\-\-\-\-\-starting" | wc -l`
	testlist=`cat python_log.log  | grep "\-\-\-\-\-\-\-\-starting" | sed 's/----------------starting://g'  | sed 's/-//g'`

	#look through numtests:
	for i in `echo $testlist`
	do
		juLog  -test=PYTHON-$i -name=Error -error=ERROR awk "/starting:$i/{flag=1;next}/finished/{flag=0} flag{print}" python_log.log
		juLog  -test=PYTHON-$i -name=Failure -error=FAILURE awk "/starting:$i/{flag=1;next}/finished/{flag=0} flag{print}" python_log.log
	done
fi
if [ $EXAMPLES_TEST -eq 1 ];
then
	# Inexplicably, there are backspace chars in the error output that are causing issues
	$SED -i.bak 's///g' matlab_log_examples.log

	numtests=`cat matlab_log_examples.log  | grep "starting: " | wc -l`
	testlist=`cat matlab_log_examples.log   | grep "starting: " | sed 's/starting: //'`

	echo "Processing: $numtests"
	for i in `echo $testlist`
	do
		juLog  -test=Example-$i -name=Error -error=FAILURE awk "/starting: $i/{flag=1;next}/finished: $i/{flag=0} flag{print}" matlab_log_examples.log
	done
fi
#}}}
