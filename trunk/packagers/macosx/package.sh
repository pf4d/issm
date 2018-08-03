#!/bin/bash

echo "modify generic" 
cd $ISSM_DIR/bin
cat generic_static.m | sed -e "s/generic_static/generic/g" > generic.m
echo "move mpiexec to bin" 
cp ../externalpackages/mpich/install/bin/mpiexec .
cp ../externalpackages/mpich/install/bin/hydra_pmi_proxy .

#Check that test101 runs
cd $ISSM_DIR/test/NightlyRun
rm matlab.log
/Applications/MATLAB_R2015b.app/bin/matlab -nojvm -nosplash -r "try, addpath $ISSM_DIR/bin $ISSM_DIR/lib; runme('id',101);exit; catch me,fprintf('%s',getReport(me)); exit; end" -logfile matlab.log

if [[ $(cat matlab.log | grep -c SUCCESS) -lt 10 ]]; then
	echo "test101 FAILED"
	exit 1;
else
	echo "test101 passed"
fi

#Package using the Package Maker from OSX, driven by command line.
tarball_name='issm-mac-static_build.tar.gz'

echo "Cleanup first" 
cd $ISSM_DIR
rm $tarball_name

echo "Creating tarball: ${tarball_name}"
cd $ISSM_DIR
rm -rf trunk
mkdir trunk
#Need script to download data
cp -rf bin lib test examples scripts trunk/
tar -czf $tarball_name trunk
ls -lah $tarball_name

echo "Shipping binaries to website"

# We're using public key authentication method to upload the tarball The
# following lines check to see if the SSH Agent is running. If not, then it is
# started and relevant information is forwarded to a script.
pgrep "ssh-agent" > /dev/null
if [ $? -ne 0 ]; then
	echo "SSH Agent is not running. Starting it..."
	ssh-agent > ~/.ssh/agent.sh
else
	echo "SSH Agent is running..."
fi

source ~/.ssh/agent.sh
ssh-add ~/.ssh/macosx-bins_richese-to-ross

scp $tarball_name ross.ics.uci.edu:/var/www/html/$tarball_name

if [ $? -ne 0 ]; then
	echo "The upload failed."
	echo "Perhaps the SSH Agent was started by some other means."
	echo "Try killing the agent and running again."
fi
