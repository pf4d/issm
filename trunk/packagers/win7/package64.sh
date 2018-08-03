#!/bin/bash

#get brand new project: 
rm -rf ISSM.aip  ISSM-SetupFiles ISSM.msi
cp ISSM.initial.aip ISSM.aip

#recover version: 
version=`svn info | grep Revision | awk '{printf("%s\n",$2);}'`

#get windows style path to files
export ISSM_DIR_WIN=`cygpath -m "$ISSM_DIR"`

echo "----------------------------------"
rm $ISSM_DIR/bin/*.m
find $ISSM_DIR/src/m -name '*.m' | xargs cp -t $ISSM_DIR/bin/
ls $ISSM_DIR/bin
echo "----------------------------------"

#build list of files to put into the installer: 
rm -rf ISSM.aic 
cat << END > ISSM.aic
;aic
SetVersion "$version"
SetPackageName "ISSM.msi"
END

ls $ISSM_DIR_WIN/scripts/*.bat startup.m | awk '{printf("AddFile APPDIR %s\n",$1);}' | sed 's/\//\\/g' >> ISSM.aic

cat << END >> ISSM.aic
AddFolder PersonalFolder $ISSM_DIR_WIN\test
AddFolder PersonalFolder $ISSM_DIR_WIN\examples
AddFolder APPDIR $ISSM_DIR_WIN\bin
AddFolder APPDIR $ISSM_DIR_WIN\lib
NewEnvironment -name ISSM_TESTS -value [test_Dir]
NewEnvironment -name ISSM_DIR -value [APPDIR]
NewEnvironment -name ISSM_DIR_WIN -value [APPDIR]
Save
Rebuild
END
#Not needed anymore? 
#DelEnvironment ISSM_TESTS
#DelEnvironment ISSM_DIR
#DelEnvironment ISSM_DIR_WIN

#Run installer: 
"C:/Program Files (x86)/Caphyon/Advanced Installer 10.8/bin/x86/AdvancedInstaller.com" /execute  ./ISSM.aip ./ISSM.aic

#Get rid of temporary files: 
cp ISSM-SetupFiles/ISSM.msi ./
rm -rf ISSM.aip ISSM-SetupFiles ISSM.aic

#To upload to website: 
#scp ISSM.msi websites:/home/larour/files/ISSM64.msi
