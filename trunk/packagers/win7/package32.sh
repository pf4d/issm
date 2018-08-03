#!/bin/bash

#get brand new project: 
rm -rf ISSM.aip  ISSM-SetupFiles ISSM.msi
cp ISSM.initial.aip ISSM.aip

#get windows style path to files
export ISSM_DIR_WIN=`cygpath -m "$ISSM_DIR"`

#build list of files to put into the installer: 
rm -rf ISSM.aic 
cat << END > ISSM.aic
;aic
SetVersion "1.0"
SetPackageName "ISSM.msi"
END

ls $ISSM_DIR_WIN/scripts/*.bat startup.m  | awk '{printf("AddFile APPDIR %s\n",$1);}' | sed 's/\//\\/g' >> ISSM.aic

cat << END >> ISSM.aic
AddFolder PersonalFolder $ISSM_DIR_WIN\test
AddFolder APPDIR $ISSM_DIR_WIN\bin
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
