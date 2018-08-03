#!/bin/bash

#Some cleanup
\rm -rf install

#Add cvs repository
export CVSROOT=':pserver:cvsanon@mitgcm.org:/u/gcmpack'

#Download code from server
echo loging into MITgcm CVS server
echo enter MITgcm CVS password: cvsanon
#Login only the first time
#cvs login
cvs co -P MITgcm_code

#Move
mv MITgcm install
