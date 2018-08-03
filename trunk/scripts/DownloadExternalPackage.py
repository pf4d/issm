#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-
#
# TODO: Improve upon extension recognition by checking for mismatches in found targets
# and specified local file.
#

# imports {{{
import os,sys,re
import urllib
from HTMLParser import HTMLParser
from urllib import FancyURLopener
# }}}
class MyHTMLParser(HTMLParser): #{{{

    def __init__(self, pattern):
        HTMLParser.__init__(self)
        self.matcher = re.compile(pattern) 
        self.targets = []

    def handle_starttag(self, tag, attrs):
        for i in attrs:
            if "href" == i[0] and str(self.matcher.match(i[1])) != "None":
                self.targets.append(i[1])
#}}}
def main(argv=None): # {{{
    # Separates the URL into a directory and the file or pattern based on the
    # last appearance of '/'.
	if len(sys.argv) > 1:
	    pivot = sys.argv[1].rfind("/")
	    url = (sys.argv[1])[:pivot]
	    pivot += 1
	    find = (sys.argv[1])[pivot:]
	else:
	    print "******************************************************************************************************************************"
	    print "* Invalid input!                                                                                                             *"
	    print "*                                                                                                                            *"
	    print "* Try: 'DownloadExternalPackage.py url [localFile]'                                                                          *"
	    print "*                                                                                                                            *"
	    print "* Where 'URL' is the URL with an explicit package name or the URL followed by the truncated package name. And 'localFile' is *"
	    print "* the file name (including extension) that you would like to save as.                                                        *"
	    print "*                                                                                                                            *"
	    print "* Examples:                                                                                                                  *" 
	    print "*                                                                                                                            *"
	    print "* DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/petsc-2.3.2-p3.tar.gz' 'petsc-2.3.2-p3.tar.gz' *"
	    print "*                                                                                                                            *"
	    print "*     This is the old style and the safest way to download a package.                                                        *"
	    print "*                                                                                                                            *"
	    print "* DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/libtool' 'libtool.tar.gz'                      *"
	    print "*                                                                                                                            *"
	    print "*     This is the new style. For packages like 'Libtool', which we never expect to be using multiple versions, this will     *"
	    print "*     download the most recent version and save it as the generic 'libtool.tar.gz'.                                          *"
	    print "*                                                                                                                            *"
	    print "* DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/gsl-1.' 'gsl-1.15.tar.gz'                      *"
	    print "*                                                                                                                            *"
	    print "*     This is the new style. This is a demonstration of how this script can be used to disambiguate a package name if there  *"
	    print "*     are more than once package matching 'gsl-'.                                                                            *"
	    print "*                                                                                                                            *"
	    print "* DownloadExternalPackage.py 'http://issm.jpl.nasa.gov/files/externalpackages/libtool'                                       *"
	    print "*                                                                                                                            *"
	    print "*     This is the new style. This will download a package with 'libtool' as a prefix and save it as its canonical name.      *"
	    print "*                                                                                                                            *"
	    print "*                                                                                                                            *"
	    print "******************************************************************************************************************************"
	
	if len(sys.argv) > 2:
	    localFile=sys.argv[2]
	    print "Downloaded file will be saved as: " + localFile
	else:
	    localFile = None
	    print "Downloaded file will saved with the same file name."
	
	
	print "Looking for: " + find
	
	# As an extra precaution, if no extension is given for a particular package
	# such as '.../libtool', then ensure that files found are of appropriate
	# file extensions. 
	#
	# WARNING: The external packages directory includes executable binaries with
	# '.exe' extensions. As such, '.exe' is an acceptable suffix, but this is 
	# inherently dangerous since this script can be used to download from any
	# valid website. Furthermore, if an individual attempts a "man-in-the-middle"  
	# attack, then the user would be capable of downloading executables from 
	# an untrusted source.
	pattern = find + "[\w.-]*(\.tar\.gz|tar\.gz2|tgz|zip|exe)?"
	parser = MyHTMLParser(pattern)
	
	# Creates a 'FancyURL' which allows the script to fail gracefully by catching
	# HTTP error codes 30X and several 40X(where 'X' is a natural number).
	urlObject = FancyURLopener()
	obj = urlObject.open(url)
	parser.feed(obj.read())
	
	# If a file pattern was used to describe the file that should be downloaded,
	# then there is the potential for multiple file matches. Currently, the script
	# will detect this ambiguity and print out all the matches, while informing 
	# the user that he must refine his search.
	#
	# TODO: Prompt the user to select from a list his/her preferred target.
	if len(parser.targets) > 1:
	    print "Could not resolve your download due to the number of hits."
	    print "Refine your search."
	    for i in parser.targets:
	        print i
	
	elif len(parser.targets) == 1:
	    print "Found: " + parser.targets[0]
	    url += "/" + parser.targets[0]
	
	    if localFile is None:
	        if os.path.exists(parser.targets[0]): 
	            print "File " + parser.targets[0] + " already exists and will not be downloaded..."
	        else:
	            urllib.urlretrieve(url, parser.targets[0])
	            print "File saved as: " + parser.targets[0]
	    else:
	        if os.path.exists(localFile): 
	            print "File "+ localFile +" already exists and will not be downloaded..."
	        else:
	            if parser.targets[0] == localFile:
	                print "File found and destination match."
	            elif parser.matcher.match(localFile) != "None":
	                print "File found matches destination pattern."
	            else:
	                print "WARNING: the file found \'" + parser.targets[0] + "\' does not match \'" + localFile + "\'"
	                print "Ensure the downloaded version is suitable."
	
	            urllib.urlretrieve(url, localFile)
	            print "File saved as: " + localFile
	
	else:
	    print "No matches found!"
	
	obj.close()
# End 'main' function. }}}

if __name__ == "__main__":
    main()
