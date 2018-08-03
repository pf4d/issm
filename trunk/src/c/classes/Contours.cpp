/*
 * \file Contours.cpp
 * \brief: Implementation of Contours class, derived from DataSet class.
 */

/*Headers: {{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./Contours.h"
#include "../shared/shared.h"
#include "./Contour.h"

using namespace std;
/*}}}*/

/*Object constructors and destructor*/
Contours::Contours(){/*{{{*/
	enum_type=ContoursEnum;
	return;
}
/*}}}*/
Contours::~Contours(){/*{{{*/
	return;
}
/*}}}*/

/*Numerics: */
int ExpWrite(Contours* contours,char* domainname){/*{{{*/

	/*I/O: */
	FILE* fid=NULL;
	Contour<double>* contour = NULL;

	/*open domain outline file for writing: */
	if((fid=fopen(domainname,"w"))==NULL) _error_("could not open domain file " << domainname); 

	for(int counter=0;counter<contours->Size();counter++){
		contour=(Contour<double>*)contours->GetObjectByOffset(counter);

		/*Write header: */
		fprintf(fid,"## Name:%s\n",domainname);
		fprintf(fid,"## Icon:0\n");
		fprintf(fid,"# Points Count	Value\n");
		fprintf(fid,"%u %s\n",contour->nods  ,"1.");
		fprintf(fid,"# X pos	Y pos\n");

		/*Write vertices: */
		for(int i=0;i<contour->nods;i++){
			fprintf(fid,"%lf\t%lf\n",contour->x[i],contour->y[i]);
		}

		/*Write blank line: */
		if(counter<contours->Size()-1) fprintf(fid,"\n");
	}

	/*close Exp file: */
	fclose(fid);

	return 1;
}/*}}}*/
