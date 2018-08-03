#include "../../c/main/issm.h"
#include <cstddef>
#include <stdio.h>

//Android specific header includes: 
#include <jni.h>
#include <android/log.h>
#include <android/log.h>

//iOS specific header includes: 

namespace gov_nasa_jpl_issm
{
	/*Global variables{{{*/
	FemModel *fm;
	double* xyz; /*keep vertices information here*/
	/*}}}*/
	jint Initialize(JNIEnv *env, jclass clazz, jstring jsolution_type, jstring jabsfile, jstring jrelfile) /*{{{*/
	{
		
		/*arguments to constructor: */
		int argc; //arguments to constructor.
		char** argv = NULL;
		const char* issmname = "issm.exe";
		char *solution_type = NULL;
		char *absfile = NULL;
		char *relfile = NULL;
		ISSM_MPI_Comm    comm=1;

		/*log:*/
		__android_log_print(ANDROID_LOG_INFO, "Native","Initializing FemModel");

		/*retrieve from java machine: */
		solution_type = (char*)env->GetStringUTFChars(jsolution_type,0);
		absfile = (char*)env->GetStringUTFChars(jabsfile,0);
		relfile = (char*)env->GetStringUTFChars(jrelfile,0);

		/*creat arguments to call constructor for FemModel: */
		argc=4; 
		argv=(char**)malloc(argc*sizeof(char*));
		argv[0]=(char*)issmname;
		argv[1]=solution_type;
		argv[2]=absfile;
		argv[3]=relfile;
		
		/*call Model constructor passing in infile as File Descriptor parameter.*/
		fm  = new FemModel(argc,argv,comm);

		/*we'll need the toolkits activated right away to use matrices: */
		ToolkitsOptionsFromAnalysis(fm->parameters,NoneAnalysisEnum);

		/*release strings: */
		env->ReleaseStringUTFChars(jsolution_type, solution_type); //must realease the char*
		env->ReleaseStringUTFChars(jabsfile, absfile); //must realease the char*
		env->ReleaseStringUTFChars(jrelfile, relfile); //must realease the char*

		/*figure out size of solution: */
		__android_log_print(ANDROID_LOG_INFO, "Native","Number of elements");
		jint size = (jint) fm->elements->NumberOfElements();

		/*retrieve vertices x,y and z coordinates: */
		__android_log_print(ANDROID_LOG_INFO, "Native","Retrieving vertices");
		xyz=fm->vertices->ToXYZ();
		
		/*log: */
		__android_log_print(ANDROID_LOG_INFO, "Native","Done Initializing FemModel");

		return size;

	}
	/*}}}*/
	void Solve(JNIEnv *env, jclass clazz , jdouble alpha, jobject buf){ /*{{{*/

		int i,count;
		double x1,y1,z1,vel1;
		double x2,y2,z2,vel2;
		double x3,y3,z3,vel3;
		int    v1,v2,v3,eid;
		Patch* patch=NULL;
		
		/*log:*/
		__android_log_print(ANDROID_LOG_INFO, "Native","Solving ");

		/*retrieve buffer: */
		jdouble *dBuf = (jdouble *)env->GetDirectBufferAddress(buf);

		/*reset basal friction to what it was before: */
		__android_log_print(ANDROID_LOG_INFO, "Native","alpha %g ",alpha);
		
		__android_log_print(ANDROID_LOG_INFO, "Native","ok-1");

		InputDuplicatex(fm->elements,fm->nodes,fm->vertices,fm->loads,fm->materials,fm->parameters,AndroidFrictionCoefficientEnum,FrictionCoefficientEnum);
		__android_log_print(ANDROID_LOG_INFO, "Native","ok0");

		/*now scale friction by alpha: */
		InputScalex(fm->elements,fm->nodes,fm->vertices,fm->loads,fm->materials,fm->parameters,FrictionCoefficientEnum,alpha/100);
		__android_log_print(ANDROID_LOG_INFO, "Native","ok1");

		/*solve: */
		fm -> Solve();
		__android_log_print(ANDROID_LOG_INFO, "Native","ok2");

		/*retrieve results: */
		__android_log_print(ANDROID_LOG_INFO, "Native","Retrieving results ");
		//fm->elements->ProcessResultsUnits(); we are now in SI units
		patch=fm->elements->ResultsToPatch();

		/*sort out the velocities: */
		for(i=0;i<patch->numrows;i++){
			if ((patch->values[i*patch->numcols+0])==VelEnum){

				/*Each row of the Patch object is made of the following information: 
				  - the result enum_type, 
				  - the step and time, 
				  - the id of the element, 
				  - the interpolation type, 
				  - the vertices ids, 
				  - and the values at the nodes (could be different from the vertices)
				*/
				eid=(int)patch->values[i*patch->numcols+3]-1;
				v1=(int)patch->values[i*patch->numcols+5]; 
				x1=xyz[3*(v1-1)+0]; y1=xyz[3*(v1-1)+1]; z1=xyz[3*(v1-1)+2];
				
				v2=(int)patch->values[i*patch->numcols+6];
				x2=xyz[3*(v2-1)+0]; y2=xyz[3*(v2-1)+1]; z2=xyz[3*(v2-1)+2];
				
				v3=(int)patch->values[i*patch->numcols+7];
				x3=xyz[3*(v3-1)+0]; y3=xyz[3*(v3-1)+1]; z3=xyz[3*(v3-1)+2];

				vel1=patch->values[i*patch->numcols+8]; 
				vel2=patch->values[i*patch->numcols+9]; 
				vel3=patch->values[i*patch->numcols+10]; 

				/*plug into dBuf: */
				/*vertex 1: */
				dBuf[12*eid+0]=x1;
				dBuf[12*eid+1]=y1;
				dBuf[12*eid+2]=z1;

				/*vertex 2: */
				dBuf[12*eid+3]=x2;
				dBuf[12*eid+4]=y2;
				dBuf[12*eid+5]=z2;
			
				/*vertex 3: */
				dBuf[12*eid+6]=x3;
				dBuf[12*eid+7]=y3;
				dBuf[12*eid+8]=z3;
				
				/*values at 3 vertices: */
				dBuf[12*eid+9]=vel1;
				dBuf[12*eid+10]=vel2;
				dBuf[12*eid+11]=vel3;

			}
		}

		/*for(i=0;i<148;i++){
		__android_log_print(ANDROID_LOG_INFO, "Native","%g %g %g | %g %g %g | %g %g %g | %g %g %g\n",
				dBuf[12*i+0],dBuf[12*i+1],dBuf[12*i+2],
				dBuf[12*i+3],dBuf[12*i+4],dBuf[12*i+5],
				dBuf[12*i+6],dBuf[12*i+7],dBuf[12*i+8],
				dBuf[12*i+9],dBuf[12*i+10],dBuf[12*i+11]);
		}*/

		/*delete temporary data:*/
		delete patch;

	}/*}}}*/
	static JNINativeMethod method_table[] = /*{{{*/
	{
			{"createISSMModel"   ,"(Ljava/lang/String;Ljava/lang/String;Ljava/lang/String;)I"  , (void *) Initialize},
			{"solveISSMModel", "(DLjava/nio/DoubleBuffer;)V", (void *) Solve}
	};
	/*}}}*/
}

using namespace gov_nasa_jpl_issm;
extern "C" jint JNI_OnLoad(JavaVM* vm, void* reserved) /*{{{*/
{
    JNIEnv* env;
    if (vm->GetEnv(reinterpret_cast<void**>(&env), JNI_VERSION_1_6) != JNI_OK) {
        return -1;
    }
    else
    {
    	jclass clazz = env->FindClass("gov/nasa/jpl/issm/IssmJni");
    	if(clazz)
    	{
    		env->RegisterNatives(clazz, method_table, 3);
    		env->DeleteLocalRef(clazz);
    		return JNI_VERSION_1_6;
    	}
    	else return -1;
    }
}
/*}}}*/
