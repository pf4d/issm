/*
 * FemModel.h: 
 */

#ifndef _FEMMODEL_H_
#define _FEMMODEL_H_

/*Headers:*/
/*{{{*/
#include "../toolkits/toolkits.h"
class DataSet;
class Parameters;
class NodeSets;
class Nodes;
class Vertices;
class Results;
class Constraints;
class Loads;
class Materials;
class Profiler;
class Elements;
#ifdef _HAVE_NEOPZ_
#include "./AdaptiveMeshRefinement.h"
#endif
/*}}}*/

class FemModel {

	/*no private members, as we need access to these datasets quite often!:*/

	public:

		int          analysis_counter;     //counter into analysis_type_list
		int         *analysis_type_list;   //list of analyses this femmodel is going to carry out
		int          nummodels;
		int          solution_type;

		Profiler*    profiler;             //keep time, cpu and mem statistics while we are running.

		Constraints *constraints;          //one set of constraints. each constraint knows which analysis_type it handles
		Elements    *elements;             //elements (one set for all analyses)
		Loads       *loads;                //one set of constraints. each constraint knows which analysis_type it handles
		Materials   *materials;            //one set of materials, for each element
		Nodes       *nodes;                //one set of nodes
		Parameters  *parameters;           //one set of parameters, independent of the analysis_type
		Results     *results;              //results that cannot be fit into the elements 
		Vertices    *vertices;             //one set of vertices

		#ifdef _HAVE_NEOPZ_
		AdaptiveMeshRefinement *amr;		  //adaptive mesh refinement object. It keeps coarse mesh and execute refinement process
		#endif

		/*constructors, destructors: */
		FemModel(int argc,char** argv,ISSM_MPI_Comm comm_init,bool trace=false);
		FemModel(char* rootpath, char* inputfilename, char* outputfilename, char* toolkitsfilename, char* lockfilename, char* restartfilename, ISSM_MPI_Comm incomm, int solution_type,IssmPDouble* X);
		~FemModel();

		/*Methods:*/
		void CheckPoint(void);
		void CleanUp(void);
		FemModel* copy();
		void Echo();
		void InitFromFiles(char* rootpath, char* inputfilename, char* outputfilename, char* petscfilename, char* lockfilename, char* restartfilename, const int solution_type,bool trace,IssmPDouble* X=NULL);
		void InitFromFids(char* rootpath, FILE* IOMODEL, FILE* toolkitsoptionsfid, int in_solution_type, bool trace, IssmPDouble* X=NULL);
		void Marshall(char** pmarshalled_data, int* pmarshalled_data_size, int marshall_direction);
		void Restart(void);
		void SetCurrentConfiguration(int configuration_type);
		void SetCurrentConfiguration(int configuration_type,int analysis_type);
		int  Size(void);
		void SolutionAnalysesList(int** panalyses,int* pnumanalyses,IoModel* iomodel,int solutiontype);
		void Solve(void);

		/*Modules*/ 
		void BalancethicknessMisfitx(IssmDouble* pV);
		void CalvingRateDevx();
		void CalvingRateLevermannx();
		void DeviatoricStressx();
		void Divergencex(IssmDouble* pdiv);
		void ElementOperationx(void (Element::*function)(void));
		void ElementResponsex(IssmDouble* presponse,int response_enum);
		void FloatingAreax(IssmDouble* pV);
		void GetInputLocalMinMaxOnNodesx(IssmDouble** pmin,IssmDouble** pmax,IssmDouble* ug);
		void GroundedAreax(IssmDouble* pV);
		void IceMassx(IssmDouble* pV);
		void IceVolumex(IssmDouble* pV);
		void IceVolumeAboveFloatationx(IssmDouble* pV);
		void MassFluxx(IssmDouble* presponse);
		void MaxAbsVxx(IssmDouble* presponse);
		void MaxAbsVyx(IssmDouble* presponse);
		void MaxAbsVzx(IssmDouble* presponse);
		void MaxDivergencex(IssmDouble* pdiv);
		void MaxVelx(IssmDouble* presponse);
		void MaxVxx(IssmDouble* presponse);
		void MaxVyx(IssmDouble* presponse);
		void MaxVzx(IssmDouble* presponse);
		void MinVelx(IssmDouble* presponse);
		void MinVxx(IssmDouble* presponse);
		void MinVyx(IssmDouble* presponse);
		void MinVzx(IssmDouble* presponse);
		void ResetLevelset();
		void StrainRateparallelx();
		void StrainRateperpendicularx();
		void StressIntensityFactorx();
		void TotalFloatingBmbx(IssmDouble* pFbmb);
		void TotalGroundedBmbx(IssmDouble* pGbmb);
		void TotalSmbx(IssmDouble* pSmb);
		#ifdef  _HAVE_DAKOTA_
		void DakotaResponsesx(double* d_responses,char** responses_descriptors,int numresponsedescriptors,int d_numresponses);
		#endif
		void CostFunctionx(IssmDouble* pJ,IssmDouble** pJlist,int* pn);
		void OutputControlsx(Results **presults);
		void RequestedDependentsx(void);
		void RequestedOutputsx(Results **presults,char** requested_outputs, int numoutputs,bool save_results=true);
		void RequestedOutputsx(Results **presults,int* requested_outputs, int numoutputs,bool save_results=true);
		void Responsex(IssmDouble* presponse,int response_descriptor_enum);
		void Responsex(IssmDouble* presponse,const char* response_descriptor);
		void SurfaceAbsMisfitx( IssmDouble* pJ);
		void ThicknessAbsGradientx( IssmDouble* pJ);
		void ThicknessPositivex(IssmDouble* pJ);
		#ifdef _HAVE_GIAIVINS_
		void Deflection(Vector<IssmDouble>* wg,Vector<IssmDouble>* dwgdt, IssmDouble* x, IssmDouble* y);
		#endif
		#ifdef _HAVE_ESA_
		void EsaGeodetic2D(Vector<IssmDouble>* pUp, Vector<IssmDouble>* pNorth, Vector<IssmDouble>* pEast, IssmDouble* xx, IssmDouble* yy); 
		void EsaGeodetic3D(Vector<IssmDouble>* pUp, Vector<IssmDouble>* pNorth, Vector<IssmDouble>* pEast, IssmDouble* latitude, IssmDouble* longitude, IssmDouble* radius, IssmDouble* xx, IssmDouble* yy, IssmDouble* zz); 
		#endif
		#ifdef _HAVE_SEALEVELRISE_
		void SealevelriseEustatic(Vector<IssmDouble>* pSgi, IssmDouble* peustatic, IssmDouble* latitude, IssmDouble* longitude, IssmDouble* radius);
		void SealevelriseNonEustatic(Vector<IssmDouble>* pSgo, Vector<IssmDouble>* pSg_old, IssmDouble* latitude, IssmDouble* longitude, IssmDouble* radius,bool verboseconvolution);
		void SealevelriseRotationalFeedback(Vector<IssmDouble>* pSgo_rot, Vector<IssmDouble>* pSg_old, IssmDouble* latitude, IssmDouble* longitude, IssmDouble* radius);
		void SealevelriseGeodetic(Vector<IssmDouble>* pUp, Vector<IssmDouble>* pNorth, Vector<IssmDouble>* pEast, Vector<IssmDouble>* pSg_old, IssmDouble* latitude, IssmDouble* longitude, IssmDouble* radius, IssmDouble* xx, IssmDouble* yy, IssmDouble* zz); 
		IssmDouble SealevelriseOceanAverage(Vector<IssmDouble>* Sg);
		#endif
		void HydrologyEPLupdateDomainx(IssmDouble* pEplcount);
		void TimeAdaptx(IssmDouble* pdt);
		void UpdateConstraintsExtrudeFromBasex();
		void UpdateConstraintsExtrudeFromTopx();
		void UpdateConstraintsL2ProjectionEPLx(IssmDouble* pL2count);
		void UpdateConstraintsx(void);
		int  UpdateVertexPositionsx(void);

		#ifdef _HAVE_JAVASCRIPT_
		FemModel(IssmDouble* buffer, int buffersize, char* toolkits, char* solution, char* modelname,ISSM_MPI_Comm incomm, bool trace=false);
		void CleanUpJs(char** poutput, size_t* psize);
		void InitFromBuffers(char* buffer, int buffersize, char* toolkits, int solution_type,bool trace,IssmPDouble* X=NULL);
		#endif

		#ifdef _HAVE_NEOPZ_
		/*Adaptive mesh refinement methods*/
		void InitializeAdaptiveRefinement(void);
		void ReMesh(void);
		void BedrockFromMismipPlus(void);
		void AdjustBaseThicknessAndMask(void);
		void GetMesh(Vertices* femmodel_vertices,Elements* femmodel_elements,IssmDouble** px, IssmDouble** py, IssmDouble** pz, int** pelementslist);
		int GetElementsWidth(){return 3;};//just tria elements in this first version
		void ExecuteRefinement(int &numberofvertices,int &numberofelements,IssmDouble** px,IssmDouble** py,IssmDouble** pz,int** pelementslist);
		void GetGroundediceLevelSet(IssmDouble** pmasklevelset);
		void CreateVertices(int newnumberofvertices,int newnumberofelements,int elementswidth,int* newelementslist,int* my_vertices,IssmDouble* newx,IssmDouble* newy,IssmDouble* newz,Vertices* vertices);
		void CreateElements(int newnumberofelements,int elementswidth,int* newelementslist,bool* my_elements,Elements* elements);
		void CreateMaterials(int newnumberofelements,bool* my_elements,Materials* materials);
		void CreateNodes(int newnumberofvertices,int* my_vertices,int nodecounter,int analysis_enum,Nodes* nodes);
		void CreateConstraints(int newnumberofvertices,int newnumberofelements,int nodecounter,int constraintcounter,IssmDouble* newx,IssmDouble* newy,int* my_vertices,Constraints* constraints);
		void InterpolateInputs(Vertices* newfemmodel_vertices,Elements* newfemmodel_elements);
		void UpdateElements(int newnumberofelements,int* newelementslist,bool* my_elements,int nodecounter,int analysis_counter,Elements* newelements);
		void ElementsAndVerticesPartitioning(int& newnumberofvertices,int& newnumberofelements,int& elementswidth,int* newelementslist,bool** pmy_elements,int** pmy_vertices);
		void WriteMeshInResults(void);
		void SetRefPatterns(void);
		#endif
};
		

#endif
