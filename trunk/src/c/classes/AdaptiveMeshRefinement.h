#ifndef ADAPTIVEMESHREFINEMENT
#define ADAPTIVEMESHREFINEMENT

/*Includes*/
/*{{{*/
/*Common includes*/
#include <iostream>
#include <fstream>
#include <string>

/*NeoPZ includes*/
/*REAL and STATE definitions, NeoPZ variables itapopo should be read by NeoPZ's config.h*/
#ifndef REFPATTERNDIR
	# define REFPATTERNDIR "/home/santos/trunk-jpl/externalpackages/neopz/install/include/refpatterns"
#endif

#ifndef REALdouble
	#define REALdouble
#endif

#ifndef STATEdouble
	#define STATEdouble
#endif

#include <pzreal.h>
#include <pzsave.h>
#include <pzgmesh.h>
#include <pzvec.h>
#include <pzeltype.h>

#include <TPZRefPatternTools.h>
#include <TPZRefPatternDataBase.h>
#include <TPZRefPattern.h>

#include <tpzchangeel.h>
#include <TPZGeoElement.h>
#include <pzreftriangle.h>
#include <tpzgeoelrefpattern.h>
/*}}}*/

class AdaptiveMeshRefinement : public TPZSaveable {

public:

	/*Public methods*/
	/* Constructor, destructor etc*/
	AdaptiveMeshRefinement();																// Default constructor
	AdaptiveMeshRefinement(const AdaptiveMeshRefinement &cp); 					// Copy constructor
	AdaptiveMeshRefinement & operator= (const AdaptiveMeshRefinement &cp);	// Operator of copy
	virtual ~AdaptiveMeshRefinement();													// Destructor

    /*Savable methods*/
	virtual int ClassId() const;                                            // ClassId to save the class
   virtual void Read(TPZStream &buf, void *context);								// Read this class
   virtual void Write(TPZStream &buf, int withclassid);                    // Write this class, using ClassId to identify
    
	/*General methods*/
	void CleanUp();																			// Clean all attributes
	void Initialize();																		// Initialize the attributes with NULL and values out of usually range
	void SetLevelMax(int &h);                                               // Define the max level of refinement
   void SetRegions(double &D1,double Dhmax);										// Define the regions which will be refined
	void SetElementWidth(int &width);                                       // Define elements width
	void ExecuteRefinement(int &type_process,double *vx,double *vy,double *masklevelset,int &nvertices,int &nelements,int &nsegments,double** px,double** py,double** pz,int** pelements,int** psegments=NULL);					// A new mesh will be created and refined. This returns the new mesh
	void CreateInitialMesh(int &nvertices,int &nelements,int &nsegments,int &width,double* x,double* y,double* z,int* elements,int* segments=NULL); // Create a NeoPZ geometric mesh by coords and elements
	TPZGeoMesh* CreateRefPatternMesh(TPZGeoMesh* gmesh);
	void CheckMesh(int &nvertices,int &nelements,int &nsegments,int &width,double** px,double** py,double** pz,int** pelements,int** psegments=NULL); // Check the consistency of the mesh

private:

	/*Private attributes*/
   int elementswidth;                                                      // Geometric nodes for element: 3 == Tria, 4 == Tetra, 6 == Penta
   int levelmax;                                                           // Max level of refinement
	double regionlevel1;																		// Region which will be refined with level 1
	double regionlevelmax;																	// Region which will be refined with level max
	TPZGeoMesh *fathermesh;																	// Father Mesh is the entire mesh without refinement
	TPZGeoMesh *previousmesh;																// Previous mesh is a refined mesh of last step

	/*Private methods*/
   void RefinementProcess(TPZGeoMesh *gmesh,std::vector<TPZVec<REAL> > &GLvec);  // Start the refinement process
	void RefineMesh(TPZGeoMesh *gmesh, std::vector<int> &ElemVec); 					// Refine the elements in ElemVec
   void RefineMeshToAvoidHangingNodes(TPZGeoMesh *gmesh);                        // Refine the elements to avoid hanging nodes
	void SetElementsToRefine(TPZGeoMesh *gmesh,std::vector<TPZVec<REAL> > &GLvec,int &hlevel, std::vector<int> &ElemVec); 	//Define wich elements will be refined
   void TagAllElements(TPZGeoMesh *gmesh,std::vector<int> &ElemVec);				 // This tag all elements to be refined, that is, refine all elements
   void TagElementsNearGroundingLine(TPZGeoMesh *gmesh,std::vector<TPZVec<REAL> > &GLvec,int &hlevel,std::vector<int> &ElemVec);    // This tag elements near the grounding line
   void CalcGroundingLinePosition(double *masklevelset,std::vector<TPZVec<REAL> > &GLvec);	// Calculate the grounding line position using previous mesh
	void GetMesh(TPZGeoMesh *gmesh,int &nvertices,int &nelements,int &nsegments,double** px,double** py,double** pz,int** pelements,int** psegments=NULL); // Return coords and elements in ISSM data structure
   inline int GetElemMaterialID(){return 1;}                               // Return element material ID
   inline int GetBoundaryMaterialID(){return 2;}                           // Return segment (2D boundary) material ID
};

#endif
