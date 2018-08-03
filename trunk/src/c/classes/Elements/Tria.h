/*! \file Tria.h 
 *  \brief: header file for tria object
 */

#ifndef _TRIA_H_
#define _TRIA_H_

/*Headers:*/
/*{{{*/
#include "./Element.h"
#include "./ElementHook.h"
#include "./TriaRef.h"
class Parameters;
class Inputs;
class IoModel;
class Results;
class Node;
class Material;
class Matpar;
class Seg;
class ElementMatrix;
class ElementVector;
class Vertex;
class GaussTria;

#include "../../shared/Exceptions/exceptions.h"
#include "../../shared/Enum/Enum.h"
/*}}}*/

class Tria: public Element,public ElementHook,public TriaRef{

	public:

		/*Tria constructors, destructors {{{*/
		Tria(){};
		Tria(int tria_id,int tria_sid,int i, IoModel* iomodel,int nummodels);
		~Tria();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Object *copy();
		void    Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction);
		int     ObjectEnum();
		/*}}}*/
		/*Update virtual functions resolution: {{{*/
		#ifdef _HAVE_DAKOTA_
		void  InputUpdateFromMatrixDakota(IssmDouble* matrix, int nows, int ncols, int name, int type);
		void  InputUpdateFromVectorDakota(IssmDouble* vector, int name, int type);
		#endif
		void  InputUpdateFromIoModel(int index, IoModel* iomodel);
		void  InputUpdateFromVector(IssmDouble* vector, int name, int type);
		/*}}}*/
		/*Element virtual functions definitions: {{{*/
		void        AverageOntoPartition(Vector<IssmDouble>* partition_contributions,Vector<IssmDouble>* partition_areas,IssmDouble* vertex_response,IssmDouble* qmu_part);
		void			CalvingRateDev();
		void			CalvingRateLevermann();
		IssmDouble  CharacteristicLength(void);
		void        ComputeBasalStress(Vector<IssmDouble>* sigma_b);
		void        ComputeDeviatoricStressTensor();
		void        ComputeEsaStrainAndVorticity();
		void        ComputeSigmaNN();
		void        ComputeStressTensor();
		void        ComputeSurfaceNormalVelocity();
		void        Configure(Elements* elements,Loads* loads,Nodes* nodesin,Vertices* verticesin,Materials* materials,Parameters* parameters);
		void        ControlInputSetGradient(IssmDouble* gradient,int enum_type,int control_index);
		void        ControlToVectors(Vector<IssmPDouble>* vector_control, Vector<IssmPDouble>* vector_gradient,int control_enum);
		int         EdgeOnBaseIndex();
		void        EdgeOnBaseIndices(int* pindex1,int* pindex);
		int         EdgeOnSurfaceIndex();
		void        EdgeOnSurfaceIndices(int* pindex1,int* pindex);
		void        ElementResponse(IssmDouble* presponse,int response_enum);
		void        ElementSizes(IssmDouble* hx,IssmDouble* hy,IssmDouble* hz);
		int         FiniteElement(void);
		IssmDouble  FloatingArea(void);
		void        FSContactMigration(Vector<IssmDouble>* vertexgrounded,Vector<IssmDouble>* vertexfloating);
		Element*    GetBasalElement(void){_error_("not implemented yet");};
		void        GetLevelsetPositivePart(int* point1,IssmDouble* fraction1,IssmDouble* fraction2, bool* mainlynegative,IssmDouble* levelsetvalues);
		void        GetGroundedPart(int* point1,IssmDouble* fraction1, IssmDouble* fraction2,bool* mainlyfloating);
		IssmDouble  GetGroundedPortion(IssmDouble* xyz_list);
		void	      GetIcefrontCoordinates(IssmDouble** pxyz_front,IssmDouble* xyz_list,int levelsetenum);
		void	      GetLevelCoordinates(IssmDouble** pxyz_front,IssmDouble* xyz_list,int levelsetenum,IssmDouble level);
		int         GetNodeIndex(Node* node);
		int         GetNumberOfNodes(void);
		int         GetNumberOfNodes(int enum_type);
		int         GetNumberOfVertices(void);
		void        GetSolutionFromInputsOneDof(Vector<IssmDouble>* solution,int enum_type);
		Element*    GetUpperElement(void){_error_("not implemented yet");};
		void        GetVectorFromControlInputs(Vector<IssmDouble>* gradient,int control_enum,int control_index,const char* data,bool onsid);
		void        GetVerticesCoordinatesBase(IssmDouble** pxyz_list);
		void        GetVerticesCoordinatesTop(IssmDouble** pxyz_list);
		IssmDouble  GroundedArea(void);
		bool        HasEdgeOnBase();
		bool        HasEdgeOnSurface();
		IssmDouble  IceMass(void);
		IssmDouble  IceVolume(void);
		IssmDouble  IceVolumeAboveFloatation(void);
		void        InputControlUpdate(IssmDouble scalar,bool save_parameter);
		void        InputDepthAverageAtBase(int enum_type,int average_enum_type);
		void        InputExtrude(int enum_type,int start){_error_("not implemented"); /*For penta only*/};
		void        InputScale(int enum_type,IssmDouble scale_factor);
		bool	   	IsFaceOnBoundary(void);
		bool	   	IsIcefront(void);
		bool        IsNodeOnShelfFromFlags(IssmDouble* flags);
		bool        IsOnBase();
		bool        IsOnSurface();
		bool        IsZeroLevelset(int levelset_enum);
		IssmDouble  Masscon(IssmDouble* levelset);
		IssmDouble  MassFlux(IssmDouble* segment);
		IssmDouble  MassFlux(IssmDouble x1,IssmDouble y1, IssmDouble x2, IssmDouble y2,int segment_id);
		void        MaterialUpdateFromTemperature(void){_error_("not implemented yet");};
		IssmDouble  Misfit(int modelenum,int observationenum,int weightsenum);
		IssmDouble  MisfitArea(int weightsenum);
		int         NodalValue(IssmDouble* pvalue, int index, int natureofdataenum);
		int         NumberofNodesPressure(void);
		int         NumberofNodesVelocity(void);
		void        PotentialUngrounding(Vector<IssmDouble>* potential_sheet_ungrounding);
		int         PressureInterpolation();
		void        ReduceMatrices(ElementMatrix* Ke,ElementVector* pe);
		void        ResetFSBasalBoundaryCondition(void);
		void        ResetHooks();
		void        ResetLevelsetFromSegmentlist(IssmDouble* segments,int numsegments);
		void        SetControlInputsFromVector(IssmDouble* vector,int control_enum,int control_index);
		void        SetCurrentConfiguration(Elements* elements,Loads* loads,Nodes* nodes,Materials* materials,Parameters* parameters);
	   Element*    SpawnBasalElement(void);
		Element*    SpawnTopElement(void);
		void			StrainRateparallel();
		void			StrainRateperpendicular();
		void        StressIntensityFactor(void){_error_("not implemented yet");};
		IssmDouble  SurfaceArea(void);
		int         TensorInterpolation();
		IssmDouble  TimeAdapt();
		IssmDouble  TotalFloatingBmb(void);
		IssmDouble  TotalGroundedBmb(void);
		IssmDouble  TotalSmb(void);
		void        Update(int index, IoModel* iomodel,int analysis_counter,int analysis_type,int finitelement);
		int         UpdatePotentialUngrounding(IssmDouble* vertices_potentially_ungrounding,Vector<IssmDouble>* vec_nodes_on_iceshelf,IssmDouble* nodes_on_iceshelf);
		void        ValueP1DerivativesOnGauss(IssmDouble* dvalue,IssmDouble* values,IssmDouble* xyz_list,Gauss* gauss);
		void        ValueP1OnGauss(IssmDouble* pvalue,IssmDouble* values,Gauss* gauss);
		int         VelocityInterpolation();
		int         VertexConnectivity(int vertexindex);
		void        VerticalSegmentIndices(int** pindices,int* pnumseg){_error_("not implemented yet");};
		void        VerticalSegmentIndicesBase(int** pindices,int* pnumseg){_error_("not implemented yet");};
		void			WriteLevelsetSegment(DataSet* segments);
		void        ZeroLevelsetCoordinates(IssmDouble** pxyz_zero,IssmDouble* xyz_list,int levelsetenum);

		#ifdef _HAVE_GIAIVINS_
		void   GiaDeflection(Vector<IssmDouble>* wg,Vector<IssmDouble>* dwgdt,IssmDouble* x,IssmDouble* y);
		#endif
		#ifdef _HAVE_ESA_
		void    EsaGeodetic2D(Vector<IssmDouble>* pUp,Vector<IssmDouble>* pNorth,Vector<IssmDouble>* pEast,IssmDouble* xx,IssmDouble* yy);
		void    EsaGeodetic3D(Vector<IssmDouble>* pUp,Vector<IssmDouble>* pNorth,Vector<IssmDouble>* pEast,IssmDouble* latitude,IssmDouble* longitude,IssmDouble* radius,IssmDouble* xx,IssmDouble* yy,IssmDouble* zz,IssmDouble eartharea);
		#endif
		#ifdef _HAVE_SEALEVELRISE_
		IssmDouble OceanArea(void);
		IssmDouble OceanAverage(IssmDouble* Sg);
		void    SealevelriseMomentOfInertia(IssmDouble* dI_list,IssmDouble* Sg_old,IssmDouble eartharea); 
		void    SealevelriseEustatic(Vector<IssmDouble>* pSgi,IssmDouble* peustatic,IssmDouble* latitude,IssmDouble* longitude,IssmDouble* radius,IssmDouble oceanarea,IssmDouble eartharea);
		void    SealevelriseNonEustatic(Vector<IssmDouble>* pSgo,IssmDouble* Sg_old,IssmDouble* latitude,IssmDouble* longitude,IssmDouble* radius,IssmDouble eartharea);
		void    SealevelriseGeodetic(Vector<IssmDouble>* pUp,Vector<IssmDouble>* pNorth,Vector<IssmDouble>* pEast,IssmDouble* Sg,IssmDouble* latitude,IssmDouble* longitude,IssmDouble* radius,IssmDouble* xx,IssmDouble* yy,IssmDouble* zz,IssmDouble eartharea);
		#endif
		/*}}}*/
		/*Tria specific routines:{{{*/
		void           AddBasalInput(int input_enum, IssmDouble* values, int interpolation_enum);
		void           AddInput(int input_enum, IssmDouble* values, int interpolation_enum);
		IssmDouble     GetArea(void);
		IssmDouble 	   GetArea3D(void);
		IssmDouble 	   GetAreaIce(void);
		IssmDouble 	   GetAreaSpherical(void);
		void           GetAreaCoordinates(IssmDouble *area_coordinates,IssmDouble* xyz_zero,IssmDouble* xyz_list,int numpoints);
		int            GetElementType(void);
		void           GetInputValue(IssmDouble* pvalue,Node* node,int enumtype);
		void		GetLevelsetIntersection(int** pindices, int* pnumiceverts, IssmDouble* fraction, int levelset_enum, IssmDouble level);
		void           GetMaterialInputValue(IssmDouble* pvalue,Node* node,int enumtype);
		Node*          GetNode(int node_number);
		void	         InputUpdateFromSolutionOneDof(IssmDouble* solution,int enum_type);
		void	         InputUpdateFromSolutionOneDofCollapsed(IssmDouble* solution,int enum_type){_error_("not implemented yet");};
		void           JacobianDeterminant(IssmDouble*  pJdet, IssmDouble* xyz_list,Gauss* gauss);
		void           JacobianDeterminantBase(IssmDouble* pJdet,IssmDouble* xyz_list_base,Gauss* gauss);
		void           JacobianDeterminantLine(IssmDouble* Jdet, IssmDouble* xyz_list,Gauss* gauss){_error_("not implemented yet");};
		void           JacobianDeterminantSurface(IssmDouble*  pJdet, IssmDouble* xyz_list,Gauss* gauss);
		void           JacobianDeterminantTop(IssmDouble* pJdet,IssmDouble* xyz_list_base,Gauss* gauss);
		IssmDouble     MinEdgeLength(IssmDouble* xyz_list){_error_("not implemented yet");};
		Gauss*         NewGauss(void);
		Gauss*         NewGauss(int order);
      Gauss*         NewGauss(IssmDouble* xyz_list, IssmDouble* xyz_list_front,int order);
      Gauss*         NewGauss(int point1,IssmDouble fraction1,IssmDouble fraction2,bool mainlyfloating,int order);
      Gauss*         NewGauss(IssmDouble* xyz_list, IssmDouble* xyz_list_front,int order_horiz,int order_vert);
		Gauss*         NewGaussBase(int order);
		Gauss*         NewGaussLine(int vertex1,int vertex2,int order){_error_("not implemented yet");};
		Gauss*         NewGaussTop(int order);
		void           NodalFunctions(IssmDouble* basis,Gauss* gauss);
		void           NodalFunctionsDerivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss);
		void           NodalFunctionsDerivativesVelocity(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss);
		void           NodalFunctionsMINIDerivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss){_error_("not implemented yet");};
		void           NodalFunctionsPressure(IssmDouble* basis,Gauss* gauss);
		void           NodalFunctionsP1(IssmDouble* basis,Gauss* gauss);
		void           NodalFunctionsP1Derivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss);
		void           NodalFunctionsP2(IssmDouble* basis,Gauss* gauss);
		void           NodalFunctionsTensor(IssmDouble* basis,Gauss* gauss);
		void           NodalFunctionsVelocity(IssmDouble* basis,Gauss* gauss);
		void           NormalBase(IssmDouble* normal,IssmDouble* xyz_list);
		void           NormalSection(IssmDouble* normal,IssmDouble* xyz_list);
		void           NormalTop(IssmDouble* normal,IssmDouble* xyz_list);
		void	         SetClone(int* minranks);
		void           SetTemporaryElementType(int element_type_in){_error_("not implemented yet");};
		Seg*	         SpawnSeg(int index1,int index2);
		IssmDouble     StabilizationParameter(IssmDouble u, IssmDouble v, IssmDouble w, IssmDouble diameter, IssmDouble kappa){_error_("not implemented yet");};
		void           UpdateConstraintsExtrudeFromBase(void);
		void           UpdateConstraintsExtrudeFromTop(void);
		void           ViscousHeating(IssmDouble* pphi,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vz_input){_error_("not implemented yet");};
		/*}}}*/

};
#endif  /* _TRIA_H */
