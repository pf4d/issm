/*! \file Seg.h 
 *  \brief: header file for seg object
 */

#ifndef _SEG_H_
#define _SEG_H_

/*Headers:*/
/*{{{*/
#include "./Element.h"
#include "./ElementHook.h"
#include "./SegRef.h"
class Parameters;
class Inputs;
class IoModel;
class Results;
class Node;
class Material;
class Matpar;
class ElementMatrix;
class ElementVector;
class Vertex;

#include "../../shared/Exceptions/exceptions.h"
#include "../../shared/Enum/Enum.h"
/*}}}*/

class Seg: public Element,public ElementHook,public SegRef{

	public:

		/*Seg constructors, destructors {{{*/
		Seg(){};
		Seg(int seg_id,int seg_sid,int i, IoModel* iomodel,int nummodels);
		~Seg();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Object *copy();
		void    Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction);
		int     ObjectEnum();
		/*}}}*/
		/*Element virtual functions definitions: {{{*/
		void        AddBasalInput(int input_enum, IssmDouble* values, int interpolation_enum){_error_("not implemented yet");};
		void        AddInput(int input_enum, IssmDouble* values, int interpolation_enum){_error_("not implemented yet");};
		void        AverageOntoPartition(Vector<IssmDouble>* partition_contributions,Vector<IssmDouble>* partition_areas,IssmDouble* vertex_response,IssmDouble* qmu_part){_error_("not implemented yet");};
		void        CalvingRateLevermann(void){_error_("not implemented yet");};
		IssmDouble  CharacteristicLength(void);
		void        ComputeBasalStress(Vector<IssmDouble>* sigma_b){_error_("not implemented yet");};
		void        ComputeDeviatoricStressTensor(){_error_("not implemented yet");};
		void        ComputeEsaStrainAndVorticity(){_error_("not implemented yet!");};
		void        ComputeSigmaNN(){_error_("not implemented yet");};
		void        ComputeStressTensor(){_error_("not implemented yet");};
		void        Configure(Elements* elements,Loads* loads,Nodes* nodesin,Vertices* verticesin,Materials* materials,Parameters* parameters){_error_("not implemented yet");};
		void        ControlInputSetGradient(IssmDouble* gradient,int enum_type,int control_index){_error_("not implemented yet");};
		void        ControlToVectors(Vector<IssmPDouble>* vector_control, Vector<IssmPDouble>* vector_gradient,int control_enum){_error_("not implemented yet");};
		void        ElementResponse(IssmDouble* presponse,int response_enum){_error_("not implemented yet");};
		void        ElementSizes(IssmDouble* hx,IssmDouble* hy,IssmDouble* hz){_error_("not implemented yet");};
		int         FiniteElement(void);
		IssmDouble  FloatingArea(void){_error_("not implemented yet");};
		void        FSContactMigration(Vector<IssmDouble>* vertexgrounded,Vector<IssmDouble>* vertexfloating){_error_("not implemented yet");};
		Element*    GetBasalElement(void){_error_("not implemented yet");};
		int         GetElementType(void){_error_("not implemented yet");};
		void        GetGroundedPart(int* point1,IssmDouble* fraction1, IssmDouble* fraction2,bool* mainlyfloating){_error_("not implemented yet");};
		IssmDouble  GetGroundedPortion(IssmDouble* xyz_list){_error_("not implemented yet");};
		void		   GetIcefrontCoordinates(IssmDouble** pxyz_front,IssmDouble* xyz_list,int levelsetenum);
		void        GetInputValue(IssmDouble* pvalue,Node* node,int enumtype){_error_("not implemented yet");};
		void		   GetLevelCoordinates(IssmDouble** pxyz_front,IssmDouble* xyz_list,int levelsetenum,IssmDouble level){_error_("not implemented");};
		void        GetLevelsetPositivePart(int* point1,IssmDouble* fraction1,IssmDouble* fraction2, bool* mainlynegative,IssmDouble* levelsetvalues){_error_("not implemented yet");};
		Node*       GetNode(int node_number){_error_("Not implemented");};
		int         GetNodeIndex(Node* node){_error_("not implemented yet");};
		int         GetNumberOfNodes(void);
		int         GetNumberOfNodes(int enum_type){_error_("not implemented yet");};
		int         GetNumberOfVertices(void);
		void        GetSolutionFromInputsOneDof(Vector<IssmDouble>* solution,int enum_type){_error_("not implemented yet");};
		Element*    GetUpperElement(void){_error_("not implemented yet");};
		void        GetVectorFromControlInputs(Vector<IssmDouble>* gradient,int control_enum,int control_index,const char* data,bool onsid){_error_("not implemented yet");};
		void        GetVerticesCoordinates(IssmDouble** pxyz_list);
		void        GetVerticesCoordinatesBase(IssmDouble** pxyz_list){_error_("not implemented yet");};
		void        GetVerticesCoordinatesTop(IssmDouble** pxyz_list){_error_("not implemented yet");};
		IssmDouble  GroundedArea(void){_error_("not implemented yet");};
		IssmDouble  IceMass(void){_error_("not implemented yet");};
		IssmDouble  IceVolume(void){_error_("not implemented yet");};
		IssmDouble  IceVolumeAboveFloatation(void){_error_("not implemented yet");};
		void        InputControlUpdate(IssmDouble scalar,bool save_parameter){_error_("not implemented yet");};
		void        InputDepthAverageAtBase(int enum_type,int average_enum_type){_error_("not implemented yet");};
		void        InputExtrude(int enum_type,int start){_error_("not implemented"); /*For penta only*/};
		void        InputScale(int enum_type,IssmDouble scale_factor){_error_("not implemented yet");};
		void        InputUpdateFromIoModel(int index, IoModel* iomodel){_error_("not implemented yet");};
		void        InputUpdateFromSolutionOneDof(IssmDouble* solution,int inputenum){_error_("not implemented yet");};
		void        InputUpdateFromSolutionOneDofCollapsed(IssmDouble* solution,int inputenum){_error_("not implemented yet");};
		void        InputUpdateFromVector(IssmDouble* vector, int name, int type){_error_("not implemented yet");};
		bool        IsFaceOnBoundary(void){_error_("not implemented yet");};
		bool		   IsIcefront(void);
		bool        IsNodeOnShelfFromFlags(IssmDouble* flags){_error_("not implemented yet");};
		bool        IsOnBase(){_error_("not implemented yet");};
		bool        IsOnSurface(){_error_("not implemented yet");};
		bool        IsZeroLevelset(int levelset_enum){_error_("not implemented");};
		void        JacobianDeterminant(IssmDouble*  Jdet, IssmDouble* xyz_list,Gauss* gauss);
		void        JacobianDeterminantBase(IssmDouble* pJdet,IssmDouble* xyz_list_base,Gauss* gauss){_error_("not implemented yet");};
		void        JacobianDeterminantLine(IssmDouble* Jdet, IssmDouble* xyz_list,Gauss* gauss){_error_("not implemented yet");};
		void        JacobianDeterminantSurface(IssmDouble*  pJdet, IssmDouble* xyz_list,Gauss* gauss);
		void        JacobianDeterminantTop(IssmDouble* pJdet,IssmDouble* xyz_list_base,Gauss* gauss){_error_("not implemented yet");};
		IssmDouble  Masscon(IssmDouble* levelset){_error_("not implemented yet");};
		IssmDouble  MassFlux(IssmDouble* segment){_error_("not implemented yet");};
		IssmDouble  MassFlux(IssmDouble x1,IssmDouble y1, IssmDouble x2, IssmDouble y2,int segment_id){_error_("not implemented yet");}
		void        MaterialUpdateFromTemperature(void){_error_("not implemented yet");};
		IssmDouble  MinEdgeLength(IssmDouble* xyz_list){_error_("not implemented yet");};
		IssmDouble  Misfit(int modelenum,int observationenum,int weightsenum){_error_("not implemented yet");};
		IssmDouble  MisfitArea(int weightsenum){_error_("not implemented yet");};
		Gauss*      NewGauss(void);
		Gauss*      NewGauss(int order);
      Gauss*      NewGauss(IssmDouble* xyz_list, IssmDouble* xyz_list_front,int order);
      Gauss*      NewGauss(IssmDouble* xyz_list, IssmDouble* xyz_list_front,int order_horiz,int order_vert){_error_("not implemented yet");};
      Gauss*      NewGauss(int point1,IssmDouble fraction1,IssmDouble fraction2,bool mainlyfloating,int order){_error_("not implemented yet");};
		Gauss*      NewGaussBase(int order){_error_("not implemented yet");};
		Gauss*      NewGaussLine(int vertex1,int vertex2,int order){_error_("not implemented yet");};
		Gauss*      NewGaussTop(int order){_error_("not implemented yet");};
		int         NodalValue(IssmDouble* pvalue, int index, int natureofdataenum){_error_("not implemented yet");};
		void        NodalFunctions(IssmDouble* basis,Gauss* gauss);
		void        NodalFunctionsDerivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss);
		void        NodalFunctionsDerivativesVelocity(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss){_error_("not implemented yet");};
		void        NodalFunctionsMINIDerivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss){_error_("not implemented yet");};
		void        NodalFunctionsPressure(IssmDouble* basis,Gauss* gauss){_error_("not implemented yet");};
		void        NodalFunctionsP1(IssmDouble* basis,Gauss* gauss);
		void        NodalFunctionsP1Derivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss){_error_("not implemented yet");};
		void        NodalFunctionsP2(IssmDouble* basis,Gauss* gauss);
		void        NodalFunctionsTensor(IssmDouble* basis,Gauss* gauss){_error_("not implemented yet");};
		void        NodalFunctionsVelocity(IssmDouble* basis,Gauss* gauss){_error_("not implemented yet");};
		void        NormalBase(IssmDouble* normal,IssmDouble* xyz_list){_error_("not implemented yet");};
		void        NormalSection(IssmDouble* normal,IssmDouble* xyz_list);
		void        NormalTop(IssmDouble* normal,IssmDouble* xyz_list){_error_("not implemented yet");};
		int         NumberofNodesPressure(void){_error_("not implemented yet");};
		int         NumberofNodesVelocity(void){_error_("not implemented yet");};
		void        PotentialUngrounding(Vector<IssmDouble>* potential_sheet_ungrounding){_error_("not implemented yet");};
		int         PressureInterpolation(void){_error_("not implemented yet");};
		void        ReduceMatrices(ElementMatrix* Ke,ElementVector* pe){_error_("not implemented yet");};
		void        ResetFSBasalBoundaryCondition(void){_error_("not implemented yet");};
		void        ResetHooks(){_error_("not implemented yet");};
		void        SetControlInputsFromVector(IssmDouble* vector,int control_enum,int control_index){_error_("not implemented yet");};
		void        SetCurrentConfiguration(Elements* elements,Loads* loads,Nodes* nodes,Materials* materials,Parameters* parameters){_error_("not implemented yet");};
		void        SetTemporaryElementType(int element_type_in){_error_("not implemented yet");};
	   Element*    SpawnBasalElement(void){_error_("not implemented yet");};
		Element*    SpawnTopElement(void){_error_("not implemented yet");};
		IssmDouble  StabilizationParameter(IssmDouble u, IssmDouble v, IssmDouble w, IssmDouble diameter, IssmDouble kappa){_error_("not implemented yet");};
		void        StrainRateparallel(void){_error_("not implemented yet");};
		void        StrainRateperpendicular(void){_error_("not implemented yet");};
		void        StressIntensityFactor(void){_error_("not implemented yet");};
		IssmDouble  SurfaceArea(void){_error_("not implemented yet");};
		int         TensorInterpolation(void){_error_("not implemented yet");};
		IssmDouble  TimeAdapt(){_error_("not implemented yet");};
		IssmDouble  TotalFloatingBmb(void){_error_("not implemented yet");};
		IssmDouble  TotalGroundedBmb(void){_error_("not implemented yet");};
		IssmDouble  TotalSmb(void){_error_("not implemented yet");};
		void        Update(int index, IoModel* iomodel,int analysis_counter,int analysis_type,int finitelement){_error_("not implemented yet");};
		void        UpdateConstraintsExtrudeFromBase(){_error_("not implemented");};
		void        UpdateConstraintsExtrudeFromTop(){_error_("not implemented");};
		int         UpdatePotentialUngrounding(IssmDouble* vertices_potentially_ungrounding,Vector<IssmDouble>* vec_nodes_on_iceshelf,IssmDouble* nodes_on_iceshelf){_error_("not implemented yet");};
		void        ValueP1DerivativesOnGauss(IssmDouble* dvalue,IssmDouble* values,IssmDouble* xyz_list,Gauss* gauss){_error_("not implemented yet");};
		void        ValueP1OnGauss(IssmDouble* pvalue,IssmDouble* values,Gauss* gauss){_error_("not implemented yet");};
		int         VelocityInterpolation(void){_error_("not implemented yet");};
		int         VertexConnectivity(int vertexindex){_error_("not implemented yet");};
		void        VerticalSegmentIndices(int** pindices,int* pnumseg){_error_("not implemented yet");};
		void        VerticalSegmentIndicesBase(int** pindices,int* pnumseg){_error_("not implemented yet");};
		void        ViscousHeating(IssmDouble* pphi,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vz_input){_error_("not implemented yet");};
		void        ZeroLevelsetCoordinates(IssmDouble** pxyz_zero,IssmDouble* xyz_list,int levelsetenum){_error_("not implemented");};
		IssmDouble     GetArea3D(void){_error_("not implemented yet!");};
		IssmDouble     GetAreaSpherical(void){_error_("not implemented yet!");};

#ifdef _HAVE_GIAIVINS_
		void        GiaDeflection(Vector<IssmDouble>* wg,Vector<IssmDouble>* dwgdt,IssmDouble* x,IssmDouble* y){_error_("not implemented yet");};
#endif
#ifdef _HAVE_ESA_
		void    EsaGeodetic2D(Vector<IssmDouble>* pUp,Vector<IssmDouble>* pNorth,Vector<IssmDouble>* pEast,IssmDouble* xx,IssmDouble* yy){_error_("not implemented yet!");};
		void    EsaGeodetic3D(Vector<IssmDouble>* pUp,Vector<IssmDouble>* pNorth,Vector<IssmDouble>* pEast,IssmDouble* latitude,IssmDouble* longitude,IssmDouble* radius,IssmDouble* xx,IssmDouble* yy,IssmDouble* zz,IssmDouble eartharea){_error_("not implemented yet!");};
#endif
#ifdef _HAVE_SEALEVELRISE_
		void    SealevelriseMomentOfInertia(IssmDouble* dI_list,IssmDouble* Sg_old,IssmDouble eartharea){_error_("not implemented yet!");};
		void    SealevelriseEustatic(Vector<IssmDouble>* pSgi,IssmDouble* peustatic,IssmDouble* latitude,IssmDouble* longitude,IssmDouble* radius,IssmDouble oceanarea,IssmDouble eartharea){_error_("not implemented yet!");};
		void    SealevelriseNonEustatic(Vector<IssmDouble>* pSgo,IssmDouble* Sg_old,IssmDouble* latitude,IssmDouble* longitude,IssmDouble* radius,IssmDouble eartharea){_error_("not implemented yet!");};
		void    SealevelriseGeodetic(Vector<IssmDouble>* pUp,Vector<IssmDouble>* pNorth,Vector<IssmDouble>* pEast,IssmDouble* Sg,IssmDouble* latitude,IssmDouble* longitude,IssmDouble* radius,IssmDouble* xx,IssmDouble* yy,IssmDouble* zz,IssmDouble eartharea){_error_("not implemented yet!");};
		IssmDouble    OceanArea(void){_error_("not implemented yet!");};
		IssmDouble    OceanAverage(IssmDouble* Sg){_error_("not implemented yet!");};
#endif

#ifdef _HAVE_DAKOTA_
		void        InputUpdateFromMatrixDakota(IssmDouble* matrix, int nows, int ncols, int name, int type){_error_("not implemented yet");};
		void        InputUpdateFromVectorDakota(IssmDouble* vector, int name, int type){_error_("not implemented yet");};
#endif
		/*}}}*/
};
#endif  /* _SEG_H */
