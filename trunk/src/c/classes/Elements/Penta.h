/*! \file Penta.h 
 *  \brief: header file for penta object
 */

#ifndef _PENTA_H_
#define _PENTA_H_

/*Headers:*/
/*{{{*/
#include "./Element.h"
#include "./ElementHook.h"
#include "./PentaRef.h"
class Object;
class Parameters;
class Results;
class Inputs;
class Input;
class IoModel;
class Node;
class Material;
class Matpar;
class Tria;
class ElementMatrix;
class ElementVector;
class GaussPenta;
#include "../../shared/Exceptions/exceptions.h"
#include "../../shared/Enum/Enum.h"
/*}}}*/

class Penta: public Element,public ElementHook,public PentaRef{

	public:

		Penta      **verticalneighbors;           // 2 neighbors: first one under, second one above

		/*Penta constructors and destructor: {{{*/
		Penta(){};
		Penta(int penta_id,int penta_sid,int i, IoModel* iomodel,int nummodels);
		~Penta();
		/*}}}*/
		/*Object virtual functions definitions: {{{*/
		Object *copy();
		void    Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction);
		int     ObjectEnum();
		/*}}}*/
		/*Penta routines:{{{*/
		void           AddBasalInput(int input_enum, IssmDouble* values, int interpolation_enum);
		void           AddInput(int input_enum, IssmDouble* values, int interpolation_enum);
		void           AverageOntoPartition(Vector<IssmDouble>* partition_contributions,Vector<IssmDouble>* partition_areas,IssmDouble* vertex_response,IssmDouble* qmu_part);
		void           BasalNodeIndices(int* pnumindices,int** pindices,int finiteelement);
		void           CalvingRateDev();
		void           CalvingRateLevermann();
		IssmDouble     CharacteristicLength(void){_error_("not implemented yet");};
		void           ComputeBasalStress(Vector<IssmDouble>* sigma_b);
		void           ComputeDeviatoricStressTensor();
		void           ComputeEsaStrainAndVorticity(){_error_("not implemented yet!");};
		void           ComputeSigmaNN(){_error_("not implemented yet");};
		void           ComputeStressTensor();
		void           Configure(Elements* elements,Loads* loads,Nodes* nodes,Vertices* vertices,Materials* materials,Parameters* parameters);
		void           ControlInputSetGradient(IssmDouble* gradient,int enum_type,int control_index);
		void           ControlToVectors(Vector<IssmPDouble>* vector_control, Vector<IssmPDouble>* vector_gradient,int control_enum);
		ElementMatrix* CreateBasalMassMatrix(void);
		void           ElementResponse(IssmDouble* presponse,int response_enum);
		void           ElementSizes(IssmDouble* hx,IssmDouble* hy,IssmDouble* hz);
		int            FiniteElement(void);
		IssmDouble     FloatingArea(void);
		void           FSContactMigration(Vector<IssmDouble>* vertexgrounded,Vector<IssmDouble>* vertexfloating);
		IssmDouble     GetArea3D(void){_error_("not implemented yet!");};
		IssmDouble     GetAreaSpherical(void){_error_("not implemented yet!");};
		void           GetAreaCoordinates(IssmDouble *area_coordinates,IssmDouble* xyz_zero,IssmDouble* xyz_list,int numpoints);
		Element*       GetBasalElement(void);
		Penta*         GetBasalPenta(void);
		int            GetElementType(void);
		void           GetGroundedPart(int* point1,IssmDouble* fraction1, IssmDouble* fraction2,bool* mainlyfloating);
		IssmDouble     GetGroundedPortion(IssmDouble* xyz_list);
		void           GetIcefrontCoordinates(IssmDouble** pxyz_front,IssmDouble* xyz_list,int levelsetenum);
		void           GetInputValue(IssmDouble* pvalue,Node* node,int enumtype);
		void           GetLevelCoordinates(IssmDouble** pxyz_front,IssmDouble* xyz_list,int levelsetenum,IssmDouble level){_error_("not implemented yet");};
		void           GetLevelsetPositivePart(int* point1,IssmDouble* fraction1,IssmDouble* fraction2, bool* mainlynegative,IssmDouble* levelsetvalues){_error_("not implemented yet");};
		Node*          GetNode(int node_number);
		int            GetNodeIndex(Node* node);
		int            GetNumberOfNodes(void);
		int            GetNumberOfNodes(int enum_type);
		int            GetNumberOfVertices(void);
		Penta*         GetLowerPenta(void);
		void           GetSolutionFromInputsOneDof(Vector<IssmDouble>* solution,int enum_type);
		Penta*         GetSurfacePenta(void);
		Element*       GetUpperElement(void);
		Penta*         GetUpperPenta(void);
		void           GetVectorFromControlInputs(Vector<IssmDouble>* gradient,int control_enum,int control_index,const char* data,bool onsid);
		void           GetVerticesCoordinatesBase(IssmDouble** pxyz_list);
		void           GetVerticesCoordinatesTop(IssmDouble** pxyz_list);
		IssmDouble     GroundedArea(void);
		IssmDouble     IceMass(void);
		IssmDouble     IceVolume(void);
		IssmDouble     IceVolumeAboveFloatation(void);
		void           InputControlUpdate(IssmDouble scalar,bool save_parameter);
		void           InputDepthAverageAtBase(int enum_type,int average_enum_type);
		void	         InputExtrude(int enum_type,int start);
		void           InputScale(int enum_type,IssmDouble scale_factor);
		void           InputUpdateFromIoModel(int index, IoModel* iomodel);
		void           InputUpdateFromSolutionOneDof(IssmDouble* solutiong,int enum_type);
		void           InputUpdateFromSolutionOneDofCollapsed(IssmDouble* solutiong,int enum_type);
		void           InputUpdateFromVector(IssmDouble* vector, int name, int type);
		bool           IsFaceOnBoundary(void){_error_("not implemented yet");};
		bool           IsIcefront(void);
		bool           IsNodeOnShelfFromFlags(IssmDouble* flags);
		bool	         IsOnBase(void);
		bool	         IsOnSurface(void);
		bool           IsZeroLevelset(int levelset_enum);
		void           JacobianDeterminant(IssmDouble*  Jdet, IssmDouble* xyz_list,Gauss* gauss);
		void           JacobianDeterminantBase(IssmDouble* pJdet,IssmDouble* xyz_list_base,Gauss* gauss);
		void           JacobianDeterminantLine(IssmDouble* Jdet, IssmDouble* xyz_list,Gauss* gauss);
		void           JacobianDeterminantSurface(IssmDouble*  pJdet, IssmDouble* xyz_list,Gauss* gauss);
		void           JacobianDeterminantTop(IssmDouble* pJdet,IssmDouble* xyz_list_base,Gauss* gauss);
		IssmDouble     Masscon(IssmDouble* levelset){_error_("not implemented yet");};
		IssmDouble     MassFlux(IssmDouble* segment);
		IssmDouble     MassFlux(IssmDouble x1,IssmDouble y1, IssmDouble x2, IssmDouble y2,int segment_id);
		IssmDouble     MinEdgeLength(IssmDouble* xyz_list);
		IssmDouble     Misfit(int modelenum,int observationenum,int weightsenum){_error_("not implemented yet");};
		IssmDouble     MisfitArea(int weightsenum){_error_("not implemented yet");};
		Gauss*         NewGauss(void);
		Gauss*         NewGauss(int order);
		Gauss*         NewGauss(IssmDouble* xyz_list, IssmDouble* xyz_list_front,int order){_error_("not implemented yet");};
		Gauss*         NewGauss(IssmDouble* xyz_list, IssmDouble* xyz_list_front,int order_horiz,int order_vert);
		Gauss*         NewGauss(int point1,IssmDouble fraction1,IssmDouble fraction2,bool mainlyfloating,int order);
		Gauss*         NewGaussBase(int order);
		Gauss*         NewGaussLine(int vertex1,int vertex2,int order);
		Gauss*         NewGaussTop(int order);
		void           NodalFunctions(IssmDouble* basis,Gauss* gauss);
		void           NodalFunctionsDerivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss);
		void           NodalFunctionsDerivativesVelocity(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss);
		void           NodalFunctionsMINIDerivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss);
		void           NodalFunctionsPressure(IssmDouble* basis,Gauss* gauss);
		void           NodalFunctionsP1(IssmDouble* basis,Gauss* gauss);
		void           NodalFunctionsP1Derivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss);
		void           NodalFunctionsP2(IssmDouble* basis,Gauss* gauss);
		void           NodalFunctionsTensor(IssmDouble* basis,Gauss* gauss);
		void           NodalFunctionsVelocity(IssmDouble* basis,Gauss* gauss);
		void	         NormalBase(IssmDouble* bed_normal, IssmDouble* xyz_list);
		void           NormalSection(IssmDouble* normal,IssmDouble* xyz_list);
		void	         NormalTop(IssmDouble* bed_normal, IssmDouble* xyz_list);
		int            NodalValue(IssmDouble* pvalue, int index, int natureofdataenum);
		int            NumberofNodesPressure(void);
		int            NumberofNodesVelocity(void);
		void           PotentialUngrounding(Vector<IssmDouble>* potential_sheet_ungrounding);
		int            PressureInterpolation();
		void           ReduceMatrices(ElementMatrix* Ke,ElementVector* pe);
		void           ResetFSBasalBoundaryCondition(void);
		void           ResetHooks();
		void           ResetLevelsetFromSegmentlist(IssmDouble* segments,int numsegments);
		void	         SetClone(int* minranks);
		void           SetControlInputsFromVector(IssmDouble* vector,int control_enum,int control_index);
		void           SetCurrentConfiguration(Elements* elements,Loads* loads,Nodes* nodes,Materials* materials,Parameters* parameters);
		void           SetTemporaryElementType(int element_type_in);
	   Element*       SpawnBasalElement(void);
		Element*       SpawnTopElement(void);
		Tria*	         SpawnTria(int index1,int index2,int index3);
		IssmDouble     StabilizationParameter(IssmDouble u, IssmDouble v, IssmDouble w, IssmDouble diameter, IssmDouble kappa);
		void           StressIntensityFactor();
		void           StrainRateparallel();
		void           StrainRateperpendicular();
		IssmDouble     SurfaceArea(void);
		int            TensorInterpolation(){_error_("not implemented yet");};
		IssmDouble     TimeAdapt();
		IssmDouble     TotalFloatingBmb(void);
		IssmDouble     TotalGroundedBmb(void);
		IssmDouble     TotalSmb(void);
		void           Update(int index, IoModel* iomodel,int analysis_counter,int analysis_type,int finitelement);
		void           UpdateConstraintsExtrudeFromBase(void);
		void           UpdateConstraintsExtrudeFromTop(void);
		int            UpdatePotentialUngrounding(IssmDouble* potential_sheet_ungrounding,Vector<IssmDouble>* vec_nodes_on_iceshelf,IssmDouble* nodes_on_iceshelf);
		void           ValueP1DerivativesOnGauss(IssmDouble* dvalue,IssmDouble* values,IssmDouble* xyz_list,Gauss* gauss);
		void           ValueP1OnGauss(IssmDouble* pvalue,IssmDouble* values,Gauss* gauss);
		int            VelocityInterpolation();
		int            VertexConnectivity(int vertexindex);
		void           VerticalSegmentIndices(int** pindices,int* pnumseg);
		void           VerticalSegmentIndicesBase(int** pindices,int* pnumseg);
		void           ViscousHeating(IssmDouble* pphi,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vz_input);
		void           ZeroLevelsetCoordinates(IssmDouble** pxyz_zero,IssmDouble* xyz_list,int levelsetenum);

		#ifdef _HAVE_DAKOTA_
		void           InputUpdateFromMatrixDakota(IssmDouble* matrix, int nows, int ncols, int name, int type);
		void           InputUpdateFromVectorDakota(IssmDouble* vector, int name, int type);
		#endif

		#ifdef _HAVE_GIAIVINS_
		void           GiaDeflection(Vector<IssmDouble>* wg,Vector<IssmDouble>* dwgdt,IssmDouble* x,IssmDouble* y);
		#endif
		#ifdef _HAVE_ESA_
		void    EsaGeodetic2D(Vector<IssmDouble>* pUp,Vector<IssmDouble>* pNorth,Vector<IssmDouble>* pEast,IssmDouble* xx,IssmDouble* yy){_error_("not implemented yet!");};
		void    EsaGeodetic3D(Vector<IssmDouble>* pUp,Vector<IssmDouble>* pNorth,Vector<IssmDouble>* pEast,IssmDouble* latitude,IssmDouble* longitude,IssmDouble* radius,IssmDouble* xx,IssmDouble* yy,IssmDouble* zz,IssmDouble eartharea){_error_("not implemented yet!");};
		#endif
		#ifdef _HAVE_SEALEVELRISE_
		IssmDouble    OceanArea(void){_error_("not implemented yet!");};
		IssmDouble    OceanAverage(IssmDouble* Sg){_error_("not implemented yet!");};
		void    SealevelriseMomentOfInertia(IssmDouble* dI_list,IssmDouble* Sg_old,IssmDouble eartharea){_error_("not implemented yet!");};
		void    SealevelriseEustatic(Vector<IssmDouble>* pSgi,IssmDouble* peustatic,IssmDouble* latitude,IssmDouble* longitude,IssmDouble* radius,IssmDouble oceanarea,IssmDouble eartharea){_error_("not implemented yet!");};
		void    SealevelriseNonEustatic(Vector<IssmDouble>* pSgo,IssmDouble* Sg_old,IssmDouble* latitude,IssmDouble* longitude,IssmDouble* radius,IssmDouble eartharea){_error_("not implemented yet!");};
		void    SealevelriseGeodetic(Vector<IssmDouble>* pUp,Vector<IssmDouble>* pNorth,Vector<IssmDouble>* pEast,IssmDouble* Sg,IssmDouble* latitude,IssmDouble* longitude,IssmDouble* radius,IssmDouble* xx,IssmDouble* yy,IssmDouble* zz,IssmDouble eartharea){_error_("not implemented yet!");};
		#endif

		/*}}}*/
};
#endif  /* _PENTA_H */
