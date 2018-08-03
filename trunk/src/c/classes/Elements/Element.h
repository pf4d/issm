/*!\file:  Element.h
 * \brief abstract class for Element object
 * This class is a place holder for the Tria and the Penta elements. 
 * It is derived from Element, so DataSets can contain them.
 */ 

#ifndef _ELEMENT_H_
#define _ELEMENT_H_

/*Headers:*/
/*{{{*/
#include "../../datastructures/datastructures.h"
#include "../../toolkits/toolkits.h"
#include "../Update.h"
class DataSet;
class Parameters;
class Parameter;
class Elements;
class Loads;
class Nodes;
class Node;
class Vertices;
class Vertex;
class Materials;
class Material;
class Matpar;
class Inputs;
class Input;
class Gauss;
class ElementVector;
template <class doublematrix> class Matrix;
template <class doubletype> class Vector;
class ElementMatrix;
class ElementVector;
/*}}}*/

class Element: public Object,public Update{

	public:
		int          id;
		int          sid;
		Inputs      *inputs;
		Node       **nodes;
		Vertex     **vertices;
		Material    *material;
		Matpar      *matpar;
		Parameters  *parameters;

		int* element_type_list;
		int  element_type;

	public: 
		/*Constructors/Destructores*/
		Element();
		~Element();

		/*Functions*/
		void               AddInput(Input* input_in);
		/*bool               AllActive(void);*/
		/*bool               AnyActive(void);*/
		void               ComputeLambdaS(void);
		void               ComputeNewDamage();
		void               ComputeStrainRate();
		void               CoordinateSystemTransform(IssmDouble** ptransform,Node** nodes,int numnodes,int* cs_array);
		void               DeepEcho();
		void               DeleteInput(int input_enum);
		void               DeleteMaterials(void);
		void               Delta18oParameterization(void);
		void               Delta18opdParameterization(void);
		IssmDouble         Divergence(void);
		void               dViscositydBFS(IssmDouble* pdmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vz_input);
		void               dViscositydBHO(IssmDouble* pdmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input);
		void               dViscositydBSSA(IssmDouble* pdmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input);
		void               dViscositydDSSA(IssmDouble* pdmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input);
		void               Echo();
		IssmDouble         EnthalpyDiffusionParameter(IssmDouble enthalpy,IssmDouble pressure);
		IssmDouble         EnthalpyDiffusionParameterVolume(int numvertices,IssmDouble* enthalpy,IssmDouble* pressure);
		void               EnthalpyToThermal(IssmDouble* ptemperature,IssmDouble* pwaterfraction,IssmDouble enthalpy,IssmDouble pressure);
		void               FindParam(bool* pvalue,int paramenum);
		void               FindParam(int* pvalue,int paramenum);
		void               FindParam(IssmDouble* pvalue,int paramenum);
		void               FindParam(int** pvalues,int* psize,int paramenum);
		void	             GetDofList(int** pdoflist,int approximation_enum,int setenum);
		void	             GetDofListPressure(int** pdoflist,int setenum);
		void	             GetDofListVelocity(int** pdoflist,int setenum);
		Input*             GetInput(int inputenum);
		void               GetInputListOnNodes(IssmDouble* pvalue,int enumtype);
		void               GetInputListOnNodes(IssmDouble* pvalue,int enumtype,IssmDouble defaultvalue);
		void               GetInputListOnNodesVelocity(IssmDouble* pvalue,int enumtype);
		void               GetInputListOnVertices(IssmDouble* pvalue,int enumtype);
		void               GetInputListOnVertices(IssmDouble* pvalue,int enumtype,IssmDouble defaultvalue);
		void               GetInputLocalMinMaxOnNodes(IssmDouble* min,IssmDouble* max,IssmDouble* ug);
		void               GetInputValue(bool* pvalue,int enum_type);
		void               GetInputValue(int* pvalue,int enum_type);
		void               GetInputValue(IssmDouble* pvalue,int enum_type);
		void               GetInputValue(IssmDouble* pvalue,Gauss* gauss,int enum_type);
		void               GetInputsInterpolations(Vector<IssmDouble>* interps);
		IssmDouble         GetMaterialParameter(int enum_in);
		void               GetNodesLidList(int* lidlist);
		void               GetNodesSidList(int* sidlist);
		void               GetPhi(IssmDouble* phi, IssmDouble*  epsilon, IssmDouble viscosity);
		void               GetVectorFromInputs(Vector<IssmDouble>* vector, int name_enum, int type);
		void	             GetVertexPidList(int* pidlist);
		void               GetVerticesConnectivityList(int* connectivitylist);
		void               GetVerticesCoordinates(IssmDouble** xyz_list);
		void               GetVerticesSidList(int* sidlist);
		IssmDouble         GetXcoord(IssmDouble* xyz_list,Gauss* gauss);
		IssmDouble         GetYcoord(IssmDouble* xyz_list,Gauss* gauss);
		IssmDouble         GetZcoord(IssmDouble* xyz_list,Gauss* gauss);
		void               GradientIndexing(int* indexing,int control_index,bool onsid=false);
		bool               HasNodeOnBase();
		bool               HasNodeOnSurface();
		int                Id();
		void               InputChangeName(int enum_type,int enum_type_old);
		void               InputCreate(IssmDouble* vector,IoModel* iomodel,int M,int N,int vector_type,int vector_enum,int code);
		void               InputDuplicate(int original_enum,int new_enum);
		void               InputUpdateFromConstant(IssmDouble constant, int name);
		void               InputUpdateFromConstant(int constant, int name);
		void               InputUpdateFromConstant(bool constant, int name);
		bool               IsFloating(); 
		bool               IsIceInElement();
		bool	             IsInput(int name);
		bool               IsLandInElement();
		bool               IsWaterInElement();
		void               LinearFloatingiceMeltingRate(); 
		void               MantlePlumeGeothermalFlux(); 
		void               MarshallElement(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction,int numanalyses);
		void               MigrateGroundingLine(IssmDouble* sheet_ungrounding);
		void               MismipFloatingiceMeltingRate(); 
		void               MungsmtpParameterization(void);
		ElementMatrix*     NewElementMatrix(int approximation_enum=NoneApproximationEnum);
		ElementMatrix*     NewElementMatrixCoupling(int number_nodes,int approximation_enum=NoneApproximationEnum);
		ElementVector*     NewElementVector(int approximation_enum=NoneApproximationEnum);
		void               PositiveDegreeDay(IssmDouble* pdds,IssmDouble* pds,IssmDouble signorm,bool ismungsm);
		IssmDouble         PureIceEnthalpy(IssmDouble pressure);
		void               ResultInterpolation(int* pinterpolation,int*nodesperelement,int* parray_size, int output_enum);
		void               ResultToPatch(IssmDouble* values,int nodesperelement,int output_enum);
		void               ResultToMatrix(IssmDouble* values,int ncols,int output_enum);
		void               ResultToVector(Vector<IssmDouble>* vector,int output_enum);
		void               SetwiseNodeConnectivity(int* d_nz,int* o_nz,Node* node,bool* flags,int* flagsindices,int set1_enum,int set2_enum);
		int                Sid();
		void               SmbGemb();
		void               StrainRateESA(IssmDouble* epsilon,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input);
		void               StrainRateFS(IssmDouble* epsilon,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vz_input);
		void               StrainRateHO(IssmDouble* epsilon,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input);
		void               StrainRateHO2dvertical(IssmDouble* epsilon,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input);
		void               StrainRateSSA(IssmDouble* epsilon,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input);
		void               StrainRateSSA1d(IssmDouble* epsilon,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input);
		void               StressMaxPrincipalCreateInput(void);
		void               ThermalToEnthalpy(IssmDouble* penthalpy,IssmDouble temperature,IssmDouble waterfraction,IssmDouble pressure);
		IssmDouble         TMeltingPoint(IssmDouble pressure);
		void               TransformInvStiffnessMatrixCoord(ElementMatrix* Ke,int cs_enum);
		void               TransformInvStiffnessMatrixCoord(ElementMatrix* Ke,Node** nodes,int numnodes,int cs_enum);
		void               TransformInvStiffnessMatrixCoord(ElementMatrix* Ke,Node** nodes,int numnodes,int* cs_array);
		void               TransformLoadVectorCoord(ElementVector* pe,int cs_enum);
		void               TransformLoadVectorCoord(ElementVector* pe,int* cs_array);
		void               TransformLoadVectorCoord(ElementVector* pe,Node** nodes,int numnodes,int cs_enum);
		void               TransformLoadVectorCoord(ElementVector* pe,Node** nodes,int numnodes,int* cs_array);
		void               TransformLoadVectorCoord(ElementVector* pe,int numnodes,int transformenum){_error_("not implemented yet");};/*Tiling only*/
		void               TransformLoadVectorCoord(ElementVector* pe,int numnodes,int* transformenum_list){_error_("not implemented yet");};/*Tiling only*/
		void               TransformSolutionCoord(IssmDouble* solution,int cs_enum);
		void               TransformSolutionCoord(IssmDouble* solution,int* cs_array);
		void               TransformSolutionCoord(IssmDouble* solution,int numnodes,int cs_enum);
		void               TransformSolutionCoord(IssmDouble* solution,int numnodes,int* cs_array);
		void               TransformSolutionCoord(IssmDouble* solution,Node** nodes,int numnodes,int cs_enum);
		void               TransformSolutionCoord(IssmDouble* solution,Node** nodes,int numnodes,int* cs_array);
		void               TransformStiffnessMatrixCoord(ElementMatrix* Ke,int cs_enum);
		void               TransformStiffnessMatrixCoord(ElementMatrix* Ke,int* cs_array);
		void               TransformStiffnessMatrixCoord(ElementMatrix* Ke,Node** nodes,int numnodes,int cs_enum);
		void               TransformStiffnessMatrixCoord(ElementMatrix* Ke,Node** nodes,int numnodes,int* cs_array);
		void               TransformStiffnessMatrixCoord(ElementMatrix* Ke,int numnodes,int* transformenum_list){_error_("not implemented yet");};/*Tiling only*/
		void               ViscousHeatingCreateInput(void);

		/*Virtual functions*/
		virtual void       AddBasalInput(int input_enum, IssmDouble* values, int interpolation_enum)=0;
		virtual void       AddInput(int input_enum, IssmDouble* values, int interpolation_enum)=0;
		virtual void       AverageOntoPartition(Vector<IssmDouble>* partition_contributions,Vector<IssmDouble>* partition_areas,IssmDouble* vertex_response,IssmDouble* qmu_part)=0;
		virtual void		 BasalNodeIndices(int* pnumindices,int** pindices,int finiteelement){_error_("not implemented yet");};
		virtual void       CalvingRateDev(void){_error_("not implemented yet");};
		virtual void	    CalvingRateLevermann(void)=0;
		virtual IssmDouble CharacteristicLength(void)=0;
		virtual void       ComputeBasalStress(Vector<IssmDouble>* sigma_b)=0;
		virtual void       ComputeDeviatoricStressTensor(void)=0;
		virtual void       ComputeSigmaNN(void)=0;
		virtual void       ComputeStressTensor(void)=0;
		virtual void       ComputeEsaStrainAndVorticity(void)=0;
		virtual void       Configure(Elements* elements,Loads* loads,Nodes* nodes,Vertices* vertices,Materials* materials,Parameters* parameters)=0;
		virtual void       ControlInputSetGradient(IssmDouble* gradient,int enum_type,int control_index)=0;
		virtual void       ControlToVectors(Vector<IssmPDouble>* vector_control, Vector<IssmPDouble>* vector_gradient,int control_enum)=0;
		virtual void       ElementResponse(IssmDouble* presponse,int response_enum)=0;
		virtual void       ElementSizes(IssmDouble* phx,IssmDouble* phy,IssmDouble* phz)=0;
		virtual int        FiniteElement(void)=0;
		virtual IssmDouble FloatingArea(void)=0;
		virtual void       FSContactMigration(Vector<IssmDouble>* vertexgrounded,Vector<IssmDouble>* vertexfloating)=0;
		virtual Element*   GetBasalElement(void)=0;
		virtual int        GetElementType(void)=0;
		virtual void       GetGroundedPart(int* point1,IssmDouble* fraction1,IssmDouble* fraction2, bool* mainlyfloating)=0;
		virtual IssmDouble GetGroundedPortion(IssmDouble* xyz_list)=0;
		virtual void       GetIcefrontCoordinates(IssmDouble** pxyz_front,IssmDouble* xyz_list,int levelsetenum)=0;
		virtual void       GetInputValue(IssmDouble* pvalue,Node* node,int enumtype)=0;
		virtual void       GetLevelCoordinates(IssmDouble** pxyz_front,IssmDouble* xyz_list,int levelsetenum,IssmDouble level)=0;
		virtual void       GetLevelsetPositivePart(int* point1,IssmDouble* fraction1,IssmDouble* fraction2, bool* mainlynegative,IssmDouble* levelsetvalues)=0;
		virtual Node*      GetNode(int node_number)=0;
		virtual int        GetNodeIndex(Node* node)=0;
		virtual int        GetNumberOfNodes(void)=0;
		virtual int        GetNumberOfNodes(int enum_type)=0;
		virtual int        GetNumberOfVertices(void)=0;
		virtual void       GetSolutionFromInputsOneDof(Vector<IssmDouble>* solution,int solutionenum)=0;
		virtual Element*   GetUpperElement(void)=0;
		virtual void       GetVectorFromControlInputs(Vector<IssmDouble>* gradient,int control_enum,int control_index,const char* data,bool onsid)=0;
		virtual void       GetVerticesCoordinatesBase(IssmDouble** xyz_list)=0;
		virtual void       GetVerticesCoordinatesTop(IssmDouble** xyz_list)=0;
		virtual IssmDouble GroundedArea(void)=0;
		virtual IssmDouble IceMass(void)=0;
		virtual IssmDouble IceVolume(void)=0;
		virtual IssmDouble IceVolumeAboveFloatation(void)=0;
		virtual void       InputControlUpdate(IssmDouble scalar,bool save_parameter)=0;
		virtual void       InputDepthAverageAtBase(int enum_type,int average_enum_type)=0;
		virtual void       InputExtrude(int input_enum,int start)=0;
		virtual void       InputScale(int enum_type,IssmDouble scale_factor)=0;
		virtual void       InputUpdateFromSolutionOneDofCollapsed(IssmDouble* solution,int inputenum)=0;
		virtual void       InputUpdateFromSolutionOneDof(IssmDouble* solution,int inputenum)=0;
		virtual bool       IsFaceOnBoundary(void)=0;
		virtual bool       IsIcefront(void)=0;
		virtual bool       IsNodeOnShelfFromFlags(IssmDouble* flags)=0; 
		virtual bool       IsOnBase()=0;
		virtual bool       IsOnSurface()=0;
		virtual bool       IsZeroLevelset(int levelset_enum)=0;
		virtual void       JacobianDeterminant(IssmDouble*  Jdet, IssmDouble* xyz_list,Gauss* gauss)=0;
		virtual void       JacobianDeterminantBase(IssmDouble* Jdet,IssmDouble* xyz_list_base,Gauss* gauss)=0;
		virtual void       JacobianDeterminantLine(IssmDouble* Jdet, IssmDouble* xyz_list,Gauss* gauss)=0;
		virtual void       JacobianDeterminantSurface(IssmDouble* Jdet, IssmDouble* xyz_list,Gauss* gauss)=0;
		virtual void       JacobianDeterminantTop(IssmDouble* Jdet,IssmDouble* xyz_list_base,Gauss* gauss)=0;
		virtual void       Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction)=0;
		virtual IssmDouble Masscon(IssmDouble* levelset)=0;
		virtual IssmDouble MassFlux(IssmDouble* segment)=0;
		virtual IssmDouble MassFlux(IssmDouble x1,IssmDouble y1, IssmDouble x2, IssmDouble y2,int segment_id)=0;
		virtual IssmDouble MinEdgeLength(IssmDouble* xyz_list)=0;
		virtual IssmDouble Misfit(int modelenum,int observationenum,int weightsenum)=0;
		virtual IssmDouble MisfitArea(int weightsenum)=0;
		virtual Gauss*     NewGauss(void)=0;
		virtual Gauss*     NewGauss(int order)=0;
      virtual Gauss*     NewGauss(IssmDouble* xyz_list, IssmDouble* xyz_list_front,int order)=0;
      virtual Gauss*     NewGauss(IssmDouble* xyz_list, IssmDouble* xyz_list_front,int order_horiz,int order_vert)=0;
      virtual Gauss*     NewGauss(int point1,IssmDouble fraction1,IssmDouble fraction2,bool mainlyfloating,int order)=0;
		virtual Gauss*     NewGaussBase(int order)=0;
		virtual Gauss*     NewGaussLine(int vertex1,int vertex2,int order)=0;
		virtual Gauss*     NewGaussTop(int order)=0;
		virtual void       NodalFunctions(IssmDouble* basis,Gauss* gauss)=0;
		virtual void       NodalFunctionsDerivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss)=0;
		virtual void       NodalFunctionsDerivativesVelocity(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss)=0;
		virtual void       NodalFunctionsMINIDerivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss)=0;
		virtual void       NodalFunctionsPressure(IssmDouble* basis, Gauss* gauss)=0;
		virtual void       NodalFunctionsP1(IssmDouble* basis,Gauss* gauss)=0;
		virtual void       NodalFunctionsP1Derivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss)=0;
		virtual void       NodalFunctionsP2(IssmDouble* basis,Gauss* gauss)=0;
		virtual void       NodalFunctionsVelocity(IssmDouble* basis, Gauss* gauss)=0;
		virtual void       NodalFunctionsTensor(IssmDouble* basis, Gauss* gauss)=0;
		virtual int        NodalValue(IssmDouble* pvalue, int index, int natureofdataenum)=0;
		virtual void       NormalBase(IssmDouble* normal,IssmDouble* xyz_list)=0;
		virtual void       NormalSection(IssmDouble* normal,IssmDouble* xyz_list)=0;
		virtual void       NormalTop(IssmDouble* normal,IssmDouble* xyz_list)=0;
		virtual int        NumberofNodesPressure(void)=0;
		virtual int        NumberofNodesVelocity(void)=0;
		virtual void       PotentialUngrounding(Vector<IssmDouble>* potential_sheet_ungrounding)=0;
		virtual int        PressureInterpolation()=0;
		virtual void       ReduceMatrices(ElementMatrix* Ke,ElementVector* pe)=0;
		virtual void       ResetFSBasalBoundaryCondition()=0;
		virtual void       ResetHooks()=0;
		virtual void       ResetLevelsetFromSegmentlist(IssmDouble* segments,int numsegments){_error_("not implemented yet");};
		virtual void       SetControlInputsFromVector(IssmDouble* vector,int control_enum,int control_index)=0;
		virtual void       SetCurrentConfiguration(Elements* elements,Loads* loads,Nodes* nodes,Materials* materials,Parameters* parameters)=0;
		virtual void       SetTemporaryElementType(int element_type_in)=0;
	   virtual Element*   SpawnBasalElement(void)=0;
		virtual Element*   SpawnTopElement(void)=0;
		virtual IssmDouble StabilizationParameter(IssmDouble u, IssmDouble v, IssmDouble w, IssmDouble diameter, IssmDouble kappa)=0;
		virtual void	    StrainRateparallel(void)=0;
		virtual void	    StrainRateperpendicular(void)=0;
		virtual void	    StressIntensityFactor(void)=0;
		virtual IssmDouble SurfaceArea(void)=0;
		virtual int        TensorInterpolation()=0;
		virtual IssmDouble TimeAdapt()=0;
		virtual IssmDouble TotalFloatingBmb(void)=0;
		virtual IssmDouble TotalGroundedBmb(void)=0;
		virtual IssmDouble TotalSmb(void)=0;
		virtual void       Update(int index, IoModel* iomodel,int analysis_counter,int analysis_type,int finite_element)=0;
		virtual void       UpdateConstraintsExtrudeFromBase(void)=0;
		virtual void       UpdateConstraintsExtrudeFromTop(void)=0;
		virtual int        UpdatePotentialUngrounding(IssmDouble* potential_sheet_ungrounding,Vector<IssmDouble>* vec_nodes_on_iceshelf,IssmDouble* nodes_on_iceshelf)=0;
		virtual void       ValueP1DerivativesOnGauss(IssmDouble* dvalue,IssmDouble* values,IssmDouble* xyz_list,Gauss* gauss)=0;
		virtual void       ValueP1OnGauss(IssmDouble* pvalue,IssmDouble* values,Gauss* gauss)=0;
		virtual int        VelocityInterpolation()=0;
		virtual int        VertexConnectivity(int vertexindex)=0;
		virtual void       VerticalSegmentIndices(int** pindices,int* pnumseg)=0;
		virtual void       VerticalSegmentIndicesBase(int** pindices,int* pnumseg)=0;
		virtual void       ViscousHeating(IssmDouble* pphi,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vz_input)=0;
		virtual void       WriteLevelsetSegment(DataSet* segments){_error_("not implemented yet");};
		virtual void       ZeroLevelsetCoordinates(IssmDouble** pxyz_zero,IssmDouble* xyz_list,int levelsetenum)=0;

		#ifdef _HAVE_GIAIVINS_
		virtual void       GiaDeflection(Vector<IssmDouble>* wg,Vector<IssmDouble>* dwgdt,IssmDouble* x,IssmDouble* y)=0;
		#endif
		#ifdef _HAVE_ESA_
		virtual void          EsaGeodetic2D(Vector<IssmDouble>* pUp,Vector<IssmDouble>* pNorth,Vector<IssmDouble>* pEast,IssmDouble* xx,IssmDouble* yy)=0;
		virtual void          EsaGeodetic3D(Vector<IssmDouble>* pUp,Vector<IssmDouble>* pNorth,Vector<IssmDouble>* pEast,IssmDouble* latitude,IssmDouble* longitude,IssmDouble* radius,IssmDouble* xx,IssmDouble* yy,IssmDouble* zz,IssmDouble eartharea)=0;
		#endif
		#ifdef _HAVE_SEALEVELRISE_
		virtual IssmDouble    GetArea3D(void)=0;
		virtual IssmDouble    GetAreaSpherical(void)=0;
		virtual IssmDouble    OceanAverage(IssmDouble* Sg)=0;
		virtual IssmDouble    OceanArea(void)=0;
		virtual void          SealevelriseMomentOfInertia(IssmDouble* dI_list,IssmDouble* Sg_old,IssmDouble eartharea)=0; 
		virtual void          SealevelriseEustatic(Vector<IssmDouble>* pSgi,IssmDouble* peustatic,IssmDouble* latitude,IssmDouble* longitude,IssmDouble* radius,IssmDouble oceanarea,IssmDouble eartharea)=0;
		virtual void          SealevelriseNonEustatic(Vector<IssmDouble>* pSgo,IssmDouble* Sg_old,IssmDouble* latitude,IssmDouble* longitude,IssmDouble* radius,IssmDouble eartharea)=0;
		virtual void          SealevelriseGeodetic(Vector<IssmDouble>* pUp,Vector<IssmDouble>* pNorth,Vector<IssmDouble>* pEast,IssmDouble* Sg,IssmDouble* latitude,IssmDouble* longitude,IssmDouble* radius,IssmDouble* xx,IssmDouble* yy,IssmDouble* zz,IssmDouble eartharea)=0;
		#endif

};
#endif
