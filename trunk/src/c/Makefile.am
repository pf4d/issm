AM_CPPFLAGS = @NEOPZINCL@ @DAKOTAINCL@ @SHAPELIBINCL@ @PETSCINCL@ @SLEPCINCL@ @AMPIINCL@ @MPIINCL@ @METISINCL@ @CHACOINCL@ @SCOTCHINCL@ @PLAPACKINCL@ @BLASLAPACKINCL@ @MKLINCL@ @MUMPSINCL@ @TRIANGLEINCL@ @SPAIINCL@ @HYPREINCL@ @PROMETHEUSINCL@ @SUPERLUINCL@ @SPOOLESINCL@ @PASTIXINCL@ @MLINCL@ @TAOINCL@ @ADIC2INCL@ @ADOLCINCL@ @GSLINCL@ @BOOSTINCL@ @ANDROID_NDKINCL@ @METEOIOINCL@ @SNOWPACKINCL@ @PROJ4INCL@

AUTOMAKE_OPTIONS = subdir-objects

EXEEXT=$(ISSMEXT)

#Library declaration {{{
if !WINDOWS
lib_LTLIBRARIES = libISSMCore.la libISSMOverload.la 
if WRAPPERS
lib_LTLIBRARIES += libISSMModules.la
endif
else
noinst_LTLIBRARIES = libISSMCore.la libISSMOverload.la 
if WRAPPERS
noinst_LTLIBRARIES += libISSMModules.la
endif
endif
#}}}

#Core sources
#BAMG sources  {{{
issm_sources = 
if BAMG
issm_sources += ./bamg/BamgGeom.cpp\
					 ./bamg/BamgMesh.cpp\
					 ./bamg/BamgOpts.cpp\
					 ./bamg/CrackedEdge.cpp\
					 ./bamg/Curve.cpp\
					 ./bamg/Edge.cpp\
					 ./bamg/GeomEdge.cpp\
					 ./bamg/GeomSubDomain.cpp\
					 ./bamg/GeomVertex.cpp\
					 ./bamg/Geometry.cpp\
					 ./bamg/ListofIntersectionTriangles.cpp\
					 ./bamg/EigenMetric.cpp\
					 ./bamg/Metric.cpp\
					 ./bamg/BamgQuadtree.cpp\
					 ./bamg/SetOfE4.cpp\
					 ./bamg/SubDomain.cpp\
					 ./bamg/AdjacentTriangle.cpp\
					 ./bamg/Triangle.cpp\
					 ./bamg/BamgVertex.cpp\
					 ./bamg/VertexOnEdge.cpp\
					 ./bamg/VertexOnGeom.cpp\
					 ./bamg/VertexOnVertex.cpp\
					 ./bamg/Mesh.cpp\
					 ./shared/Bamg/BigPrimeNumber.cpp\
					 ./modules/Bamgx/Bamgx.cpp\
					 ./modules/BamgConvertMeshx/BamgConvertMeshx.cpp\
					 ./modules/BamgTriangulatex/BamgTriangulatex.cpp
endif
#}}}
#Core sources{{{
issm_sources += ./datastructures/DataSet.cpp\
					./classes/gauss/GaussSeg.cpp\
					./classes/gauss/GaussTria.cpp\
					./classes/gauss/GaussTetra.cpp\
					./classes/gauss/GaussPenta.cpp\
					./classes/IoModel.cpp\
					./classes/FemModel.cpp\
					./classes/Loads/Friction.cpp\
					./classes/Inputs/TransientInput.cpp\
					./classes/Constraints/SpcTransient.cpp\
					./classes/DependentObject.cpp\
					./classes/DofIndexing.cpp\
					./classes/Contours.cpp\
					./classes/Nodes.cpp\
					./classes/Vertices.cpp\
					./classes/Node.cpp\
					./classes/Vertex.cpp\
					./classes/Hook.cpp\
					./classes/ExternalResults/Results.cpp\
					./classes/Elements/Element.cpp\
					./classes/Elements/Elements.cpp\
					./classes/Elements/ElementHook.cpp\
					./classes/Elements/Seg.cpp\
					./classes/Elements/SegRef.cpp\
					./classes/Elements/Tria.cpp\
					./classes/Elements/TriaRef.cpp\
					./classes/Elements/Tetra.cpp\
					./classes/Elements/TetraRef.cpp\
					./classes/Elements/Penta.cpp\
					./classes/Elements/PentaRef.cpp\
					./classes/Inputs/Inputs.cpp\
					./classes/Inputs/SegInput.cpp\
					./classes/Inputs/TriaInput.cpp\
					./classes/Inputs/BoolInput.cpp\
					./classes/Inputs/IntInput.cpp\
					./classes/Inputs/DoubleInput.cpp\
					./classes/Inputs/DoubleArrayInput.cpp\
					./classes/Inputs/DatasetInput.cpp\
					./classes/Materials/Materials.cpp\
					./classes/Materials/Matice.cpp\
					./classes/Materials/Matestar.cpp\
					./classes/Materials/Matpar.cpp\
					./classes/Constraints/Constraints.cpp\
					./classes/Constraints/SpcStatic.cpp\
					./classes/Constraints/SpcDynamic.cpp\
					./classes/Loads/Loads.cpp\
					./classes/Loads/Penpair.cpp\
					./classes/Loads/Pengrid.cpp\
					./classes/Loads/Moulin.cpp\
					./classes/Loads/Numericalflux.cpp\
					./classes/Loads/Neumannflux.cpp\
					./classes/matrix/ElementMatrix.cpp\
					./classes/matrix/ElementVector.cpp\
					./classes/Params/Parameters.cpp\
					./classes/Params/BoolParam.cpp\
					./classes/Params/IntParam.cpp\
					./classes/Params/IntVecParam.cpp\
					./classes/Params/IntMatParam.cpp\
					./classes/Params/DoubleParam.cpp\
					./classes/Params/FileParam.cpp\
					./classes/Params/StringArrayParam.cpp\
					./classes/Params/DoubleMatParam.cpp\
					./classes/Params/DoubleTransientMatParam.cpp\
					./classes/Params/DoubleMatArrayParam.cpp\
					./classes/Params/DoubleVecParam.cpp\
					./classes/Params/StringParam.cpp\
					./classes/Params/MatrixParam.cpp\
					./classes/Params/VectorParam.cpp\
					./classes/Params/TransientParam.cpp\
					./classes/Params/DataSetParam.cpp\
					./classes/Profiler.cpp\
					./shared/MemOps/MemOps.cpp\
					./shared/Matrix/MatrixUtils.cpp\
					./shared/io/Disk/pfopen.cpp\
					./shared/io/Disk/pfclose.cpp\
					./shared/io/Disk/WriteLockFile.cpp\
					./shared/io/Print/PrintfFunction.cpp\
					./shared/io/Comm/IssmComm.cpp\
					./shared/io/Marshalling/IoCodeConversions.cpp \
					./shared/LatLong/Ll2xyx.cpp\
					./shared/LatLong/Xy2llx.cpp\
					./shared/FSanalyticals/fsanalyticals.cpp\
					./shared/Enum/EnumToStringx.cpp\
					./shared/Enum/StringToEnumx.cpp\
					./shared/Numerics/Verbosity.cpp\
					./shared/Numerics/GaussPoints.cpp\
					./shared/Numerics/cross.cpp\
					./shared/Numerics/isnan.cpp\
					./shared/Numerics/cubic.cpp\
					./shared/Numerics/NewtonSolveDnorm.cpp\
					./shared/Numerics/extrema.cpp\
					./shared/Numerics/legendre.cpp\
					./shared/Numerics/XZvectorsToCoordinateSystem.cpp\
					./shared/Exceptions/Exceptions.cpp\
					./shared/Sorting/binary_search.cpp\
					./shared/Elements/Cuffey.cpp\
					./shared/Elements/BuddJacka.cpp\
					./shared/Elements/CuffeyTemperate.cpp\
					./shared/Elements/StressIntensityIntegralWeight.cpp\
					./shared/Elements/Paterson.cpp\
					./shared/Elements/Arrhenius.cpp\
					./shared/Elements/LliboutryDuval.cpp \
					./shared/Elements/PrintArrays.cpp\
					./shared/Elements/PddSurfaceMassBalance.cpp\
					./shared/Elements/ComputeDelta18oTemperaturePrecipitation.cpp\
					./shared/Elements/ComputeMungsmTemperaturePrecipitation.cpp\
					./shared/Elements/ComputeD18OTemperaturePrecipitationFromPD.cpp\
					./shared/Elements/DrainageFunctionWaterfraction.cpp\
					./shared/Elements/EstarComponents.cpp\
					./shared/String/DescriptorIndex.cpp\
					./toolkits/issm/IssmToolkitUtils.cpp\
					./toolkits/issm/IssmSolver.cpp\
					./toolkits/mpi/issmmpi.cpp\
					./toolkits/mpi/commops/DetermineLocalSize.cpp\
					./toolkits/mpi/commops/DetermineGlobalSize.cpp\
					./toolkits/mpi/commops/DetermineRowRankFromLocalSize.cpp\
					./toolkits/mpi/commops/GetOwnershipBoundariesFromRange.cpp\
					./toolkits/ToolkitOptions.cpp\
					./modules/ModelProcessorx/ModelProcessorx.cpp\
					./modules/ModelProcessorx/ElementsAndVerticesPartitioning.cpp\
					./modules/ModelProcessorx/NodesPartitioning.cpp\
					./modules/ModelProcessorx/EdgesPartitioning.cpp\
					./modules/ModelProcessorx/FacesPartitioning.cpp\
					./modules/ModelProcessorx/CreateParameters.cpp\
					./modules/ModelProcessorx/Autodiff/CreateParametersAutodiff.cpp\
					./modules/ModelProcessorx/CreateFaces.cpp\
					./modules/ModelProcessorx/CreateEdges.cpp\
					./modules/ModelProcessorx/CreateSingleNodeToElementConnectivity.cpp\
					./modules/ModelProcessorx/CreateNumberNodeToElementConnectivity.cpp\
					./modules/ModelProcessorx/CreateElementsVerticesAndMaterials.cpp\
					./modules/ModelProcessorx/CreateNodes.cpp\
					./modules/ParseToolkitsOptionsx/ParseToolkitsOptionsx.cpp\
					./modules/NodesDofx/NodesDofx.cpp\
					./modules/NodalValuex/NodalValuex.cpp\
					./modules/VerticesDofx/VerticesDofx.cpp\
					./modules/VertexCoordinatesx/VertexCoordinatesx.cpp\
					./modules/OutputResultsx/OutputResultsx.cpp\
					./modules/InputDepthAverageAtBasex/InputDepthAverageAtBasex.cpp\
					./modules/InputDuplicatex/InputDuplicatex.cpp\
					./modules/InputExtrudex/InputExtrudex.cpp\
					./modules/SurfaceAreax/SurfaceAreax.cpp\
					./modules/AllocateSystemMatricesx/AllocateSystemMatricesx.cpp\
					./modules/CreateJacobianMatrixx/CreateJacobianMatrixx.cpp\
					./modules/SystemMatricesx/SystemMatricesx.cpp\
					./modules/CreateNodalConstraintsx/CreateNodalConstraintsx.cpp\
					./modules/UpdateDynamicConstraintsx/UpdateDynamicConstraintsx.cpp\
					./modules/IoModelToConstraintsx/IoModelToConstraintsx.cpp\
					./modules/SetActiveNodesLSMx/SetActiveNodesLSMx.cpp\
					./modules/InputUpdateFromConstantx/InputUpdateFromConstantx.cpp\
					./modules/InputUpdateFromSolutionx/InputUpdateFromSolutionx.cpp\
					./modules/GeothermalFluxx/GeothermalFluxx.cpp\
					./modules/GetSolutionFromInputsx/GetSolutionFromInputsx.cpp\
					./modules/GetVectorFromInputsx/GetVectorFromInputsx.cpp\
					./modules/InputUpdateFromVectorx/InputUpdateFromVectorx.cpp\
					./modules/FloatingiceMeltingRatex/FloatingiceMeltingRatex.cpp\
					./modules/ConfigureObjectsx/ConfigureObjectsx.cpp\
					./modules/SpcNodesx/SpcNodesx.cpp\
					./modules/SurfaceMassBalancex/SurfaceMassBalancex.cpp\
					./modules/SurfaceMassBalancex/Gembx.cpp\
					./modules/Reducevectorgtofx/Reducevectorgtofx.cpp\
					./modules/Reduceloadx/Reduceloadx.cpp\
					./modules/ConstraintsStatex/ConstraintsStatex.cpp\
					./modules/ResetConstraintsx/ResetConstraintsx.cpp\
					./modules/ResetFSBasalBoundaryConditionx/ResetFSBasalBoundaryConditionx.cpp\
					./modules/Solverx/Solverx.cpp\
					./modules/VecMergex/VecMergex.cpp\
					./modules/Mergesolutionfromftogx/Mergesolutionfromftogx.cpp\
					./cores/ProcessArguments.cpp\
					./cores/ResetBoundaryConditions.cpp\
					./cores/WrapperCorePointerFromSolutionEnum.cpp\
					./cores/CorePointerFromSolutionEnum.cpp\
					./cores/ad_core.cpp\
					./cores/adgradient_core.cpp\
					./main/EnvironmentInit.cpp\
					./main/EnvironmentFinalize.cpp\
					./analyses/EnumToAnalysis.cpp\
					./solutionsequences/solutionsequence_la.cpp\
					./solutionsequences/solutionsequence_la_theta.cpp\
					./solutionsequences/solutionsequence_linear.cpp\
					./solutionsequences/solutionsequence_nonlinear.cpp\
					./solutionsequences/solutionsequence_newton.cpp\
					./solutionsequences/solutionsequence_fct.cpp\
					./solutionsequences/convergence.cpp\
					./classes/Options/Options.cpp\
					./classes/Options/OptionUtilities.cpp\
					./classes/RiftStruct.cpp\
					./modules/ModelProcessorx/Transient/UpdateElementsTransient.cpp \
					./cores/transient_core.cpp\
					./cores/steadystate_core.cpp\
					./cores/masstransport_core.cpp\
					./cores/depthaverage_core.cpp\
					./cores/extrudefrombase_core.cpp\
					./cores/extrudefromtop_core.cpp\
					./cores/thermal_core.cpp\
					./cores/smb_core.cpp\
					./solutionsequences/solutionsequence_thermal_nonlinear.cpp\
					./modules/ControlInputSetGradientx/ControlInputSetGradientx.cpp\
					./modules/GetVectorFromControlInputsx/GetVectorFromControlInputsx.cpp\
					./modules/SetControlInputsFromVectorx/SetControlInputsFromVectorx.cpp\
					./modules/ModelProcessorx/Control/CreateParametersControl.cpp\
					./modules/ModelProcessorx/Control/UpdateElementsAndMaterialsControl.cpp\
					./modules/SurfaceAbsVelMisfitx/SurfaceAbsVelMisfitx.cpp\
					./modules/SurfaceRelVelMisfitx/SurfaceRelVelMisfitx.cpp\
					./modules/SurfaceLogVelMisfitx/SurfaceLogVelMisfitx.cpp\
					./modules/SurfaceLogVxVyMisfitx/SurfaceLogVxVyMisfitx.cpp\
					./modules/SurfaceAverageVelMisfitx/SurfaceAverageVelMisfitx.cpp\
					./modules/ThicknessAbsMisfitx/ThicknessAbsMisfitx.cpp\
					./modules/Gradjx/Gradjx.cpp\
					./modules/DragCoefficientAbsGradientx/DragCoefficientAbsGradientx.cpp\
					./modules/ThicknessAlongGradientx/ThicknessAlongGradientx.cpp\
					./modules/ThicknessAcrossGradientx/ThicknessAcrossGradientx.cpp\
					./modules/RheologyBbarAbsGradientx/RheologyBbarAbsGradientx.cpp\
					./modules/RheologyBAbsGradientx/RheologyBAbsGradientx.cpp\
					./classes/Inputs/ControlInput.cpp\
					./shared/Numerics/BrentSearch.cpp\
					./cores/control_core.cpp\
					./cores/controltao_core.cpp\
					./cores/controlad_core.cpp\
					./cores/controlm1qn3_core.cpp\
					./cores/controlvalidation_core.cpp\
					./cores/adjointstressbalance_core.cpp\
					./cores/adjointbalancethickness_core.cpp\
					./cores/adjointbalancethickness2_core.cpp\
					./cores/AdjointCorePointerFromSolutionEnum.cpp\
					./solutionsequences/solutionsequence_adjoint_linear.cpp\
					./cores/hydrology_core.cpp\
					./solutionsequences/solutionsequence_hydro_nonlinear.cpp\
					./cores/stressbalance_core.cpp\
					./solutionsequences/solutionsequence_stokescoupling_nonlinear.cpp\
					./cores/balancethickness_core.cpp \
					./cores/balancethickness2_core.cpp \
					./cores/balancevelocity_core.cpp \
					./cores/dummy_core.cpp\
					./cores/surfaceslope_core.cpp\
					./cores/bedslope_core.cpp\
					./cores/meshdeformation_core.cpp\
					./cores/damage_core.cpp\
					./cores/levelsetfunctionslope_core.cpp\
					./cores/movingfront_core.cpp\
					./modules/GroundinglineMigrationx/GroundinglineMigrationx.cpp\
					./classes/Loads/Riftfront.cpp\
					./modules/ConstraintsStatex/RiftConstraintsState.cpp\
					./modules/ModelProcessorx/CreateOutputDefinitions.cpp\
					./modules/OutputDefinitionsResponsex/OutputDefinitionsResponsex.cpp\
					./modules/InterpFromMeshToMesh2dx/InterpFromMeshToMesh2dx.cpp\
					./classes/Inputs/PentaInput.cpp\
					./classes/Inputs/TetraInput.cpp
#}}}
#DAKOTA sources  {{{
if DAKOTA
issm_sources +=  ./classes/Dakota/IssmDirectApplicInterface.h\
					  ./classes/Dakota/IssmParallelDirectApplicInterface.cpp\
					  ./modules/InputUpdateFromDakotax/InputUpdateFromDakotax.cpp\
					  ./modules/InputUpdateFromVectorDakotax/InputUpdateFromVectorDakotax.cpp\
					  ./modules/InputUpdateFromMatrixDakotax/InputUpdateFromMatrixDakotax.cpp\
					  ./modules/AverageOntoPartitionx/AverageOntoPartitionx.cpp\
					  ./modules/ModelProcessorx/Dakota/CreateParametersDakota.cpp\
					  ./modules/ModelProcessorx/Dakota/UpdateElementsAndMaterialsDakota.cpp\
					  ./cores/dakota_core.cpp
endif
#}}}
#Petsc sources  {{{
if PETSC
issm_sources += ./toolkits/petsc\
					./toolkits/petsc/patches\
					./toolkits/petsc/patches/VecToMPISerial.cpp\
					./toolkits/petsc/patches/MatToSerial.cpp\
					./toolkits/petsc/patches/VecMerge.cpp\
					./toolkits/petsc/patches/NewVec.cpp\
					./toolkits/petsc/patches/PetscOptionsDetermineSolverType.cpp\
					./toolkits/petsc/patches/PetscOptionsInsertMultipleString.cpp\
					./toolkits/petsc/patches/NewMat.cpp\
					./toolkits/petsc/patches/VecFree.cpp\
					./toolkits/petsc/patches/KSPFree.cpp\
					./toolkits/petsc/patches/MatFree.cpp\
					./toolkits/petsc/patches/MatMultPatch.cpp\
					./toolkits/petsc/patches/ISSMToPetscMatrixType.cpp\
					./toolkits/petsc/patches/ISSMToPetscInsertMode.cpp\
					./toolkits/petsc/patches/ISSMToPetscNormMode.cpp\
					./toolkits/petsc/objects/PetscMat.cpp\
					./toolkits/petsc/objects/PetscVec.cpp\
					./toolkits/petsc/objects/PetscSolver.cpp
endif
#}}}
#Mumps sources  {{{
if MUMPS
issm_sources += ./toolkits/mumps/MumpsSolve.cpp
endif
#}}}
#Gsl sources  {{{
if GSL
issm_sources += ./toolkits/gsl/DenseGslSolve.cpp
endif
#}}}
#proj.4 sources  {{{
if PROJ4
issm_sources += ./modules/CoordinateSystemTransformx/CoordinateSystemTransformx.cpp
endif
#}}}
#Analyses{{{
if ADJOINTBALANCETHICKNESS
issm_sources += ./analyses/AdjointBalancethicknessAnalysis.cpp
endif
if ADJOINTBALANCETHICKNESS2
issm_sources += ./analyses/AdjointBalancethickness2Analysis.cpp
endif
if ADJOINTHORIZ
issm_sources += ./analyses/AdjointHorizAnalysis.cpp
endif
if BALANCETHICKNESS
issm_sources += ./analyses/BalancethicknessAnalysis.cpp
endif
if BALANCETHICKNESS2
issm_sources += ./analyses/Balancethickness2Analysis.cpp
endif
if BALANCETHICKNESSSOFT
issm_sources += ./analyses/BalancethicknessSoftAnalysis.cpp
endif
if BALANCEVELOCITY
issm_sources += ./analyses/BalancevelocityAnalysis.cpp
endif
if L2PROJECTIONBASE
issm_sources += ./analyses/L2ProjectionBaseAnalysis.cpp
endif
if DAMAGEEVOLUTION
issm_sources += ./analyses/DamageEvolutionAnalysis.cpp
issm_sources += ./modules/Damagex/Damagex.cpp
endif
if STRESSBALANCE
issm_sources += ./analyses/StressbalanceAnalysis.cpp
endif
if UZAWAPRESSURE
issm_sources += ./analyses/UzawaPressureAnalysis.cpp
endif
if STRESSBALANCESIA
issm_sources += ./analyses/StressbalanceSIAAnalysis.cpp
endif
if STRESSBALANCEVERTICAL
issm_sources += ./analyses/StressbalanceVerticalAnalysis.cpp
endif
if ENTHALPY
issm_sources += ./analyses/EnthalpyAnalysis.cpp
endif
if HYDROLOGYSHREVE
issm_sources += ./analyses/HydrologyShreveAnalysis.cpp
endif
if HYDROLOGYSOMMERS
issm_sources += ./analyses/HydrologySommersAnalysis.cpp
endif
if HYDROLOGYDCINEFFICIENT
issm_sources += ./analyses/HydrologyDCInefficientAnalysis.cpp
endif
if HYDROLOGYDCEFFICIENT
issm_sources += ./analyses/HydrologyDCEfficientAnalysis.cpp
endif
if L2PROJECTIONEPL
issm_sources += ./analyses/L2ProjectionEPLAnalysis.cpp
endif
if MELTING
issm_sources += ./analyses/MeltingAnalysis.cpp
endif
if MASSTRANSPORT
issm_sources += ./analyses/MasstransportAnalysis.cpp
endif
if SMB
issm_sources += ./analyses/SmbAnalysis.cpp
endif

if FREESURFACEBASE
issm_sources += ./analyses/FreeSurfaceBaseAnalysis.cpp
endif
if FREESURFACETOP
issm_sources += ./analyses/FreeSurfaceTopAnalysis.cpp
endif
if EXTRUDEFROMBASE
issm_sources += ./analyses/ExtrudeFromBaseAnalysis.cpp
endif
if EXTRUDEFROMTOP
issm_sources += ./analyses/ExtrudeFromTopAnalysis.cpp
endif
if DEPTHAVERAGE
issm_sources += ./analyses/DepthAverageAnalysis.cpp
endif
if THERMAL
issm_sources += ./analyses/ThermalAnalysis.cpp
endif
if SMOOTH
issm_sources += ./analyses/SmoothAnalysis.cpp
endif
if MESHDEFORMATION
issm_sources += ./analyses/MeshdeformationAnalysis.cpp
endif
if LEVELSET
issm_sources += ./analyses/LevelsetAnalysis.cpp
issm_sources += ./modules/Calvingx/Calvingx.cpp
endif
if EXTRAPOLATION
issm_sources += ./analyses/ExtrapolationAnalysis.cpp
endif
#}}}
#Gia sources  {{{
if GIAIVINS
issm_sources +=  ./cores/gia_core.cpp\
					./analyses/GiaIvinsAnalysis.cpp\
					./modules/GiaDeflectionCorex/GiaDeflectionCorex.cpp\
					./modules/GiaDeflectionCorex/distme.f\
					./modules/GiaDeflectionCorex/freed.f\
					./modules/GiaDeflectionCorex/ojrule.f\
					./modules/GiaDeflectionCorex/pwise.f\
					./modules/GiaDeflectionCorex/qwise.f\
					./modules/GiaDeflectionCorex/stot.f\
					./modules/GiaDeflectionCorex/what0.f
endif
#}}}
#Esa sources  {{{
if ESA 
issm_sources +=  ./cores/esa_core.cpp\
					./analyses/EsaAnalysis.cpp
endif
#}}}
#Slr sources  {{{
if SEALEVELRISE
issm_sources +=  ./cores/sealevelrise_core.cpp\
				 ./cores/sealevelrise_core_eustatic.cpp\
				 ./cores/sealevelrise_core_noneustatic.cpp\
				 ./analyses/SealevelriseAnalysis.cpp
endif
#}}}
#Metis sources  {{{
if METIS
issm_sources += ./toolkits/metis/patches/METIS_PartMeshNodalPatch.cpp
endif
#}}}
#Esmf sources  {{{
if ESMF
issm_sources += ./main/esmfbinders.cpp
endif
#}}}

#Wrapper sources
#Bamg sources  {{{

#}}}
#Kml sources  {{{
kml_sources = ./modules/Exp2Kmlx/Exp2Kmlx.cpp\
				  ./modules/Kml2Expx/Kml2Expx.cpp\
				  ./modules/Shp2Kmlx/Shp2Kmlx.cpp\
				  ./modules/KMLFileReadx/KMLFileReadx.cpp\
				  ./modules/KMLMeshWritex/KMLMeshWritex.cpp\
				  ./modules/KMLOverlayx/KMLOverlayx.cpp\
				  ./kml/KML_Attribute.cpp\
				  ./kml/KML_Comment.cpp\
				  ./kml/KML_ColorStyle.cpp\
				  ./kml/KML_Container.cpp\
				  ./kml/KML_Document.cpp\
				  ./kml/KML_Feature.cpp\
				  ./kml/KML_File.cpp\
				  ./kml/KML_Folder.cpp\
				  ./kml/KML_Geometry.cpp\
				  ./kml/KML_GroundOverlay.cpp\
				  ./kml/KML_Icon.cpp\
				  ./kml/KML_LatLonBox.cpp\
				  ./kml/KML_LinearRing.cpp\
				  ./kml/KML_LineString.cpp\
				  ./kml/KML_LineStyle.cpp\
				  ./kml/KML_MultiGeometry.cpp\
				  ./kml/KML_Object.cpp\
				  ./kml/KML_Overlay.cpp\
				  ./kml/KML_Point.cpp\
				  ./kml/KML_Placemark.cpp\
				  ./kml/KML_Polygon.cpp\
				  ./kml/KML_PolyStyle.cpp\
				  ./kml/KML_Style.cpp\
				  ./kml/KML_StyleSelector.cpp\
				  ./kml/KML_SubStyle.cpp\
				  ./kml/KML_Unknown.cpp\
				  ./kml/KMLFileReadUtils.cpp
#}}}
#NEOPZ sources  {{{
neopz_sources = ./classes/AdaptiveMeshRefinement.cpp
#}}}
#Modules sources{{{
modules_sources= ./shared/Threads/LaunchThread.cpp\
			./shared/Threads/PartitionRange.cpp\
			./shared/Exp/exp.cpp\
			./shared/TriMesh/AssociateSegmentToElement.cpp\
			./shared/TriMesh/GridInsideHole.cpp\
			./shared/TriMesh/OrderSegments.cpp\
			./shared/TriMesh/SplitMeshForRifts.cpp\
			./shared/TriMesh/TriMeshUtils.cpp\
			./modules/TriMeshx/TriMeshx.cpp\
			./modules/TriMeshProcessRiftsx/TriMeshProcessRiftsx.cpp\
			./modules/PointCloudFindNeighborsx/PointCloudFindNeighborsx.cpp\
			./modules/PointCloudFindNeighborsx/PointCloudFindNeighborsxt.cpp\
			./modules/InterpFromGridToMeshx/InterpFromGridToMeshx.cpp\
			./modules/InterpFromMesh2dx/InterpFromMesh2dx.cpp\
			./modules/InterpFromMesh2dx/InterpFromMesh2dxt.cpp\
			./modules/InterpFromMeshToMesh3dx/InterpFromMeshToMesh3dx.cpp\
			./modules/InterpFromMeshToGridx/InterpFromMeshToGridx.cpp\
			./modules/MeshProfileIntersectionx/MeshProfileIntersectionx.cpp\
			./modules/ContourToMeshx/ContourToMeshx.cpp\
			./modules/ContourToMeshx/ContourToMeshxt.cpp\
			./modules/ExpToLevelSetx/ExpToLevelSetx.cpp\
			./modules/ExpToLevelSetx/ExpToLevelSetxt.cpp\
			./modules/ContourToNodesx/ContourToNodesx.cpp\
			./modules/DistanceToMaskBoundaryx/DistanceToMaskBoundaryx.cpp\
			./modules/DistanceToMaskBoundaryx/DistanceToMaskBoundaryxt.cpp\
			./modules/NodeConnectivityx/NodeConnectivityx.cpp\
			./modules/ElementConnectivityx/ElementConnectivityx.cpp\
			./modules/PropagateFlagsFromConnectivityx/PropagateFlagsFromConnectivityx.cpp
if CHACO
modules_sources+= ./modules/Chacox/Chacox.cpp\
						./modules/Chacox/input_parse.cpp\
						./modules/Chacox/chaco_seconds.cpp\
						./modules/Chacox/user_params.cpp
endif
if SCOTCH
modules_sources+= ./modules/Scotchx/Scotchx.cpp
endif
#}}}

#kriging (WRAPPER and executable)
#Kriging sources  {{{
if KRIGING
issm_sources += ./classes/kriging/Observations.cpp\
					./classes/kriging/GaussianVariogram.cpp\
					./classes/kriging/ExponentialVariogram.cpp\
					./classes/kriging/SphericalVariogram.cpp\
					./classes/kriging/PowerVariogram.cpp\
					./classes/kriging/Quadtree.cpp\
					./classes/kriging/Covertree.cpp\
					./classes/kriging/Observation.cpp\
					./modules/Krigingx/pKrigingx.cpp

modules_sources +=./modules/Krigingx/Krigingx.cpp\
						./modules/Krigingx/pKrigingx.cpp
endif
#}}}
#Library flags and sources {{{
ALLCXXFLAGS= $(CXXFLAGS) $(CXXOPTFLAGS) 

libISSMCore_la_SOURCES  = $(issm_sources)
libISSMCore_la_CXXFLAGS = $(ALLCXXFLAGS) $(DAKOTAFLAGS)
libISSMCore_la_FFLAGS = $(AM_FFLAGS)

if !WINDOWS
if !STANDALONE_LIBRARIES
libISSMCore_la_LIBADD = $(PETSCLIB) $(TAOLIB) $(M1QN3LIB) $(PLAPACKLIB) $(MUMPSLIB) $(SUPERLULIB) $(SPOOLESLIB) $(SCALAPACKLIB) $(BLACSLIB) $(HYPRELIB) $(SPAILIB) $(PROMETHEUSLIB) $(PASTIXLIB) $(MLLIB) $(DAKOTALIB) $(METISLIB) $(CHACOLIB) $(SCOTCHLIB) $(BLASLAPACKLIB) $(MKLLIB) $(MPILIB) $(MATHLIB) $(GRAPHICSLIB) $(MULTITHREADINGLIB) $(OSLIBS) $(GSLLIB)   $(ADOLCLIB) $(AMPILIB) $(METEOIOLIB) $(SNOWPACKLIB)
if FORTRAN
libISSMCore_la_LIBADD += $(FLIBS) $(FORTRANLIB)
endif
endif
endif

if NEOPZ
libISSMCore_la_SOURCES += $(neopz_sources)
endif

if WRAPPERS
libISSMModules_la_SOURCES = $(modules_sources)
libISSMModules_la_SOURCES += $(bamg_sources)
if KRIGING
libISSMModules_la_SOURCES += $(kriging_sources)
endif
if KML
libISSMModules_la_SOURCES += $(kml_sources)
endif

libISSMModules_la_CXXFLAGS = $(ALLCXXFLAGS)
if !WINDOWS
if STANDALONE_LIBRARIES
libISSMModules_la_LIBADD = ./libISSMCore.la
else
libISSMModules_la_LIBADD = ./libISSMCore.la $(TRIANGLELIB) $(CHACOLIB)
endif
endif
endif

if VERSION
AM_LDFLAGS =
else
AM_LDFLAGS = -avoid-version
endif

if STANDALONE_EXECUTABLES
issm_LDFLAGS = -static
endif

if STANDALONE_LIBRARIES
libISSMCore_la_LDFLAGS = -static
libISSMOverload_la_LDFLAGS = -static
if WRAPPERS
libISSMModules_la_LDFLAGS = -static
else
libISSMModules_la_LDFLAGS = 
endif
else
libISSMCore_la_LDFLAGS = 
libISSMOverload_la_LDFLAGS =
endif

#}}}
#Overload library, to overload any non-standard symbols. {{{
libISSMOverload_la_SOURCES = ./shared/String/ApiPrintf.cpp

libISSMOverload_la_CFLAGS  = -D_C_ -fPIC $(COPTFLAGS) $(CFLAGS)
#}}}
#Executable {{{
if ANDROID
if ANDROIDEXE
bin_PROGRAMS = issm  issm_slr
else
bin_PROGRAMS = 
endif
else
if JAVASCRIPT
bin_PROGRAMS =
else
bin_PROGRAMS = issm  issm_slr
endif
endif

#Standard libraries
LDADD = ./libISSMCore.la ./libISSMOverload.la
if WRAPPERS
LDADD += ./libISSMModules.la $(TRIANGLELIB)
endif

#External packages
LDADD +=  $(NEOPZLIB) $(TAOLIB) $(M1QN3LIB) $(PLAPACKLIB) $(MUMPSLIB) $(SUPERLULIB) $(SPOOLESLIB) $(SCALAPACKLIB) $(BLACSLIB) $(PETSCLIB) $(HYPRELIB) $(SPAILIB) $(PROMETHEUSLIB) $(PASTIXLIB) $(MLLIB) $(DAKOTALIB) $(METISLIB) $(CHACOLIB) $(SCOTCHLIB) $(BLASLAPACKLIB) $(MKLLIB) $(MPILIB)  $(MATHLIB) $(GRAPHICSLIB) $(MULTITHREADINGLIB) $(OSLIBS) $(GSLLIB) $(AMPILIB) $(ADOLCLIB) $(MPILIB) $(METEOIOLIB) $(SNOWPACKLIB) $(PROJ4LIB)

if FORTRAN
LDADD += $(FLIBS) $(FORTRANLIB)
endif

issm_SOURCES = main/issm.cpp
issm_CXXFLAGS= $(CXXFLAGS) $(CXXOPTFLAGS)

issm_slr_SOURCES = main/issm_slr.cpp
issm_slr_CXXFLAGS= $(CXXFLAGS) $(CXXOPTFLAGS)

if OCEAN
bin_PROGRAMS += issm_ocean
issm_ocean_SOURCES = main/issm_ocean.cpp
issm_ocean_CXXFLAGS= $(CXXFLAGS) $(CXXOPTFLAGS)
endif

if KRIGING
bin_PROGRAMS += kriging
kriging_SOURCES = main/kriging.cpp
kriging_CXXFLAGS= $(CXXFLAGS) $(CXXOPTFLAGS) $(COPTFLAGS) 
endif

if ISSM_DAKOTA
bin_PROGRAMS += issm_dakota
issm_dakota_SOURCES = main/issm_dakota.cpp
issm_dakota_CXXFLAGS= $(CXXFLAGS) $(CXXOPTFLAGS) $(COPTFLAGS) 
issm_dakota_LDADD= $(LDADD)
endif

#}}}
#Automatic differentiation (must be done at the end) {{{
if ADIC2 
if !WINDOWS
lib_LTLIBRARIES += libAD.la libISSMRose.la
else
noinst_LTLIBRARIES += libAD.la libISSMRose.la
endif

#ADIC2 library, for automatic differentiation 
#libAD_a_SOURCES = ./mini1.ad.c
libAD_la_SOURCES = 
libAD_la_CFLAGS = -D_C_ -fPIC $(COPTFLAGS)

#test rose preprocessing 
%.r2cpp.cpp : %.cpp
	testTranslator -rose:o $@ -rose:skipfinalCompileStep -DHAVE_CONFIG_H -D_C_ -I. -I../.. $(INCLUDES) $<
libISSMRose_la_SOURCES = $(libISSMCore_a_SOURCES:.cpp=.r2cpp.cpp)
libISSMRose_la_CXXFLAGS= -D_C_ $(CXXOPTFLAGS)

#Automatic differentiation rules: 
%.ad.c: %.c
	adic2 -mforward  $< --nary

#Executable
bin_PROGRAMS +=  issmRose.exe
issmRose_exe_LDADD = ./libISSMRose.a $(LDADD)
issmRose_exe_SOURCES = main/issm.cpp
issmRose_exe_CXXFLAGS= $(CXXFLAGS) $(CXXOPTFLAGS) $(COPTFLAGS) 
LDADD +=  $(ADIC2LIB) 
endif #}}}
