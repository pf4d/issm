#To generate list of directories
#
#find . -type d -exec ls -d {} \;

bin_SCRIPTS = 
if WRAPPERS
if MATLAB
if !DEVELOPMENT
bin_SCRIPTS += ${ISSM_DIR}/src/m/qmu/*.m \
					${ISSM_DIR}/src/m/archive/*.m \
					${ISSM_DIR}/src/m/qmu/setupdesign/*.m \
					${ISSM_DIR}/src/m/qmu/plot/*.m \
					${ISSM_DIR}/src/m/qmu/examples/*.m \
					${ISSM_DIR}/src/m/kml/*.m \
					${ISSM_DIR}/src/m/dev/issmversion.m \
					${ISSM_DIR}/src/m/classes/*.m \
					${ISSM_DIR}/src/m/classes/qmu/*.m \
					${ISSM_DIR}/src/m/classes/qmu/@dakota_method/*.m \
					${ISSM_DIR}/src/m/classes/clusters/*.m \
					${ISSM_DIR}/src/m/consistency/*.m \
					${ISSM_DIR}/src/m/array/*.m \
					${ISSM_DIR}/src/m/boundaryconditions/*.m \
					${ISSM_DIR}/src/m/exp/*.m \
					${ISSM_DIR}/src/m/exp/operation/*.m \
					${ISSM_DIR}/src/m/geometry/*.m \
					${ISSM_DIR}/src/m/interp/*.m \
					${ISSM_DIR}/src/m/coordsystems/*.m \
					${ISSM_DIR}/src/m/mech/*.m \
					${ISSM_DIR}/src/m/mesh/*.m \
					${ISSM_DIR}/src/m/mesh/planet/spheretri/*.m \
					${ISSM_DIR}/src/m/mesh/planet/gmsh/*.m \
					${ISSM_DIR}/src/m/mesh/rifts/*.m \
					${ISSM_DIR}/src/m/miscellaneous/*.m \
					${ISSM_DIR}/src/m/modules/*.m \
					${ISSM_DIR}/src/m/os/*.m \
					${ISSM_DIR}/src/m/plot/*.m \
					${ISSM_DIR}/src/m/plot/colormaps/*.m \
					${ISSM_DIR}/src/m/string/*.m \
					${ISSM_DIR}/src/m/extrusion/*.m \
					${ISSM_DIR}/src/m/inversions/*.m \
					${ISSM_DIR}/src/m/io/*.m \
					${ISSM_DIR}/src/m/parameterization/*.m \
					${ISSM_DIR}/src/m/partition/*.m \
					${ISSM_DIR}/src/m/print/*.m \
					${ISSM_DIR}/src/m/regional/*.m \
					${ISSM_DIR}/src/m/solve/*.m \
					${ISSM_DIR}/src/m/solvers/*.m \
					${ISSM_DIR}/src/m/materials/*.m
endif
endif
if PYTHON
if !DEVELOPMENT
bin_SCRIPTS += ${ISSM_DIR}/src/m/archive/*.py \
					${ISSM_DIR}/src/m/classes/*.py \
					${ISSM_DIR}/src/m/classes/clusters/*.py \
					${ISSM_DIR}/src/m/consistency/*.py \
					${ISSM_DIR}/src/m/dev/issmversion.py \
					${ISSM_DIR}/src/m/boundaryconditions/*.py \
					${ISSM_DIR}/src/m/exp/*.py \
					${ISSM_DIR}/src/m/geometry/*.py \
					${ISSM_DIR}/src/m/coordsystems/*.py \
					${ISSM_DIR}/src/m/interp/*.py \
					${ISSM_DIR}/src/m/inversions/*.py \
					${ISSM_DIR}/src/m/mech/*.py \
					${ISSM_DIR}/src/m/mesh/*.py \
					${ISSM_DIR}/src/m/mesh/rifts/*.py \
					${ISSM_DIR}/src/m/miscellaneous/*.py \
					${ISSM_DIR}/src/m/modules/*.py \
					${ISSM_DIR}/src/m/os/*.py \
					${ISSM_DIR}/src/m/plot/*.py \
					${ISSM_DIR}/src/m/extrusion/*.py \
					${ISSM_DIR}/src/m/io/*.py \
					${ISSM_DIR}/src/m/parameterization/*.py \
					${ISSM_DIR}/src/m/solve/*.py \
					${ISSM_DIR}/src/m/solvers/*.py \
					${ISSM_DIR}/src/m/materials/*.py
endif
endif

if JAVASCRIPT
if !DEVELOPMENT

js_scripts =  ${ISSM_DIR}/src/m/array/arrayoperations.js \
				${ISSM_DIR}/src/m/boundaryconditions/SetIceShelfBC.js \
				${ISSM_DIR}/src/m/classes/autodiff.js \
				${ISSM_DIR}/src/m/classes/balancethickness.js \
				${ISSM_DIR}/src/m/classes/basalforcings.js \
				${ISSM_DIR}/src/m/classes/calving.js \
				${ISSM_DIR}/src/m/classes/clusters/generic.js \
				${ISSM_DIR}/src/m/classes/clusters/local.js \
				${ISSM_DIR}/src/m/classes/constants.js \
				${ISSM_DIR}/src/m/classes/damage.js \
				${ISSM_DIR}/src/m/classes/debug.js \
				${ISSM_DIR}/src/m/classes/flowequation.js \
				${ISSM_DIR}/src/m/classes/friction.js \
				${ISSM_DIR}/src/m/classes/geometry.js \
				${ISSM_DIR}/src/m/classes/giaivins.js \
				${ISSM_DIR}/src/m/classes/esa.js \
				${ISSM_DIR}/src/m/classes/groundingline.js \
				${ISSM_DIR}/src/m/classes/hydrologyshreve.js \
				${ISSM_DIR}/src/m/classes/initialization.js \
				${ISSM_DIR}/src/m/classes/inversion.js \
				${ISSM_DIR}/src/m/classes/mask.js \
				${ISSM_DIR}/src/m/classes/maskpsl.js \
				${ISSM_DIR}/src/m/classes/masstransport.js \
				${ISSM_DIR}/src/m/classes/matice.js \
				${ISSM_DIR}/src/m/classes/mesh2d.js \
				${ISSM_DIR}/src/m/classes/mesh3dsurface.js \
				${ISSM_DIR}/src/m/classes/miscellaneous.js \
				${ISSM_DIR}/src/m/classes/model.js \
				${ISSM_DIR}/src/m/classes/outputdefinition.js \
				${ISSM_DIR}/src/m/classes/pairoptions.js \
				${ISSM_DIR}/src/m/classes/plotoptions.js \
				${ISSM_DIR}/src/m/classes/priv.js \
				${ISSM_DIR}/src/m/classes/qmu.js \
				${ISSM_DIR}/src/m/classes/radaroverlay.js \
				${ISSM_DIR}/src/m/classes/rifts.js \
				${ISSM_DIR}/src/m/classes/slr.js \
				${ISSM_DIR}/src/m/classes/settings.js \
				${ISSM_DIR}/src/m/classes/SMBforcing.js \
				${ISSM_DIR}/src/m/classes/steadystate.js \
				${ISSM_DIR}/src/m/classes/stressbalance.js \
				${ISSM_DIR}/src/m/classes/template.js \
				${ISSM_DIR}/src/m/classes/thermal.js \
				${ISSM_DIR}/src/m/classes/timestepping.js \
				${ISSM_DIR}/src/m/classes/toolkits.js \
				${ISSM_DIR}/src/m/classes/levelset.js \
				${ISSM_DIR}/src/m/classes/trans.js \
				${ISSM_DIR}/src/m/classes/verbose.js \
				${ISSM_DIR}/src/m/consistency/checkfield.js \
				${ISSM_DIR}/src/m/consistency/ismodelselfconsistent.js \
				${ISSM_DIR}/src/m/geometry/FlagElements.js \
				${ISSM_DIR}/src/m/inversions/marshallcostfunctions.js \
				${ISSM_DIR}/src/m/io/fileptr.js \
				${ISSM_DIR}/src/m/io/savemodel.js \
				${ISSM_DIR}/src/m/io/loadmodel.js \
				${ISSM_DIR}/src/m/io/saveAsFile.js \
				${ISSM_DIR}/src/m/io/writetofile.js \
				${ISSM_DIR}/src/m/io/download.js \
				${ISSM_DIR}/src/m/io/writetofile.js \
				${ISSM_DIR}/src/m/materials/paterson.js \
				${ISSM_DIR}/src/m/mesh/triangle.js \
				${ISSM_DIR}/src/m/miscellaneous/colorbars.js \
				${ISSM_DIR}/src/m/miscellaneous/fielddisplay.js \
				${ISSM_DIR}/src/m/miscellaneous/rgbcolor.js \
				${ISSM_DIR}/src/m/parameterization/setflowequation.js \
				${ISSM_DIR}/src/m/parameterization/setmask.js \
				${ISSM_DIR}/src/m/plot/slider.js \
				${ISSM_DIR}/src/m/plot/applyoptions.js \
				${ISSM_DIR}/src/m/plot/checkplotoptions.js \
				${ISSM_DIR}/src/m/plot/plot_manager.js \
				${ISSM_DIR}/src/m/plot/plot_mesh.js \
				${ISSM_DIR}/src/m/plot/plot_overlay.js \
				${ISSM_DIR}/src/m/plot/plot_unit.js \
				${ISSM_DIR}/src/m/plot/plot_quiver.js \
				${ISSM_DIR}/src/m/plot/plotdoc.js \
				${ISSM_DIR}/src/m/plot/plotmodel.js \
				${ISSM_DIR}/src/m/plot/processdata.js \
				${ISSM_DIR}/src/m/plot/processmesh.js \
				${ISSM_DIR}/src/m/plot/gauge.js \
				${ISSM_DIR}/src/m/plot/webgl.js \
				${ISSM_DIR}/src/m/print/sprintf.js \
				${ISSM_DIR}/src/m/solve/loadresultsfrombuffer.js \
				${ISSM_DIR}/src/m/solve/marshall.js \
				${ISSM_DIR}/src/m/solve/parseresultsfrombuffer.js \
				${ISSM_DIR}/src/m/solve/solve.js \
				${ISSM_DIR}/src/m/solve/WriteData.js \
				${ISSM_DIR}/src/m/solvers/issmgslsolver.js 
bin_SCRIPTS +=  issm-bin.js

issm-bin.js: ${js_scripts}
	cat ${js_scripts}  > issm-bin.js

endif
endif

endif
