//MODEL class definition
//
//   Usage:
//      md = new model()

function model () {
	//methods
		this.disp = function() { //{{{
			console.log(sprintf("class model echo: "));
			console.log(sprintf("//19s: //-22s -- //s","mesh"            ,"[1x1 " + typeof(this.mesh) + "]","mesh properties"));
			console.log(sprintf("//19s: //-22s -- //s","mask"            ,"[1x1 " + typeof(this.mask) + "]","defines grounded and floating elements"));
			console.log(sprintf("//19s: //-22s -- //s","geometry"        ,"[1x1 " + typeof(this.geometry) + "]","surface elevation, bedrock topography, ice thickness,..."));
			console.log(sprintf("//19s: //-22s -- //s","constants"       ,"[1x1 " + typeof(this.constants) + "]","physical constants"));
			console.log(sprintf("//19s: //-22s -- //s","smb"             ,"[1x1 " + typeof(this.smb) + "]","surface mass balance"));
			console.log(sprintf("//19s: //-22s -- //s","basalforcings"   ,"[1x1 " + typeof(this.basalforcings) + "]","bed forcings"));
			console.log(sprintf("//19s: //-22s -- //s","materials"       ,"[1x1 " + typeof(this.materials) + "]","material properties"));
			console.log(sprintf("//19s: //-22s -- //s","damage"          ,"[1x1 " + typeof(this.damage) + "]","parameters for damage evolution solution"));
			console.log(sprintf("//19s: //-22s -- //s","friction"        ,"[1x1 " + typeof(this.friction) + "]","basal friction/drag properties"));
			console.log(sprintf("//19s: //-22s -- //s","flowequation"    ,"[1x1 " + typeof(this.flowequation) + "]","flow equations"));
			console.log(sprintf("//19s: //-22s -- //s","timestepping"    ,"[1x1 " + typeof(this.timestepping) + "]","time stepping for trans models"));
			console.log(sprintf("//19s: //-22s -- //s","initialization"  ,"[1x1 " + typeof(this.initialization) + "]","initial guess/state"));
			console.log(sprintf("//19s: //-22s -- //s","rifts"           ,"[1x1 " + typeof(this.rifts) + "]","rifts properties"));
			console.log(sprintf("//19s: //-22s -- //s","slr"             ,"[1x1 " + typeof(this.slr) + "]","slr forcings"));
			console.log(sprintf("//19s: //-22s -- //s","debug"           ,"[1x1 " + typeof(this.debug) + "]","debugging tools (valgrind, gprof)"));
			console.log(sprintf("//19s: //-22s -- //s","verbose"         ,"[1x1 " + typeof(this.verbose) + "]","verbosity level in solve"));
			console.log(sprintf("//19s: //-22s -- //s","settings"        ,"[1x1 " + typeof(this.settings) + "]","settings properties"));
			console.log(sprintf("//19s: //-22s -- //s","toolkits"        ,"[1x1 " + typeof(this.toolkits) + "]","PETSc options for each solution"));
			console.log(sprintf("//19s: //-22s -- //s","cluster"         ,"[1x1 " + typeof(this.cluster) + "]","cluster parameters (number of cpus...)"));
			console.log(sprintf("//19s: //-22s -- //s","balancethickness","[1x1 " + typeof(this.balancethickness) + "]","parameters for balancethickness solution"));
			console.log(sprintf("//19s: //-22s -- //s","stressbalance"   ,"[1x1 " + typeof(this.stressbalance) + "]","parameters for stressbalance solution"));
			console.log(sprintf("//19s: //-22s -- //s","groundingline"   ,"[1x1 " + typeof(this.groundingline) + "]","parameters for groundingline solution"));
			console.log(sprintf("//19s: //-22s -- //s","hydrology"       ,"[1x1 " + typeof(this.hydrology) + "]","parameters for hydrology solution"));
			console.log(sprintf("//19s: //-22s -- //s","masstransport"   ,"[1x1 " + typeof(this.masstransport) + "]","parameters for masstransport solution"));
			console.log(sprintf("//19s: //-22s -- //s","thermal"         ,"[1x1 " + typeof(this.thermal) + "]","parameters for thermal solution"));
			console.log(sprintf("//19s: //-22s -- //s","steadystate"     ,"[1x1 " + typeof(this.steadystate) + "]","parameters for steadystate solution"));
			console.log(sprintf("//19s: //-22s -- //s","trans"       ,"[1x1 " + typeof(this.trans) + "]","parameters for trans solution"));
			console.log(sprintf("//19s: //-22s -- //s","levelset"        ,"[1x1 " + typeof(this.levelset) + "]","parameters for moving boundaries (level-set method)"));
			console.log(sprintf("//19s: //-22s -- //s","calving"         ,"[1x1 " + typeof(this.calving) + "]","parameters for calving"));
			console.log(sprintf("//19s: //-22s -- //s","gia"             ,"[1x1 " + typeof(this.gia) + "]","parameters for gia solution"));
			console.log(sprintf("//19s: //-22s -- //s","autodiff"        ,"[1x1 " + typeof(this.autodiff) + "]","automatic differentiation parameters"));
			console.log(sprintf("//19s: //-22s -- //s","inversion"       ,"[1x1 " + typeof(this.inversion) + "]","parameters for inverse methods"));
			console.log(sprintf("//19s: //-22s -- //s","qmu"             ,"[1x1 " + typeof(this.qmu) + "]","dakota properties"));
			console.log(sprintf("//19s: //-22s -- //s","amr"             ,"[1x1 " + typeof(this.amr) + "]","adaptive mesh refinement properties"));
			console.log(sprintf("//19s: //-22s -- //s","outputdefinition","[1x1 " + typeof(this.outputdefinition) + "]","output definition"));
			console.log(sprintf("//19s: //-22s -- //s","results"         ,"[1x1 " + typeof(this.results) + "]","model results"));
			console.log(sprintf("//19s: //-22s -- //s","radaroverlay"    ,"[1x1 " + typeof(this.radaroverlay) + "]","radar image for plot overlay"));
			console.log(sprintf("//19s: //-22s -- //s","miscellaneous"   ,"[1x1 " + typeof(this.miscellaneous) + "]","miscellaneous fields"));
		} //}}}
		this.setdefaultparameters = function () { // {{{

			//initialize subclasses
			this.mesh             = new mesh2d();
			this.mask             = new mask();
			this.constants        = new constants();
			this.geometry         = new geometry();
			this.initialization   = new initialization();
			this.smb              = new SMBforcing();
			this.basalforcings    = new basalforcings();
			this.friction         = new friction();
			this.rifts            = new rifts();
			this.slr              = new slr();
			this.timestepping     = new timestepping();
			this.groundingline    = new groundingline();
			this.materials        = new matice();
			this.damage           = new damage();
			this.flowequation     = new flowequation();
			this.debug            = new debug();
			this.verbose          = new verbose();
			this.settings         = new settings();
			this.toolkits         = new toolkits();
			this.cluster          = new local();
			this.balancethickness = new balancethickness();
			this.stressbalance    = new stressbalance();
			this.hydrology        = new hydrologyshreve();
			this.masstransport    = new masstransport();
			this.thermal          = new thermal();
			this.steadystate      = new steadystate();
			this.trans            = new trans();
			this.levelset		  = new levelset();
			this.calving          = new calving();
			this.gia              = new giaivins();
			this.autodiff         = new autodiff();
			this.inversion        = new inversion();
			this.qmu              = new qmu();
			this.amr					 = new amr();
			this.radaroverlay     = new radaroverlay();
			this.results          = {};
			this.outputdefinition = new outputdefinition();
			this.miscellaneous    = new miscellaneous();
			this.priv             = new priv();
		} //}}}
		this.checkmessage = function(string){ //{{{
			console.log('model not consistent: ' + string);
			this.priv.isconsistent=false;
		} //}}}
		this.fix = function(){ //{{{

			for (var field in this){

				//Some properties do not need to be fixed
				if (field == 'results' | field =='radaroverlay' | field == 'toolkits' | field =='cluster' | field == 'priv') continue;

				//Check that current field is a class
				if(typeof this[field] == 'function'){
					continue;
				}

				//Fix current object
				this[field].fix(this);
			}

		} //}}}
        this.extrude = function(md) { //{{{
        //EXTRUDE - vertically extrude a 2d mesh
        //
        //   vertically extrude a 2d mesh and create corresponding 3d mesh.
        //   The vertical distribution can:
        //    - follow a polynomial law
        //    - follow two polynomial laws, one for the lower part and one for the upper part of the mesh
        //    - be discribed by a list of coefficients (between 0 and 1)
        //   
        //
        //   Usage:
        //      md=extrude(md,numlayers,extrusionexponent);
        //      md=extrude(md,numlayers,lowerexponent,upperexponent);
        //      md=extrude(md,listofcoefficients);
        //
        //   Example:
        //      md=extrude(md,15,1.3);
        //      md=extrude(md,15,1.3,1.2);
        //      md=extrude(md,[0 0.2 0.5 0.7 0.9 0.95 1]);
        //
        //   See also: MODELEXTRACT, COLLAPSE

        //some checks on list of arguments
            var argc = arguments.length;

            if ((argc > 4) | (argc < 2)) {
                throw "extrude error message";
            }

            //Extrude the mesh
            if (argc==2) { //list of coefficients
                clist=arguments[0];

                if (clist.some(function(x) { return x < 0; }) | clist.some(function(x) { return x > 1; }))
                    throw 'extrusioncoefficients must be between 0 and 1';
                clistcopy.push.apply(clistcopy, [0, 1]);

                extrusionlist=clistcopy.reverse().filter(function(e, i, arr) { //only keep unique elements
                    return arr.indexOf(e, i+1) === -1;
                }).sort();

                numlayers=extrusionlist.length;
            } else if (argc==3) { //one polynomial law
                if (arguments[1]<=0) {
                    //help extrude;
                    throw 'extrusionexponent must be >=0';
                }
                numlayers=arguments[0];

                extrusionlist = [];

                for (var i = 0; i < numlayers; i++) {
                    extrusionlist.push(Math.pow(i/(numlayers-1), arguments[1]));
                };

            } else if (argc==4) { //two polynomial laws
                numlayers=arguments[0];
                lowerexp=arguments[1];
                upperexp=arguments[2];

                if (arguments[1]<=0 | args[2]<=0) {
                    throw 'lower and upper extrusionexponents must be >=0';
                }

                lowerextrusionlist = [];
                for (var i = 0; i <= 1; i += (2/(numlayers-1))) {
                    lowerextrusionlist.push(Math.pow(i, lowerexp)/2)
                };

                upperextrusionlist = [];
                for (var i = 0; i <= 1; i += (2/(numlayers-1))) {
                    upperextrusionlist.push(Math.pow(i, upperexp)/2)
                };

                //lowerextrusionlist=[(0:2/(numlayers-1):1).^lowerexp]/2;
                //upperextrusionlist=[(0:2/(numlayers-1):1).^upperexp]/2;

                var temp = lowerextrusionlist;
                temp.push.apply(temp, [1]);

                var temp2 = upperextrusionlist.map(function(x) { return -x; });
                temp.push.apply(temp, temp2);

                clistcopy.push.apply(clistcopy, [0, 1]);
                //extrusionlist=sort(unique([lowerextrusionlist 1-upperextrusionlist]));
                extrusionlist=temp.sort();
            }

            if (numlayers<2) {
                console.error('number of layers should be at least 2');
            }

            if (md.mesh.domaintype() === '3D') {
                console.error('Cannot extrude a 3d mesh (extrude cannot be called more than once)');
            }

            //Initialize with the 2d mesh
            var mesh2d = md.mesh;
            md.mesh=new mesh3dprisms();

            md.mesh.x                           = mesh2d.x;
            md.mesh.y                           = mesh2d.y;
            md.mesh.elements                    = mesh2d.elements;
            md.mesh.numberofelements            = mesh2d.numberofelements;
            md.mesh.numberofvertices            = mesh2d.numberofvertices;

            md.mesh.lat                         = mesh2d.lat;
            md.mesh.long                        = mesh2d.long;
            md.mesh.epsg                        = mesh2d.epsg;

            md.mesh.vertexonboundary            = mesh2d.vertexonboundary;
            md.mesh.vertexconnectivity          = mesh2d.vertexconnectivity;
            md.mesh.elementconnectivity         = mesh2d.elementconnectivity;
            md.mesh.average_vertex_connectivity = mesh2d.average_vertex_connectivity;

            md.mesh.extractedvertices           = mesh2d.extractedvertices;
            md.mesh.extractedelements           = mesh2d.extractedelements;

            x3d=[]; 
            y3d=[];
            z3d=[];  //the lower node is on the bed
            thickness3d=md.geometry.thickness; //thickness and bed for these nodes
                bed3d=md.geometry.base;

            //Create the new layers
            for (var i = 1; i <= numlayers; i++) {
                x3d=[x3d, md.mesh.x]; 
                y3d=[y3d, md.mesh.y];
                //nodes are distributed between bed and surface accordingly to the given exponent
                z3d=[z3d, bed3d+thickness3d*extrusionlist(i)]; 
            };

            number_nodes3d=x3d.length; //number of 3d nodes for the non extruded part of the mesh

            //Extrude elements 
            elements3d=[];
            for (var i = 0; i < numlayers; i++) {
                elements3d.push(md.mesh.elements+(i-1)*md.mesh.numberofvertices, md.mesh.elements+i*md.mesh.numberofvertices); //Create the elements of the 3d mesh for the non extruded part
            };

            number_el3d=elements3d.length; //number of 3d nodes for the non extruded part of the mesh

            function fillArray(value, len) {
                if (len == 0) return [];
                var a = [value];
                while (a.length * 2 <= len) a = a.concat(a);
                if (a.length < len) a = a.concat(a.slice(0, len - a.length));
                return a;
            }

            //Keep a trace of lower and upper nodes
            //lowervertex=NaN*ones(number_nodes3d,1);
            //uppervertex=NaN*ones(number_nodes3d,1);

            lowervertex = fillArray(NaN, number_nodes3d);
            uppervertex = fillArray(NaN, number_nodes3d);

            //lowervertex(md.mesh.numberofvertices+1:end)=1:(numlayers-1)*md.mesh.numberofvertices;
            //uppervertex(1:(numlayers-1)*md.mesh.numberofvertices)=md.mesh.numberofvertices+1:number_nodes3d;

            for (var i = md.mesh.numberofvertices+1, k = 1; i < lowervertex.length && k <= (numlayers-1)*md.mesh.numberofvertices; i++, k++) {
                lowervertex[i] = k;
            };

            for (var i = 1, k = md.mesh.numberofvertices+1; i <= (numlayers-1)*md.mesh.numberofvertices && k <= number_nodes3d; i++, k++) {
                uppervertex[i] = k;
            };

            md.mesh.lowervertex=lowervertex;
            md.mesh.uppervertex=uppervertex;

            //same for lower and upper elements
            lowerelements = fillArray(NaN, number_el3d);
            upperelements = fillArray(NaN, number_el3d);
            //lowerelements(md.mesh.numberofelements+1:end)=1:(numlayers-2)*md.mesh.numberofelements;
            //upperelements(1:(numlayers-2)*md.mesh.numberofelements)=md.mesh.numberofelements+1:(numlayers-1)*md.mesh.numberofelements;

            for (var i = md.mesh.numberofvertices+1, k = 1; i < lowerelements.length && k <= (numlayers-2)*md.mesh.numberofelements; i++, k++) {
                lowerelements[i] = k;
            };

            for (var i = 1, k = md.mesh.numberofelements + 1; i <= (numlayers-2)*md.mesh.numberofelements && k <= (numlayers-1)*md.mesh.numberofelements; i++, k++) {
                upperelements[i] = k;
            };

            md.mesh.lowerelements=lowerelements;
            md.mesh.upperelements=upperelements;

            //Save old mesh 
            md.mesh.x2d=md.mesh.x;
            md.mesh.y2d=md.mesh.y;
            md.mesh.elements2d=md.mesh.elements;
            md.mesh.numberofelements2d=md.mesh.numberofelements;
            md.mesh.numberofvertices2d=md.mesh.numberofvertices;

            //Build global 3d mesh 
            md.mesh.elements=elements3d;
            md.mesh.x=x3d;
            md.mesh.y=y3d;
            md.mesh.z=z3d;
            md.mesh.numberofelements=number_el3d;
            md.mesh.numberofvertices=number_nodes3d;
            md.mesh.numberoflayers=numlayers;

            //Ok, now deal with the other fields from the 2d mesh:

            //bedinfo and surface info
            md.mesh.vertexonbase=project3d(md,'vector',ones(md.mesh.numberofvertices2d,1),'type','node','layer',1);
            md.mesh.vertexonsurface=project3d(md,'vector',ones(md.mesh.numberofvertices2d,1),'type','node','layer',md.mesh.numberoflayers);
            md.mesh.vertexonboundary=project3d(md,'vector',md.mesh.vertexonboundary,'type','node');

            //lat long
            md.mesh.lat=project3d(md,'vector',md.mesh.lat,'type','node');
            md.mesh.long=project3d(md,'vector',md.mesh.long,'type','node');

                //md.geometry=md.geometry.extrude(md);
                //md.friction=md.friction.extrude(md);
                //md.inversion=md.inversion.extrude(md);
                //md.smb=md.smb.extrude(md);
                //md.initialization=md.initialization.extrude(md);

                //md.flowequation=md.flowequation.extrude(md);
                //md.stressbalance=md.stressbalance.extrude(md);
                //md.thermal=md.thermal.extrude(md);
                //md.masstransport=md.masstransport.extrude(md);
                //md.levelset=md.levelset.extrude(md);
                //md.calving=md.calving.extrude(md);
                //md.hydrology = md.hydrology.extrude(md);

                ////connectivity
                ////if ~isnan(md.mesh.elementconnectivity)

                //if (md.mesh.elementconnectivity.every(function(e, i, arr) { return e !== null && typeof e !== 'undefined'; })) {
                    //var temparr = [];
                    //for (var i = 0; i < numlayers; i++) {
                        //temparr.push(md.mesh.elementconnectivity);
                    //};

                    ////md.mesh.elementconnectivity=repmat(md.mesh.elementconnectivity,numlayers-1,1);

                    //md.mesh.elementconnectivity = temparr;

                    ////md.mesh.elementconnectivity(find(md.mesh.elementconnectivity==0))=NaN;

                    //for (var i = 0; i < md.mesh.elementconnectivity.length; i++) {
                        //if (md.mesh.elementconnectivity[i] == 0) {
                            //md.mesh.elementconnectivity[i] = NaN;
                        //}
                    //};

                    //for (var i = 2; i < numlayers; i++) {
                        //for (var j = (i-1)*md.mesh.nuberofelements2d+1; i <= (i)*md.mesh.numberofelements2d; i++) {
                            //md.mesh.elementconnectivity[i] += md.mesh.numberofelements2d;
                        //};
                        ////md.mesh.elementconnectivity((i-1)*md.mesh.numberofelements2d+1:(i)*md.mesh.numberofelements2d,:)...
                            ////=md.mesh.elementconnectivity((i-1)*md.mesh.numberofelements2d+1:(i)*md.mesh.numberofelements2d,:)+md.mesh.numberofelements2d;
                    //};

                    //md.mesh.elementconnectivity.map(function(x) { return (isNaN(x)) ? 0 : x; });
                //}

                //md.materials=md.materials.extrude(md);
                //md.damage=md.damage.extrude(md);
                //md.mask=md.mask.extrude(md);
                //md.qmu=md.qmu.extrude(md);
                //md.basalforcings=md.basalforcings.extrude(md);

                ////increase connectivity if less than 25:
                //if (md.mesh.average_vertex_connectivity<=25)
                    //md.mesh.average_vertex_connectivity=100;
        } //}}}
        this.collapse = function(md) { //{{{
            /*
             *COLLAPSE - collapses a 3d mesh into a 2d mesh
             *
             *   This routine collapses a 3d model into a 2d model
             *   and collapses all the fileds of the 3d model by
             *   taking their depth-averaged values
             *
             *   Usage:
             *     md=collapse(md)
             *
             *   See also: EXTRUDE, MODELEXTRACT
             */

            // Check that the model is really a 3d model
            if (md.mesh.elementtype() !== 'Penta') {
                console.error('collapse error message: only 3d mesh can be collapsed')
            }

            // Start with changing all the fields from the 3d mesh 

            // dealing with the friction law
            // drag is limited to nodes that are on the bedrock.
            if (md.friction.classname() === 'friction') {
                md.friction.coefficient=project2d(md,md.friction.coefficient,1);
                md.friction.p=project2d(md,md.friction.p,1);
                md.friction.q=project2d(md,md.friction.q,1);
            } else if (md.friction.classname() === 'frictionhydro') {
                md.friction.q=project2d(md,md.friction.q,1);
                md.friction.C=project2d(md,md.friction.C,1);
                md.friction.As=project2d(md,md.friction.As,1);
                md.friction.effective_pressure=project2d(md,md.friction.effective_pressure,1);
            } else if (md.friction.classname() === 'frictionwaterlayer') {
                md.friction.coefficient=project2d(md,md.friction.coefficient,1);
                md.friction.p=project2d(md,md.friction.p,1);
                md.friction.q=project2d(md,md.friction.q,1);
                md.friction.water_layer=project2d(md,md.friction.water_layer,1);
            } else if (md.friction.classname() === 'frictionweertman') {
                md.friction.C=project2d(md,md.friction.C,1);
                md.friction.m=project2d(md,md.friction.m,1);
            } else if (md.friction.classname() === 'frictionweertmantemp') {
                md.friction.C=project2d(md,md.friction.C,1);
                md.friction.m=project2d(md,md.friction.m,1);
            } else {
                disp('friction type not supported');
            }

            // observations
            if (!isNaN(md.inversion.vx_obs))
                md.inversion.vx_obs=project2d(md,md.inversion.vx_obs,md.mesh.numberoflayers);

            if (!isNaN(md.inversion.vy_obs))
                md.inversion.vy_obs=project2d(md,md.inversion.vy_obs,md.mesh.numberoflayers);

            if (!isNaN(md.inversion.vel_obs))
                md.inversion.vel_obs=project2d(md,md.inversion.vel_obs,md.mesh.numberoflayers);

            if (!isNaN(md.inversion.cost_functions_coefficients))
                md.inversion.cost_functions_coefficients=project2d(md,md.inversion.cost_functions_coefficients,md.mesh.numberoflayers);

            if (numel(md.inversion.min_parameters)>1)
                md.inversion.min_parameters=project2d(md,md.inversion.min_parameters,md.mesh.numberoflayers);

            if (numel(md.inversion.max_parameters)>1) 
                md.inversion.max_parameters=project2d(md,md.inversion.max_parameters,md.mesh.numberoflayers);

            if (md.smb.classname() === 'SMBforcing' && !isNaN(md.smb.mass_balance)) {
                md.smb.mass_balance=project2d(md,md.smb.mass_balance,md.mesh.numberoflayers);
            } else if (md.smb.classname() === 'SMBhenning' && !isNaN(md.smb.smbref)) {
                md.smb.smbref=project2d(md,md.smb.smbref,md.mesh.numberoflayers);
            }

            // results
            if (!isNaN(md.initialization.vx))
                md.initialization.vx=DepthAverage(md,md.initialization.vx);
            if (!isNaN(md.initialization.vy))
                md.initialization.vy=DepthAverage(md,md.initialization.vy);
            if (!isNaN(md.initialization.vz))
                md.initialization.vz=DepthAverage(md,md.initialization.vz);
            if (!isNaN(md.initialization.vel))
                md.initialization.vel=DepthAverage(md,md.initialization.vel);
            if (!isNaN(md.initialization.temperature))
                md.initialization.temperature=DepthAverage(md,md.initialization.temperature);
            if (!isNaN(md.initialization.pressure))
                md.initialization.pressure=project2d(md,md.initialization.pressure,1);
            if (!isNaN(md.initialization.sediment_head))
                md.initialization.sediment_head=project2d(md,md.initialization.sediment_head,1);
            if (!isNaN(md.initialization.epl_head))
                md.initialization.epl_head=project2d(md,md.initialization.epl_head,1);
            if (!isNaN(md.initialization.epl_thickness))
                md.initialization.epl_thickness=project2d(md,md.initialization.epl_thickness,1);

            // giaivins
            if (!isNaN(md.gia.mantle_viscosity))
                md.gia.mantle_viscosity=project2d(md,md.gia.mantle_viscosity,1);
            if (!isNaN(md.gia.lithosphere_thickness))
                md.gia.lithosphere_thickness=project2d(md,md.gia.lithosphere_thickness,1);

            // elementstype
            if (!isNaN(md.flowequation.element_equation)) {
                md.flowequation.element_equation=project2d(md,md.flowequation.element_equation,1);
                md.flowequation.vertex_equation=project2d(md,md.flowequation.vertex_equation,1);
                md.flowequation.borderSSA=project2d(md,md.flowequation.borderSSA,1);
                md.flowequation.borderHO=project2d(md,md.flowequation.borderHO,1);
                md.flowequation.borderFS=project2d(md,md.flowequation.borderFS,1);
            }

            // boundary conditions
            md.stressbalance.spcvx=project2d(md,md.stressbalance.spcvx,md.mesh.numberoflayers);
            md.stressbalance.spcvy=project2d(md,md.stressbalance.spcvy,md.mesh.numberoflayers);
            md.stressbalance.spcvz=project2d(md,md.stressbalance.spcvz,md.mesh.numberoflayers);
            md.stressbalance.referential=project2d(md,md.stressbalance.referential,md.mesh.numberoflayers);
            md.stressbalance.loadingforce=project2d(md,md.stressbalance.loadingforce,md.mesh.numberoflayers);
            md.masstransport.spcthickness=project2d(md,md.masstransport.spcthickness,md.mesh.numberoflayers);
            if (!isNaN(md.damage.spcdamage))
                md.damage.spcdamage=project2d(md,md.damage.spcdamage,md.mesh.numberoflayers);
            md.thermal.spctemperature=project2d(md,md.thermal.spctemperature,md.mesh.numberoflayers);

            // Hydrologydc variables
            if (md.hydrology.classname() === 'hydrologydc') {
                md.hydrology.spcsediment_head=project2d(md,md.hydrology.spcsediment_head,1);
                md.hydrology.mask_eplactive_node=project2d(md,md.hydrology.mask_eplactive_node,1);
                md.hydrology.sediment_transmitivity=project2d(md,md.hydrology.sediment_transmitivity,1);
                md.hydrology.basal_moulin_input=project2d(md,md.hydrology.basal_moulin_input,1);

                if(md.hydrology.isefficientlayer==1)
                    md.hydrology.spcepl_head=project2d(md,md.hydrology.spcepl_head,1);
            }
            
            // materials
            md.materials.rheology_B=DepthAverage(md,md.materials.rheology_B);
            md.materials.rheology_n=project2d(md,md.materials.rheology_n,1);
            
            // damage: 
            if (md.damage.isdamage)
                md.damage.D=DepthAverage(md,md.damage.D);

            // special for thermal modeling:
            if (!isNaN(md.basalforcings.groundedice_melting_rate))
                md.basalforcings.groundedice_melting_rate=project2d(md,md.basalforcings.groundedice_melting_rate,1); 

            if (!isNaN(md.basalforcings.floatingice_melting_rate))
                md.basalforcings.floatingice_melting_rate=project2d(md,md.basalforcings.floatingice_melting_rate,1); 

            md.basalforcings.geothermalflux=project2d(md,md.basalforcings.geothermalflux,1); // bedrock only gets geothermal flux

            // update of connectivity matrix
            md.mesh.average_vertex_connectivity=25;

            // Collapse the mesh
            var nodes2d=md.mesh.numberofvertices2d;
            var elements2d=md.mesh.numberofelements2d;

            // parameters
            md.geometry.surface=project2d(md,md.geometry.surface,1);
            md.geometry.thickness=project2d(md,md.geometry.thickness,1);
            md.geometry.base=project2d(md,md.geometry.base,1);
            if (!isNaN(md.geometry.bed))
                md.geometry.bed=project2d(md,md.geometry.bed,1);

            if (!isNaN(md.mask.groundedice_levelset))
                md.mask.groundedice_levelset=project2d(md,md.mask.groundedice_levelset,1);

            if (!isNaN(md.mask.ice_levelset))
                md.mask.ice_levelset=project2d(md,md.mask.ice_levelset,1);

            // lat long
            if (numel(md.mesh.lat) === md.mesh.numberofvertices)
                md.mesh.lat=project2d(md,md.mesh.lat,1);
            if (numel(md.mesh.long) === md.mesh.numberofvertices)
                md.mesh.long=project2d(md,md.mesh.long,1);

            // Initialize with the 2d mesh
            var mesh = new mesh2d();
            mesh.x=md.mesh.x2d;
            mesh.y=md.mesh.y2d;
            mesh.numberofvertices=md.mesh.numberofvertices2d;
            mesh.numberofelements=md.mesh.numberofelements2d;
            mesh.elements=md.mesh.elements2d;

            if (!isNaN(md.mesh.vertexonboundary))
                mesh.vertexonboundary=project2d(md,md.mesh.vertexonboundary,1);
            if (!isNaN(md.mesh.elementconnectivity))
                mesh.elementconnectivity=project2d(md,md.mesh.elementconnectivity,1);

            md.mesh=mesh;
            md.mesh.vertexconnectivity=NodeConnectivity(md.mesh.elements,md.mesh.numberofvertices);
            md.mesh.elementconnectivity=ElementConnectivity(md.mesh.elements,md.mesh.vertexconnectivity);
            md.mesh.segments=contourenvelope(md.mesh);

            return md;
        } /*}}}*/
	//properties
	// {{{
		//Careful here: no other class should be used as default value this is a bug of matlab
		this.mesh             = 0;
		this.mask             = 0;

		this.geometry         = 0;
		this.constants        = 0;
		this.smb              = 0;
		this.basalforcings    = 0;
		this.materials        = 0;
		this.damage           = 0;
		this.friction         = 0;
		this.flowequation     = 0;
		this.timestepping     = 0;
		this.initialization   = 0;
		this.rifts            = 0;
		this.slr              = 0;

		this.debug            = 0;
		this.verbose          = 0;
		this.settings         = 0;
		this.toolkits         = 0;
		this.cluster          = 0;

		this.balancethickness = 0;
		this.stressbalance    = 0;
		this.groundingline    = 0;
		this.hydrology        = 0;
		this.masstransport    = 0;
		this.thermal          = 0;
		this.steadystate      = 0;
		this.trans            = 0;
		this.levelset         = 0;
		this.calving          = 0;
		this.gia              = 0;

		this.autodiff         = 0;
		this.inversion        = 0;
		this.qmu              = 0;
		this.amr					 = 0;

		this.results          = 0;
		this.outputdefinition = 0;
		this.radaroverlay     = 0;
		this.miscellaneous    = 0;
		this.priv             = 0;

		//set default values for fields
		this.setdefaultparameters();

		//}}}
}
