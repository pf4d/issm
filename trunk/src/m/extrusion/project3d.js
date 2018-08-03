project3d = function() {
    //PROJECT3D - vertically project a vector from 2d mesh
    //
    //   vertically project a vector from 2d mesh (split in noncoll and coll areas) into a 3d mesh.
    //   This vector can be a node vector of size (md.mesh.numberofvertices2d,N/A) or an 
    //   element vector of size (md.mesh.numberofelements2d,N/A). 
    //   arguments: 
    //      'vector': 2d vector
    //      'type': 'element' or 'node'. 
    //   options: 
    //      'layer' a layer number where vector should keep its values. If not specified, all layers adopt the 
    //             value of the 2d vector.
    //      'padding': default to 0 (value adopted by other 3d layers not being projected0
    //
    //   Egs:
    // md.extruded_vector=project3d(md,'vector',vector2d,'type','node','layer',1,'padding',null);
    // md.extruded_vector=project3d(md,'vector',vector2d,'type','element','padding',0);
    // md.extruded_vector=project3d(md,'vector',vector2d,'type','node');

    //some regular checks
    
    function remove_first_n(start, arglist) { // Using slice() on arguments is discouraged because it prevents optimizations in engines such as V8
        var args = [];

        for (var i = start; i < arglist.length; i++) {
            args.push(arglist[i]);
        }

        return args;
    }

    if (arguments.length===1 || arguments.length===0) {
        //help project3d
        console.error('project3d bad usage');
    }

    var md = arguments[0];

    if (md.mesh.elementtype() !== 'Penta') {
        console.error('input model is not 3d');
    }

    //retrieve parameters from options.
    options      = new pairoptions(remove_first_n(1, arguments)); // slice to remove md
    vector2d     = options.getfieldvalue('vector');     //mandatory
    type         = options.getfieldvalue('type');       //mandatory
    layer        = options.getfieldvalue('layer',0);    //optional (do all layers default:)
    paddingvalue = options.getfieldvalue('padding',0);  //0 by default

    if (isNaN(vector2d) || vector2d === 0 || vector2d.length === 1) { // NaN treated as length 1 in MATLAB
        projected_vector=vector2d;
    } else if (type.toLowerCase() === 'node') {

        //Initialize 3d vector
        if (vector2d.length===md.mesh.numberofvertices2d) {
                projected_vector=ones(md.mesh.numberofvertices,  vector2d[0].length).map(function(arr) { return arr.fill(paddingvalue); });
        } else if (vector2d.length===md.mesh.numberofvertices2d+1) {
            projected_vector=ones(md.mesh.numberofvertices+1,vector2d[0].length).map(function(arr) { return arr.fill(paddingvalue); });
            projected_vector[projected_vector.length-1] = projected_vector[projected_vector.length-1].map(function(element, index) { return vector2d[vector2d.length-1][index]; });
            vector2d.pop(); // Remove last array
        } else {
            console.error('vector length not supported')
        }

            //Fill in
        if (layer===0) {
            for (var i = 1; i < md.mesh.numberoflayers; ++i) {
                var vec_idx = 0;

                for (var j = (i-1)*md.mesh.numberofvertices2d+1, vec_idx = 0; j <= i * md.mesh.numberofvertices2d && vec_idx < vector2d.length; ++j, ++vec_idx) {
                    projected_vector[j].map(function(element, index) { return vector2d[vec_idx][index]; });
                }

            }
        } else {
            var vec_idx = 0;

            for (var i = (layer-1)*md.mesh.numberofvertices2d+1, vec_idx = 0; i <= layer * md.mesh.numberofvertices2d && vec_idx < vector2d.length; ++i, ++vec_idx) {
                projected_vector[i] = vector2d[vec_idx];
            }
        }
    } else if (type.toLowerCase() === 'element') {
            //Initialize 3d vector
        if (vector2d.length===md.mesh.numberofelements2d) {
                projected_vector = ones(md.mesh.numberofelements,  vector2d[0].length).map(function(arr) { return arr.fill(paddingvalue); });
        } else if (vector2d.length===md.mesh.numberofelements2d+1) {
            projected_vector = ones(md.mesh.numberofelements+1,  vector2d[0].length).map(function(arr) { return arr.fill(paddingvalue); });

            projected_vector[projected_vector.length-1].map(function(element, index) { return vector2d[vector2d.length-1][index]; });

            vector2d.pop();
        } else {
            console.error('vector length not supported')
        }

            //Fill in
        if (layer===0) {
            for (var i = 1; i < md.mesh.numberoflayers-1; ++i) {
                var vec_idx=0;
                for (var j = (i-1)*md.mesh.numberofelements2d+1, vec_idx = 0; j <= i * md.mesh.numberofelements2d && vec_idx < vector2d.length; ++j, ++vec_idx) {
                    projected_vector[j].map(function(element, index) { return vector2d[vec_idx][index]; });
                }
            }
        } else {
            var vec_idx=0;
            for (var i = (layer-1)*md.mesh.numberofelements2d+1, vec_idx = 0; i <= layer * md.mesh.numberofelements2d && vec_idx < vector2d.length; ++i, ++vec_idx) {
                projected_vector[i].map(function(element, index) { return vector2d[vec_idx][index]; });
            }
        }
    } else {
        console.error('project3d error message: unknown projection type');
    }

    return projected_vector;
};
