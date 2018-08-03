function plotdoc()
%PLOTDOC - plot documentation
%
%   Usage:
%      plotdoc()

disp(' ');
disp('   Plot usage: plotm(model,varargin)');
disp('   Options: ');
disp('       ''figure'': figure number');
disp('       ''data'' : what we want to plot');
disp('                Available values for ''data'' are: ');
disp('                  - any field of the model structure. ex: plot(md,''data'',''vel''), or plot(md,''data'',md.initialization.vel)');
disp('                  - ''basal_drag'': plot the basal drag on the bed (in kPa) based on the velocity in md.initialization');
disp('                  - ''basal_dragx'' or ''basal_dragy'' : plot a component of the basal drag on the bed (in kPa)');
disp('                  - ''boundaries'': this will draw all the segment boundaries to the model, including rifts.');
disp('                  - ''icefront'': this will show segments that are used to define the icefront of the model (Neumann boundary conditions).');
disp('                  - ''BC'': this will draw all the boundary conditions (Dirichlet and Neumann).');
disp('                  - ''deviatoricstress_tensor'': plot the components of the deviatoric stress tensor (tauxx,tauyy,tauzz,tauxy,tauxz,tauyz) if computed');
disp('                  - ''deviatoricstress_principal'': plot the deviatoricstress tensor principal axis and principal values');
disp('                  - ''deviatoricstress_principalaxis1'': arrow plot the first principal axis of the deviatoricstress tensor(replace 1 by 2 or 3 if needed)');
disp('                  - ''driving_stress'': plot the driving stress (in kPa)');
disp('                  - ''elements_type'': model used for each element');
disp('                  - ''elementnumbering'': numbering of elements');
disp('                  - ''vertexnumbering'': numbering of vertices');
disp('                  - ''highlightelements'': to highlight elements to highlight the element list');
disp('                  - ''highlightvertices'': to highlight vertices (use highlight option to enter the vertex list');
disp('                  - ''mesh'': draw mesh using trisurf');
disp('                  - ''referential'': stressbalance referential');
disp('                  - ''riftvel'': velocities along rifts');
disp('                  - ''riftrelvel'': relative velocities along rifts');
disp('                  - ''riftpenetration'': penetration levels for a fault');
disp('                  - ''riftfraction'': fill fractions for every node of the rifts');
disp('                  - ''rifts'': plot mesh with an offset so that rifts are visible');
disp('                  - ''strainrate_tensor'': plot the components of the strainrate tensor (exx,eyy,ezz,exy,exz,eyz) if computed');
disp('                  - ''strainrate_principal'': plot the strainrate tensor principal axis and principal values)');
disp('                  - ''strainrate_principalaxis1'': arrow plot the first principal axis of the strainrate tensor(replace 1 by 2 or 3 if needed)');
disp('                  - ''stress_tensor'': plot the components of stress tensor (sxx,syy,szz,sxy,sxz,syz) if computed');
disp('                  - ''stress_principal'': plot the stress tensor principal axis and principal values');
disp('                  - ''stress_principalaxis1'': arrow plot the first principal axis of the stress tensor(replace 1 by 2 or 3 if needed)');
disp('                  - ''transient_results'': this will display all the time steps of a transient run (use steps to specify the steps requested)');
disp('                  - ''transient_vel'': this will display the velocity for the time steps requested in ''steps'' of a transient run');
disp('                  - ''transient_vel'': vel can be by any field of the transient results (vx, vy, vz, vel, temperature, melting, pressure, bed, thickness, surface)');
disp('                  - ''transient_field'': dynamic plot of results. specify ''steps'' option, as fell as ''field'' (defaults are all steps, for ''Vel'' field)');
disp('                  - ''transient_movie'': this will display the time steps of a given field of a transient run');
disp('                  - ''transient_movie_field'': field to be displayed when doing  transient_movie data display');
disp('                  - ''transient_movie_output'': filename if output is desired for movie');
disp('                  - ''transient_movie_time'': time for each image (default 2 seconds)');
disp('                  - ''thermaltransient_results'': this will display all the time steps of a thermal transient run');
disp('                  - ''qmuhistnorm'': histogram normal distribution. needs option qmudata');
disp('                  - ''qmumean'': plot of mean distribution in sampling analysis with scaled response. needs option qmudata for descriptor');
disp('                  - ''qmustddev'': plot of stddev distribution in sampling analysis with scaled response. needs option qmudata for descriptor');
disp('                  - ''part_hist'': partitioning node and area histogram');
disp('                  - ''quiver'': quiver plot');

disp('       ''axis'': same as standard matlab option (''equal'',''off'',''equal on'',...)');
disp('       ''basin'': zoom on a given basin (''pineislandglacier'',''ronneiceshelf'', use isbasin to identify a basin');
disp('                 ''basindeltax'': in m');
disp('                 ''showbasins'': write lables for every existing basin name around the center of the plot');
disp('       ''caxis'': modify  colorbar range. (array of type [a b] where b>=a)');
disp('       ''backgroundcolor'': plot background color. (default is ''w'')');
disp('       ''figurebackgroundcolor'': figure background color. (default is ''none'')');
disp('       ''coord'':  ''xy'' (default) or ''latlon''');
disp('       ''colorlevels'':  N or {value1,valu2,value3,...} used if quiver, use different colors for the given number of colors or limits');
disp('       ''colorbar'': add colorbar (string ''on'' or ''off'')');
disp('       ''colorbartitle'': colorbar title (string)');
disp('       ''colorbarYlabel'': colorbar Y label (string)');
disp('       ''colorbarpos'': [x,y,dx,dy] where x,y,dx and dy are within [0 1]');
disp('       ''colorbarcornerposition'': ''West'',''North'',etc ...');
disp('       ''colorbartitlerotation'': -90, etc ...');
disp('       ''colorbarfontsize'': specify colorbar fontsize');
disp('       ''colorbarwidth'': multiplier (default 1) to the default width colorbar');
disp('       ''colorbarheight'': multiplier (default 1) to the default height colorbar');
disp('       ''colormap'': same as standard matlab option (''jet'',''hsv'',''cool'',''spring'',''gray'',''Ala'',''Rignot'',...)');
disp('       ''contourlevels'': N or {value1,valu2,value3,...} add the contours of the specified values or N contours');
disp('       ''contourticks'': ''on'' or ''off'' to display the ticks of the contours');
disp('       ''contouronly'': ''on'' or ''off'' to display the contours on a white background');
disp('       ''contourcolor'': ticks and contour color');
disp('       ''density'': density of quivers (one arrow every N nodes, N integer)');
disp('       ''inset'': add an inset (zoom) of the current figure if 1 (use ''insetx'', ''insety'' and ''insetpos'' to determine the inset position and content)');
disp('       ''insetx'': [min(x) max(x)] where min(x) and max(x) are values determining the inset content');
disp('       ''insety'': [min(y) max(y)] where min(y) and max(y) are values determining the inset content');
disp('       ''insetpos'': [x,y,dx,dy] where x,y,dx and dy are within [0 1]');
disp('       ''streamlines'': N (number of stream lines) or {[x1 y1],...} (coordinates of seed points) add streanlines on current figure');
disp('       ''edgecolor'': same as standard matlab option EdgeColor (color name: ''black'' or RGB array: [0.5 0.2 0.8])');
disp('       ''fontsize'': same as standard matlab option (10,14,...)');
disp('       ''fontweight'': same as standard matlab option (normal: ''n'',bold: ''b'',light: ''l'',demi: ''d'')');
disp('       ''fontcolor'': same as standard matlab option');
disp('       ''highlight'': highlights certain nodes or elements when using ''nodenumbering'' or ''elementnumbering'' or ''highlightnodes '' or ''highlightelements'' option');
disp('       ''resolution'': resolution used by section value (array of type [horizontal_resolution vertical_resolution])');
disp('                       horizontal_resolution must be in meter, and vertical_resolution a number of layers');
disp('       ''showsection'': show section used by ''sectionvalue'' (string ''on'' or a number of labels)');
disp('       ''sectionvalue'': give the value of data on a profile given by an Argus file (string ''Argusfile_name.exp'')');
disp('       ''profile'': give the value of data along a vertical profile ([xlocation ylocation])');
disp('       ''smooth'': smooth element data (string ''yes'' or integer)');
disp('       ''title'': same as standard matlab option');
disp('       ''view'': same as standard matlab option (ex: 2, 3 or [90 180]');
disp('       ''xlim'': same as standard matlab option (ex: [0 500])');
disp('       ''ylim'': same as standard matlab option');
disp('       ''zlim'': same as standard matlab option');
disp('       ''xlabel'': same as standard matlab option (ex:''km'')');
disp('       ''ylabel'': same as standard matlab option');
disp('       ''xticklabel'': specifiy xticklabel');
disp('       ''yticklabel'': specifiy yticklabel');
disp('       ''overlay'': yes or no. This will overlay a radar amplitude image behind');
disp('       ''overlay_image'': path to overlay image. provide overlay_xlim, overlay_ylim, overlay_xposting and overlay_yposting options also');
disp('       ''contrast'': (default 1) coefficient to add contrast to the radar amplitude image used in overlays');
disp('       ''highres'': resolution of overlayed radar amplitude image (default is 0, high resolution is 1).');
disp('       ''alpha'': transparency coefficient (the higher, the more transparent). Default is 1.5');
disp('       ''scaling'': scaling factor used by quiver plots. Default is 0.4');
disp('       ''autoscale'': set to ''off'' to have all the quivers with the same size. Default is ''on''');
disp('       ''expdisp'': plot exp file on top of a data plot. provide exp file as an argument (use a cell of strings if more than one)');
disp('       ''expstyle'': marker style for expdisp plot (use a cell of strings if more than one)');
disp('       ''linewidth'': line width for expdisp plot (use a cell of strings if more than one)');
disp('       ''border'': size of display border (in pixels). active only for overlay plots');
disp('       ''text'': print string, use a cell of strings if more than one');
disp('       ''textposition'': [x y] position of text, use a cell of strings if more than one');
disp('       ''textsize'':  same as standard ''FontSize'' matlab option applied to text, use a cell of strings if more than one');
disp('       ''textweight'':  same as standard ''FontWeight'' matlab option applied to text, use a cell of strings if more than one');
disp('       ''textcolor'':  same as standard ''color'' matlab option applied to text, use a cell of strings if more than one');
disp('       ''textrotation'':  same as standard ''Rotation'' matlab option applied to text, use a cell of strings if more than one');
disp('       ''mask'': list of flags of size numberofnodes or numberofelements. Only ''true'' values are plotted ');
disp('       ''nan'': value assigned to NaNs (convenient when plotting BC)');
disp('       ''partitionedges'': ''off'' by default. overlay plot of partition edges');
disp('       ''log'': value of log');
disp('       ''latlon'': ''on'' or {latstep lonstep [resolution [color]]} where latstep,longstep and resolution are in degrees, color is a [r g b] array');
disp('       ''latlonnumbering'': ''on'' or {latgap longap colornumber latangle lonangle} where latgap and longap are pixel gaps for the numbers,'); 
disp('       ''latlonclick'': ''on'' to click on latlon ticks positions');
disp('                   colornumber is a [r g b] array and latangle and lonangle are angles to flip the numbers');
disp('       ''northarrow'': add an arrow pointing north, ''on'' for default value or [x0 y0 length [ratio width fontsize]] where (x0,y0) are the coordinates of the base, ratio=headlength/length');
disp('       ''offset'': mesh offset used by ''rifts'', default is 500');
disp('       ''scaleruler'': add a scale ruler, ''on'' for default value or [x0 y0 length width numberofticks] where (x0,y0) are the coordinates of the lower left corner');
disp('       ''showregion'': show domain in Antarctica on an inset, use ''insetpos'' properties');
disp('       ''visible'': ''off'' to make figure unvisible, default is ''on''');
disp('       ''wrapping'': repeat ''n'' times the colormap (''n'' must be an integer)');
disp('       ''unit'': by default, in m, otherwise, ''km'' is available');
disp('       ''legend_position'': by default, ''NorthEasth''');
disp('       ''qmudata'': data for qmu  plots.');
disp('                  {dresp1   ,dresp2  ,hmin,hmax,hnint} or {samp,desc,mu,sigma,hmin,hmax,hnint}');
disp('                  where dresp1 is a structure array of responses (where we need samp and desc), ');
disp('                  dresp2 is a structure array of responses (where we only need mu and sigma)');
disp('                  hmin,hmax and hnint are the minimum, maximum and number of intervals of the histogram (optional)');
disp('       ''figposition'': position of figure: ''fullscreen'', ''halfright'', ''halfleft'', ''portrait'', ''landscape'',... (hardcoded in applyoptions.m)');
disp('       ''offsetaxispos'': offset of current axis position to get more space (ex: [-0.02 0  0.04 0])');
disp('       ''axispos'': axis position to get more space');
disp('       ''hmin'': (numeric, minimum for histogram)');
disp('       ''hmax'': (numeric, maximum for histogram)');
disp('       ''hnint'': (numeric, number of intervals for histogram)');
disp('       ''ymin1'': (numeric, minimum of histogram y-axis)');
disp('       ''ymax1'': (numeric, maximum of histogram y-axis)');
disp('       ''ymin2'': (numeric, minimum of cdf y-axis)');
disp('       ''ymax2'': (numeric, maximum of cdf y-axis)');
disp('       ''cdfplt'': (char, ''off'' to turn off cdf line plots)');
disp('       ''cdfleg'': (char, ''off'' to turn off cdf legends)');
disp('       ''segmentnumbering'': (''off'' by default)');
disp('       ''kmlgroundoverlay'': (''off'' by default)');
disp('       ''kmlfilename'': (''tempfile.kml'' by default)');
disp('       ''kmlroot'': (''./'' by default)');
disp('       ''kmlimagename'': (''tempimage'' by default)');
disp('       ''kmlimagetype'': (''png'' by default)');
disp('       ''kmlresolution'': (1 by default)');
disp('       ''kmlfolder'': (''Ground Overlay'' by default)');
disp('       ''kmlfolderdescription'': ('''' by default)');
disp('       ''kmlgroundoverlayname'': ('''' by default)');
disp('       ''kmlgroundoverlaydescription'': ('''' by default)');

disp('       any options (except ''data'') can be followed by ''#i'' where ''i'' is the subplot number, or ''#all'' if applied to all plots');
disp('  ');
disp('   Examples:');
disp('       plotmodel(md,''data'',''vel'',''data'',''mesh'',''view#2'',3,''colorbar#all'',''on'',''axis#1'',''off equal'')');
disp('       plotmodel(md,''data'',''highlightelements'',''highlight'',[1 4 10],''expdisp'',{''domain1.exp'' ''domain2.exp'' ''domain3.exp''})');
