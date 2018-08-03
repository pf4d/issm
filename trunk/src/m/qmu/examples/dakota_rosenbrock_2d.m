%
%  multi-dimensional parameter study for rosenbrock case
%  (see Users4.2.pdf, Sec. 2.4.1)
%
%  "Copyright 2009, by the California Institute of Technology.
%  ALL RIGHTS RESERVED. United States Government Sponsorship
%  acknowledged. Any commercial use must be negotiated with
%  the Office of Technology Transfer at the California Institute
%  of Technology.  (J. Schiermeier, NTR 47078)
%
%  This software may be subject to U.S. export control laws.
%  By accepting this  software, the user agrees to comply with
%  all applicable U.S. export laws and regulations. User has the
%  responsibility to obtain export licenses, or other export
%  authority as may be required before exporting such information
%  to foreign countries or providing access to foreign persons."
%

function [dout,ddat]=dakota_rosenbrock_2d()

%  define dakota variables as continuous design
%  (may use set 1 or set 2 of variables, but not both)
dvar(1).cdv(1)=continuous_design('x1',0,-2,2);
dvar(1).cdv(2)=continuous_design('x2',0,-2,2);
dvar(2).x1=continuous_design('',0,-2,2);
dvar(2).x2=continuous_design('',0,-2,2);

%  define dakota response as objective function
%  (may use set 1 or set 2 of responses, but not both)
dresp(1).of=objective_function('f');
dresp(2).f=objective_function('');

%  define dakota method and specify method-dependent parameters
dmeth=dakota_method('multidim');
dmeth=dmeth_params_set(dmeth,'partitions',[8 8]);

%  specify method-independent parameters
%  (dakota_in_params does not need to be called, but provides a template)
dparams=dakota_in_params([]);
dparams.direct=true;
dparams.analysis_driver='rosenbrock';
dparams.tabular_graphics_data=true;
dparams.tabular_graphics_file='dakota_rosenbrock_2d.dat';

%  write out dakota input file
dakota_in_write(dmeth,dvar(2),dresp(2),dparams,'dakota_rosenbrock_2d.in')

%  execute dakota
!dakota -i dakota_rosenbrock_2d.in -o dakota_rosenbrock_2d.out

%  read dakota output and tabular data files
%  (output file for parameter studies has no interesting info)
[method,dout]=dakota_out_parse('dakota_rosenbrock_2d.out');
[~     ,ddat]=dakota_out_parse('dakota_rosenbrock_2d.dat');

%  perform any desired plotting
plot_rvsv_surf(ddat,{'x1','x2'},ddat,{'f'})

end

