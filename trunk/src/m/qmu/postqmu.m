function md=postqmu(md)
%INPUT function md=postqmu(md,qmufile,qmudir)
%Deal with dakota output results in files.

%  check to see if dakota returned errors in the err file
qmuerrfile=[md.miscellaneous.name '.qmu.err'];

if exist(qmuerrfile,'file')
   fide=fopen(qmuerrfile,'r');
   fline=fgetl(fide);
   if ischar(fline)
       while ischar(fline)
           disp(sprintf('%s',fline));
           fline=fgetl(fide);
       end
       status=fclose(fide);
       cd ../
       error(['Dakota returned error in ''' qmuerrfile ' file.  ''' qmudir ''' directory retained.'])
    end
    status=fclose(fide);
end

%parse inputs and results from dakota
qmuinfile=[md.miscellaneous.name '.qmu.in'];
qmuoutfile=[md.miscellaneous.name '.qmu.out'];

%[method,dvar,dresp_in]=dakota_in_parse(qmuinfile);
%dakotaresults.method   =method;
%dakotaresults.dvar     =dvar;
%dakotaresults.dresp_in =dresp_in;

[method,dresp_out,scm,pcm,srcm,prcm]=dakota_out_parse(qmuoutfile);
dakotaresults.dresp_out=dresp_out;
dakotaresults.scm      =scm;
dakotaresults.pcm      =pcm;
dakotaresults.srcm     =srcm;
dakotaresults.prcm     =prcm;

if exist('dakota_tabular.dat','file')
    [method,dresp_dat                  ]=dakota_out_parse('dakota_tabular.dat');
    dakotaresults.dresp_dat=dresp_dat;
end

%put dakotaresults in their right location.
md.results.dakota=dakotaresults;

%  move all the individual function evalutations into zip files
if ~md.qmu.isdakota,
	system('zip -mq params.in.zip params.in.[1-9]*');
	system('zip -mq results.out.zip results.out.[1-9]*');
	system('zip -mq matlab.out.zip matlab*.out.[1-9]*');
end
