%
%  definition for the response_function class.
%
%  [rf]=response_function(varargin)
%
%  where the required varargin are:
%    descriptor    (char, description, '')
%  and the optional varargin and defaults are:
%    respl         (double vector, response levels, [])
%    probl         (double vector, probability levels, [])
%    rell          (double vector, reliability levels, [])
%    grell         (double vector, gen. reliability levels, [])
%
%  note that zero arguments constructs a default instance; one
%  argument of the class copies the instance; and one or more
%  arguments constructs a new instance from the arguments.
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
classdef response_function
    properties
        descriptor='';
        respl     =[];
        probl     =[];
        rell      =[];
        grell     =[];
    end

    methods
        function [rf]=response_function(varargin)

            switch nargin

%  create a default object

                case 0

%  copy the object or create the object from the input

                otherwise
                    if  (nargin == 1) && isa(varargin{1},'response_function')
                        rf=varargin{1};
                    else
                        asizec=num2cell(array_size(varargin{1:min(nargin,1)}));
                        rf(asizec{:})=response_function;
                        clear asizec

                        if ischar(varargin{1})
                            varargin{1}=cellstr(varargin{1});
                        end
                        for i=1:numel(rf)
                            if (numel(varargin{1}) > 1)
                                rf(i).descriptor=varargin{1}{i};
                            else
                                rf(i).descriptor=[char(varargin{1}) string_dim(rf,i,'vector')];
                            end
                        end

                        if (nargin >= 2)
                            for i=1:numel(rf)
                                rf(i).respl     =varargin{2};
                            end
                            if (nargin >= 3)
                                for i=1:numel(rf)
                                    rf(i).probl     =varargin{3};
                                end
                                if (nargin >= 4)
                                    for i=1:numel(rf)
                                        rf(i).rell      =varargin{4};
                                    end
                                    if (nargin >= 5)
                                        for i=1:numel(rf)
                                            rf(i).grell     =varargin{5};
                                        end

                                        if (nargin > 5)
                                            warning('response_function:extra_arg',...
                                                'Extra arguments for object of class ''%s''.',...
                                                class(rf));
                                        end
                                    end
                                end
                            end
                        end
                    end
            end

        end

        function []=disp(rf)

        %  display the object

            disp(sprintf('\n'));
            for i=1:numel(rf)
                disp(sprintf('class ''%s'' object ''%s%s'' = \n',...
                    class(rf),inputname(1),string_dim(rf,i)));
                disp(sprintf('    descriptor: ''%s'''  ,rf(i).descriptor));
                disp(sprintf('         respl: %s'      ,string_vec(rf(i).respl)));
                disp(sprintf('         probl: %s'      ,string_vec(rf(i).probl)));
                disp(sprintf('          rell: %s'      ,string_vec(rf(i).rell)));
                disp(sprintf('         grell: %s\n'    ,string_vec(rf(i).grell)));
            end

        end

        function [desc]  =prop_desc(rf,dstr)
            desc=cell(1,numel(rf));
            for i=1:numel(rf)
                if ~isempty(rf(i).descriptor)
                    desc(i)=cellstr(rf(i).descriptor);
                elseif ~isempty(inputname(1))
                    desc(i)=cellstr([inputname(1) string_dim(rf,i,'vector')]);
                elseif exist('dstr','var')
                    desc(i)=cellstr([dstr         string_dim(rf,i,'vector')]);
                else
                    desc(i)=cellstr(['rf'         string_dim(rf,i,'vector')]);
                end
            end
            desc=allempty(desc);
        end
        function [stype] =prop_stype(rf)
            stype={};
        end
        function [scale] =prop_scale(rf)
            scale=[];
        end
        function [weight]=prop_weight(rf)
            weight=[];
        end
        function [lower] =prop_lower(rf)
            lower=[];
        end
        function [upper] =prop_upper(rf)
            upper=[];
        end
        function [target]=prop_target(rf)
            target=[];
        end
        function [respl,probl,rell,grell]=prop_levels(rf)
            respl=cell(1,numel(rf));
            probl=cell(1,numel(rf));
            rell =cell(1,numel(rf));
            grell=cell(1,numel(rf));
            for i=1:numel(rf)
                respl(i)={rf(i).respl};
                probl(i)={rf(i).probl};
                rell (i)={rf(i).rell};
                grell(i)={rf(i).grell};
            end
            respl=allempty(respl);
            probl=allempty(probl);
            rell =allempty(rell);
            grell=allempty(grell);
        end
    end

    methods (Static)
        function [rdesc]=dakota_write(fidi,dresp,rdesc)

%  collect only the responses of the appropriate class

            rf=struc_class(dresp,'response_function');

%  write responses

            [rdesc]=rlist_write(fidi,'response_functions','response_function',rf,rdesc);
        end

        function []=dakota_rlev_write(fidi,dresp,params)

%  collect only the responses of the appropriate class

            rf=struc_class(dresp,'response_function');

%  write response levels

            rlev_write(fidi,rf,params);
        end
    end
end
