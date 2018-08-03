%
%  function to return an empty double array if all array elements are NaN
%
%  function [dout]=allnan(din)
%
function [dout]=allnan(din)

if ~nargin
    help allnan
    return
end

for i=1:numel(din)
    if ~isnan(din(i))
        dout=din;
        return
    end
end
dout=[];

end

