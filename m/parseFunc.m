function [parsedFunc] = parseFunc(func,pre_param)

if(nargin == 1)
    pre_param = 'config.nf.parameters.';
end

funcCall = func2str(func);
[beginIdx,endIdx] = regexp(funcCall,'(\(.*?)\)');

if(endIdx(1) == beginIdx(1) + 1)
    parsedFunc = '()';
    return;
end

callExpr = funcCall((beginIdx(1)+1):(endIdx(1)-1));
splitStr = regexp(callExpr,',','split');

parsedFunc = '(';

for strNum = 1:1:numel(splitStr)
    parsedFunc = strcat(parsedFunc , pre_param , splitStr{strNum} , ',' );
end

parsedFunc(end) = ')';

end

