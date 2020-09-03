function Data = myNormalization( inputData, RorC )
%MYNORMLIZATION Summary of this function goes here
%   Detailed explanation goes here
% this is a function for normalizing the data matrix that make inputData
% have unit length on its row or column
% inputData:  the input Data metrix
%RorC:  'row':: force its row to have unit length
%           'column':: force its column to have unit length
if isempty(RorC)
    RorC='column'; % default: force its column to have unit length
end

if isempty(inputData)
    disp('no Data!')
    return;
end

[H,W]=size(inputData);
if strcmp(RorC,'row')
    Data=inputData./repmat(sqrt(sum(inputData.^2,2)),[1,W]);
elseif strcmp(RorC,'column')    
    Data=inputData./repmat(sqrt(sum(inputData.^2)),[H,1]);
end

end

