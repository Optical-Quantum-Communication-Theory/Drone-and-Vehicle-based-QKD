function channelOutput = convertData(data, dataOrder, targetOrder)
%Convert data from Ohio State. Ohio clicks: [V H L R]. Our clicks : [H V R L]
%Ex map [1 0 0 0] to [0 1 0 0]
%Import data as numeric array

dataPattern = ["V", "H", "L", "R"];
targetPattern = ["H", "V", "R", "L"];

dataH = dataPattern == "H";
dataV = dataPattern == "V"; 
dataR = dataPattern == "R"; 
dataL = dataPattern == "L"; 

for i = 1:numel(dataOrder)
    dataOrder{i} = mod(floor(dataOrder{i} ./ 10 .^ (3:-1:0)), 10);
    %convert to vector format
end
formatH = [1,0,0,0];
formatV = [0,1,0,0];
formatR = [0,0,1,0];
formatL = [0,0,0,1]; 

permCols = zeros(1, numel(targetOrder));
for patInd = 1:numel(targetOrder)
    %for every target pattern, find column in data associated with data
    %pattern
    iTargetPattern = targetOrder{patInd}; %some click pattern, target
    iDataPattern = iTargetPattern([2, 1, 4, 3]); %permute to ohio format
    idx = find(cellfun(@(cell) ismember(iDataPattern, cell, 'rows'),dataOrder)); %find ohio pattern
    %fprintf('%s\n', idx);
    permCols(patInd) = idx;  %store column number
end
permRows = [3, 1, 2, 6, 4, 5, 9, 7, 8, 12, 10, 11];
channelOutput = data(permRows, permCols); 


end