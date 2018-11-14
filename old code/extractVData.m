function ev= extractVData(evaldata, k)
evaldata(evaldata(:,3)~=k,:) = [];
evaldata(:,3)= []; %removing redundant column
ev= evaldata;