function ev= extractVisitData(evaldata, k)
evaldata(evaldata.visitnum~=k,:) = [];
evaldata.visitnum= []; %removing redundant column
ev= evaldata;