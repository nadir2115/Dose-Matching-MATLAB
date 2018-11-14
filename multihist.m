% Create multiple histograms for all variables in a matrix
function multihist(data)

m= size(data,1);
n= size(data,2);
figure

if n>20 %labels for one-hot-encoding
labels= ["EQ INDEX","SIS hand","ufugm","ufugmcat","RNLIadj","NIHtot",...
    "wmftmean MAPA","wmftmean LAPA","logmean MAPA"...
    "log meanLA PA","FAS score PA","severity","FCtotalOT","DEMOedu"...
    "CSgender","Cssidehemi","onset to rand","concordance","age@rand"...
    "DEMOsmoke","site code 1","site code 2","site code 3","	Hispanic",...
    "Native/Asian/PI","African-Am.","Caucasian","Other Race",...
    "w/o hemorrhagic", "StrokeL RHem","StrokeL LHem","StrokeL Bstem",...
    "StrokeL Other","Dose hours"];

    for i= 1:n
        subplot(ceil(n/9),9,i) 
        hist(data(:,i))
        title(labels(i))
    end
    
else    %labels for just continuous variables
    labels= ["EQ INDEX","SIS hand","ufugm","RNLIadj","NIHtot",...
        "wmftmean MAPA","wmftmean LAPA","logmean MAPA","log meanLA PA",...
        "FAS score PA", "onset to rand","age@rand", "Dose hours" ];
     for i= 1:n
        subplot(ceil(n/7),7,i) 
        hist(data(:,i))
        title(labels(i))
    end
end


