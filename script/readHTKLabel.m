 function labelStruct = readHTKLabel(fname)

%%  get file id

fid = fopen(fname,'r');
fgets(fid);
i=1;
s=fgets(fid);
labelCell=cell(1,3);
while(s>0)
    
    c = textscan(s,'%n%n%s');
    if(strcmp(c{3},'sp') || strcmp(c{3},'sil'))
        s=fgets(fid);
        continue;
    end
    labelCell{1}(i,1)=c{1};
    labelCell{2}(i,1)=c{2};
    labelCell{3}(i,1)=c{3};
    s=fgets(fid);
    i=i+1;
end
fclose(fid);

startTime = labelCell{1}(1:end)/10000000;
endTime = labelCell{2}(1:end)/10000000;
phonemeLabel = labelCell{3}(1:end);

labelStruct.segment = [startTime endTime];
labelStruct.phoneme = char(phonemeLabel);
