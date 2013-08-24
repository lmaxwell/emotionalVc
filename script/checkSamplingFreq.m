emotion=strvcat('neutral' ,'happy', 'angry','sad')
pdir=[pwd '/../'];
plotFigures = 0; % If you do not need plots, set this variable 0.
for jemotion=1:4%size(emotion,1)
 
    
dir=[pdir  deblank(emotion(jemotion,:)) '/wav/'];

    for iFile=201:500
        %fprintf('read %d\n',iFile);
        [x,fs]=wavread([dir num2str(ifile) '.wav']);
        
        if (fs ~=16000)
            fprintf('fs of %s%d.wav is not 16000\n',[deblank(emotion(jemotion,:)) '/wav/'],iFile);
        end
    end
end