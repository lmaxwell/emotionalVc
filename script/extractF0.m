emotion=strvcat('neutral' ,'happy', 'angry','sad')
pdir=[pwd '/../'];
plotFigures = 0; % If you do not need plots, set this variable 0.
for jemotion=4:4%size(emotion,1)
 
    
dir=[pdir  deblank(emotion(jemotion,:))];
mkdir([ dir '/logf0']);

    for ifile=201:500
        fname=[dir '/wav/' num2str(ifile) '.wav'];

        [x,fs] = wavread(fname);

        x = x(:,1); %   Make sure it is a comum vector.
%         soundsc(x,fs) % Playback sound

        %  Extract source information
        opt.f0ceil=500;
        opt.f0floor=100;
        r = exF0candidatesTSTRAIGHTGB(x,fs,opt); % Extract F0 information



        rc = autoF0Tracking(r,x); % Clean F0 trajectory by tracking
        rc.vuv = refineVoicingDecision(x,rc);

        if plotFigures
            figure;
            plot(rc.temporalPositions,rc.f0);grid on
            set(gca,'fontsize',14);
            xlabel('time (s)')
            ylabel('fundamental frequency (Hz)');
            title('fundamental frequency')
        end;

        q = aperiodicityRatioSigmoid(x,rc,1,2,0) % aperiodicity extraction
        
        f0file=fopen([dir '/logf0/' num2str(ifile) '.f0'],'w');
        logf0=log(q.f0).*q.vuv;
          findvoiced=sign( diff([0; logf0(:)>0; 0]) );
    voicedind=[find(findvoiced>0)'; find(findvoiced<0)'-1];
   for k=1:size(voicedind,2)-1
    logf0(voicedind(2,k):voicedind(1,k+1))=linspace(logf0(voicedind(2,k)),logf0(voicedind(1,k+1)),voicedind(1,k+1)-voicedind(2,k)+1);
   end
   
   
    first=find(logf0,1,'first');
    last=find(logf0,1,'last');
    logf0(1:first-1)=logf0(first);
    logf0(last+1:end)=logf0(last);   
     
    
    
    for ilen=3:length(logf0)-2
        logf0(ilen)=median(logf0(ilen-2:ilen+2));
    end

        fprintf(f0file,'%f\n',logf0);
        fclose(f0file);
        
        if plotFigures;
            displayAperiodicityStructure(q,1);
        end;
    
    end


end
% 
% %
% emotion=strvcat('neutral' ,'happy', 'angry','sad')
% pdir=[pwd '/../'];
% plotFigures = 0; % If you do not need plots, set this variable 0.
% for jemotion=3:4%size(emotion,1)
%  
%     
% dir=[pdir  deblank(emotion(jemotion,:))];
% mkdir([ dir '/logf0']);
% 
%     for ifile=201:500
%         fname=[dir '/wav/' num2str(ifile) '.wav'];
% 
%        
%         
%         f0=load([dir '/logf0/' num2str(ifile) '.f0']);
%         logf0=zeros(length(f0),1);
%         for i=1:length(f0)
%             if(f0(i)>0)
%                 
%                 logf0(i)=log(f0(i));
%             end
%         end
%         f0file=fopen([dir '/logf0/' num2str(ifile) '.f0'],'w');
%         
%         
%         
%         fprintf(f0file,'%f\n',logf0);
%         fclose(f0file);
%         
%         if plotFigures;
%             displayAperiodicityStructure(q,1);
%         end;
%     
%     end
% 
% 
% end
