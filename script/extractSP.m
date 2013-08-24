addpath('/home/lixian/speech/MatlabCode/STRAIGHT/TandemSTRAIGHTmonolithicPackage002Test');

emotion=strvcat('neutral' ,'happy', 'angry','sad')
pdir=[pwd '/../'];
plotFigures = 0; % If you do not need plots, set this variable 0.
for jemotion=1:4%size(emotion,1)
 
    
dir=[pdir  deblank(emotion(jemotion,:))];
mkdir([ dir '/spectrum']);

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

        q = aperiodicityRatioSigmoid(x,rc,1,2,0); % aperiodicity extraction
        
        spfile=fopen([dir '/spectrum/' num2str(ifile) '.sp'],'w');
       f = exSpectrumTSTRAIGHTGB(x,fs,q);

STRAIGHTobject.waveform = x;
STRAIGHTobject.samplingFrequency = fs;
STRAIGHTobject.refinedF0Structure.temporalPositions = r.temporalPositions;
STRAIGHTobject.SpectrumStructure.spectrogramSTRAIGHT = f.spectrogramSTRAIGHT;
STRAIGHTobject.refinedF0Structure.vuv = rc.vuv;
f.spectrogramSTRAIGHT = unvoicedProcessing(STRAIGHTobject);

% sgramSTRAIGHT = 10*log10(f.spectrogramSTRAIGHT);
        fprintf(spfile,'%d ',f.spectrogramSTRAIGHT);
        fclose(spfile);
        
        if plotFigures;
            displayAperiodicityStructure(q,1);
        end;
    
    end


end
