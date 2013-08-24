  function f_morph(iFile,emotion,method,durRate,f0Rate,spRate,workDir)
%F_MORPH Summary of this function goes here
%   Detailed explanation goes here
addpath('/home/lixian/speech/MatlabCode/STRAIGHT/TandemSTRAIGHTmonolithicPackage002Test');
prjdir=[pwd '/..'];
%  iFile=201;
% emotion='sad';
% method='twotier';
% durRate=1;
% f0Rate=0;
% spRate=0;
timeRate=durRate;
fprintf('generating %s/%d_%d%d%d.wav\n',workDir,iFile,durRate,f0Rate,spRate);
%% read A
path=[prjdir '/neutral/wav/'];
file=[num2str(iFile) '.wav'];
SA=getSTRAIGHTobject(path,file,0);

STRAIGHTobject=SA;
mSubstrate = morphingSubstrateNewAP;
mSubstrate.dataDirectoryForSpeakerA = STRAIGHTobject.dataDirectory;
mSubstrate.fileNameForSpeakerA= STRAIGHTobject.dataFileName;
mSubstrate.waveformForSpeakerA= STRAIGHTobject.waveform;
mSubstrate.samplintFrequency = STRAIGHTobject.samplingFrequency;
mSubstrate.f0TimeBaseOfSpeakerA = STRAIGHTobject.refinedF0Structure.temporalPositions;
%mSubstrate.f0OfSpeakerA = STRAIGHTobject.refinedF0Structure.f0CandidatesMap(1,:)';
mSubstrate.f0OfSpeakerA = STRAIGHTobject.refinedF0Structure.f0;
mSubstrate.spectrogramTimeBaseOfSpeakerA = ...
    STRAIGHTobject.SpectrumStructure.temporalPositions;
mSubstrate.STRAIGHTspectrogramOfSpeakerA = ...
    STRAIGHTobject.SpectrumStructure.spectrogramSTRAIGHT;
mSubstrate.aperiodicityTimeBaseOfSpeakerA = ...
    STRAIGHTobject.AperiodicityStructure.temporalPositions;
mSubstrate.aperiodicityOfSpeakerA = ...
    STRAIGHTobject.AperiodicityStructure;
path=[prjdir '/neutral/lab/phone/'];
file=[num2str(iFile) '.lab'];
labelA = readHTKLabel([path file]);
mSubstrate = morphingSubstrateNewAP(mSubstrate,'set',...
    'temporaAnchorOfSpeakerA',labelA.segment(:,1));
frequencyAnchors = setBlankFrequencyAnchors(labelA);
mSubstrate.frequencyAnchorOfSpeakerA = frequencyAnchors;


%% read B 
path=[prjdir '/' emotion '/wav/'];
file=[num2str(iFile) '.wav'];
a=1;
 SB=getSTRAIGHTobject(path,file,0);
STRAIGHTobject=SB;
mSubstrate.dataDirectoryForSpeakerB = STRAIGHTobject.dataDirectory;
mSubstrate.fileNameForSpeakerB = STRAIGHTobject.dataFileName;
mSubstrate.waveformForSpeakerB= STRAIGHTobject.waveform;
mSubstrate.samplintFrequency = STRAIGHTobject.samplingFrequency;
mSubstrate.f0TimeBaseOfSpeakerB = STRAIGHTobject.refinedF0Structure.temporalPositions;
%mSubstrate.f0OfSpeakerB = STRAIGHTobject.refinedF0Structure.f0CandidatesMap(1,:)';
mSubstrate.f0OfSpeakerB = STRAIGHTobject.refinedF0Structure.f0;
mSubstrate.spectrogramTimeBaseOfSpeakerB = ...
    STRAIGHTobject.SpectrumStructure.temporalPositions;
mSubstrate.STRAIGHTspectrogramOfSpeakerB = ...
    STRAIGHTobject.SpectrumStructure.spectrogramSTRAIGHT;
mSubstrate.aperiodicityTimeBaseOfSpeakerB = ...
    STRAIGHTobject.AperiodicityStructure.temporalPositions;
mSubstrate.aperiodicityOfSpeakerB = ...
    STRAIGHTobject.AperiodicityStructure;
path=[prjdir '/' emotion '/lab/phone/'];
file=[num2str(iFile) '.lab'];
labelB=readHTKLabel([path file]);
mSubstrate = morphingSubstrateNewAP(mSubstrate,'set',...
    'temporaAnchorOfSpeakerB',labelB.segment(:,1));
frequencyAnchors = setBlankFrequencyAnchors(labelB);
mSubstrate.frequencyAnchorOfSpeakerB = frequencyAnchors;


%% morph
if(f0Rate==0)
 temp=load([workDir '/' num2str(iFile) '.phrase.f0'])+load([workDir '/' num2str(iFile) '.tone.f0']);
% temp=load([directoryBase '/logf0/' num2str(iFile) '.f0' ]);
logf0=[temp ;0*ones(length(mSubstrate.f0TimeBaseOfSpeakerA)-length(temp),1)];

   
   
    first=find(logf0,1,'first');
    last=find(logf0,1,'last');
    logf0(1:first-1)=logf0(first);
    logf0(last+1:end)=logf0(last);   
     
    
    
   

mSubstrate.f0OfSpeakerA=exp(logf0);
end
% if(strcmp(method,'linear'))
%  logf0=load([workDir '/' num2str(iFile) '.f0']);
% 
% mSubstrate.f0OfSpeakerA=exp(logf0);
% end
if(spRate==0)
    mSubstrate.STRAIGHTspectrogramOfSpeakerA=load([prjdir '/vc/to' emotion '/test/wav/neutral-' emotion '/' num2str(iFile) '.sp' ])';
end
mRateCommon = 1;
mSubstrate = morphingSubstrateNewAP(mSubstrate,'generate','morphingTimeAxis',mRateCommon);
%---- set default temporal morphing rate
% fs = mSubstrate.samplintFrequency;
% mRate.time =ones(length(mSubstrate.morphingTimeAxies),1)
% mRate.F0 = mRate.time;
% mRate.frequency = mRate.time;
% mRate.spectrum = mRate.time;
% mRate.aperiodicity = mRate.time;
% mSubstrate = morphingSubstrateNewAP(mSubstrate,'set','temporalMorphingRate',mRate);
mRate.time =timeRate+0*(mSubstrate.morphingTimeAxis(:))/mSubstrate.morphingTimeAxis(end);
zeroRate=0*(mSubstrate.morphingTimeAxis(:))/mSubstrate.morphingTimeAxis(end);
mRate.F0 = zeroRate+f0Rate;


mRate.frequency = zeroRate;
mRate.spectrum = zeroRate+spRate;
mRate.aperiodicity = zeroRate;
mSubstrate = morphingSubstrateNewAP(mSubstrate,'set','temporalMorphingRate',mRate);

morphedSignal = generateMorphedSpeechNewAP(mSubstrate);

wavwrite(morphedSignal.outputBuffer/max(abs(morphedSignal.outputBuffer))*0.99,...
    mSubstrate.samplintFrequency,[workDir '/' num2str(iFile) '_' num2str(durRate) num2str(f0Rate) num2str(spRate) '.wav']);
sound(morphedSignal.outputBuffer/max(abs(morphedSignal.outputBuffer))*0.99,...
    mSubstrate.samplintFrequency);

% 
% referenceTime = mSubstrate.morphingTimeAxis;
% mRateFields = {'spectrum' 'frequency' 'aperiodicity' 'F0' 'time'};


% knobNames = mSubstrate.knobNames;
% xdataStr = 'xdata';
% ydataStr = 'ydata';
% referenceTime = mSubstrate.morphingTimeAxis;
% for ii = 1:5
%     eval(['xdata = get(mSubstrate.' knobNames{ii} '.plot,xdataStr);']);
%     eval(['ydata = get(mSubstrate.' knobNames{ii} '.plot,ydataStr);']);
%     timeFunction = interp1(xdata,ydata,referenceTime);
%     eval(['mRate.' mRateFields{ii} '=timeFunction;']);
%     knobYdata{ii} = ydata;
% end;
% mSubstrate = morphingSubstrateNewAP(mSubstrate,'set','temporalMorphingRate',mRate);
% dBSTRAIGHTspectrogram = 10*log10(mSubstrate.STRAIGHTspectrogramOfSpeakerB);
% referenceTimeAxis = mSubstrate.spectrogramTimeBaseOfSpeakerB;
% timeAnchors = mSubstrate.temporaAnchorOfSpeakerB;
% 
% referenceTimeDefiningValue = 1;
%   
% mSubstrate = morphingSubstrateNewAP(mSubstrate,...
%     'generate','morphingTimeAxis',referenceTimeDefiningValue);
% mSubstrate.referenceTimeAxis = referenceTimeAxis;


%%
% path=[pwd '/neutral/wav/'];
% file='201.wav';
% STRAIGHTobject.creationDate = datestr(now,30);
% STRAIGHTobject.standAlone = false;
% [x,fs] = wavread([path file]);
% x = x(:,1)+randn(length(x(:,1)),1)*std(x(:,1))/100000; % safeguard 7/Msy/2012
% STRAIGHTobject.creationDate = datestr(now,30);
% STRAIGHTobject.dataDirectory = path;
% STRAIGHTobject.dataFileName = file;
% STRAIGHTobject.samplingFrequency = fs;
% STRAIGHTobject.waveform = x(:,1);
% STRAIGHTobject.standAlone = true;
% 
% STRAIGHTobject.F0extractionDate = datestr(now,30);
% if 1 == 1
% 
% f0Structure = ...
%     exF0candidatesTSTRAIGHTGB(STRAIGHTobject.waveform, ...
%     STRAIGHTobject.samplingFrequency);
% end;
% 
% 
% fs = STRAIGHTobject.samplingFrequency;
% STRAIGHTobject.f0 = f0Structure.f0;
% STRAIGHTobject.periodicityLevel = f0Structure.periodicityLevel;
% STRAIGHTobject.temporalPositions = f0Structure.temporalPositions;
% STRAIGHTobject.f0CandidatesMap = f0Structure.f0CandidatesMap;
% STRAIGHTobject.f0CandidatesScoreMap = f0Structure.f0CandidatesScoreMap;
% STRAIGHTobject.f0Structure = f0Structure;
% STRAIGHTobject.f0candidatesPowerMap = f0Structure.f0candidatesPowerMap;
% % 
% % f0Structure=autoF0Tracking(STRAIGHTobject.F0Structure,x);
% % f0Structure.vuv = refineVoicingDecision(x,f0Structure);
% 
% %%
% GUIuserData=STRAIGHTobject;
% refinedData = autoF0Tracking(GUIuserData,x);
% refinedData.vuv = refineVoicingDecision(x,refinedData);
% GUIuserData.f0 = refinedData.f0;
% GUIuserData.vuv = refinedData.vuv;
% GUIuserData.periodicityLevel = refinedData.periodicityLevel;
% 
% STRAIGHTobject.originalF0Structure = GUIuserData.f0Structure;
%     STRAIGHTobject.refinedF0Structure = GUIuserData.f0Structure;
%     STRAIGHTobject.refinedF0Structure.f0 = GUIuserData.f0;
%     STRAIGHTobject.waveform = GUIuserData.waveform;
%     if isfield(GUIuserData,'vuv')
%         STRAIGHTobject.refinedF0Structure.vuv = GUIuserData.vuv;
%     end;
%     STRAIGHTobject.refinedF0Structure.periodicityLevel = ...
%         GUIuserData.periodicityLevel;
% 
% 
% 
% 
% %%
% STRAIGHTobject.AperiodicityExtractionDate = datestr(now,30);
% 
% STRAIGHTobject.AperiodicityStructure = ...
%     aperiodicityRatioSigmoid(STRAIGHTobject.waveform, ...
%     STRAIGHTobject.refinedF0Structure,1,2,0); % originally 2
% 
% if isfield(STRAIGHTobject.refinedF0Structure,'vuv')
%     STRAIGHTobject.AperiodicityStructure.vuv = ...
%         STRAIGHTobject.refinedF0Structure.vuv;
%     disp('i am here');
% end
% 
% 
% STRAIGHTobject.SpectrumExtractionDate = datestr(now,30);
% %prmIn.compensationCoefficient = [-0.2818;0.0397]; % this is experimental 15/Oct./2009
% %prmIn.compensationCoefficient = -0.3804;
% %prmIn.exponentControl = 1/10;
% %prmIn.correctionForBlackman = 2.5;
% %prmIn.debugperiodicityShaping = 2.5; % original 2.4
% %prmIn.defaultF0 = 500;
% SpectrumStructure = ...
%     exSpectrumTSTRAIGHTGB(STRAIGHTobject.waveform, ...
%     STRAIGHTobject.samplingFrequency, ...
%     STRAIGHTobject.refinedF0Structure);%,prmIn);
% %STRAIGHTobject.SpectrumStructure = ...
% %    fixWindowEffects(SpectrumStructure,STRAIGHTobject.refinedF0Structure, ...
% %    prmIn.correctionForBlackman,2,1/2); % 18/Oct./2009
% STRAIGHTobject.SpectrumStructure = SpectrumStructure;
% processedSpectrum = unvoicedProcessing(STRAIGHTobject);
% STRAIGHTobject.SpectrumStructure.spectrogramSTRAIGHT = processedSpectrum;
% 
% 
% if save==1
%     STRAIGHTobject.lastUpdate = datestr(now);
%     outFileName = ['StrObj' datestr(now,30) '.mat'];
%     [file,path] = uiputfile(outFileName,'Save the STRAIGHT object');
%     if length(file) == 1 && length(path) == 1
%         if file == 0 || path == 0
%             %okInd = 0;
%             disp('Save is cancelled!');
%             return;
%         end;
%     end;
% 
%     %pathReg = regexprep(path,'\s','\\ ');
%     %eval(['save ' pathReg file ' STRAIGHTobject']);
%     save('test.mat','STRAIGHTobject');
% end


% end

