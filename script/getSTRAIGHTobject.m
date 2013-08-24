function [STRAIGHTobject ] = getSTRAIGHTobject( path,file,save)
%GETSTRAIGHTOBJECT Summary of this function goes here
%   Detailed explanation goes here
addpath('/home/lixian/speech/MatlabCode/STRAIGHT/TandemSTRAIGHTmonolithicPackage002Test');
STRAIGHTobject.creationDate = datestr(now,30);
STRAIGHTobject.standAlone = false;
[x,fs] = wavread([path file]);
x = x(:,1)+randn(length(x(:,1)),1)*std(x(:,1))/100000; % safeguard 7/Msy/2012
STRAIGHTobject.creationDate = datestr(now,30);
STRAIGHTobject.dataDirectory = path;
STRAIGHTobject.dataFileName = file;
STRAIGHTobject.samplingFrequency = fs;
STRAIGHTobject.waveform = x(:,1);
STRAIGHTobject.standAlone = true;
opt.f0ceil=500;
opt.f0floor=100;
STRAIGHTobject.F0extractionDate = datestr(now,30);
if 1 == 1

f0Structure = ...
    exF0candidatesTSTRAIGHTGB(STRAIGHTobject.waveform, ...
    STRAIGHTobject.samplingFrequency,opt);
end;


fs = STRAIGHTobject.samplingFrequency;
STRAIGHTobject.f0 = f0Structure.f0;
STRAIGHTobject.periodicityLevel = f0Structure.periodicityLevel;
STRAIGHTobject.temporalPositions = f0Structure.temporalPositions;
STRAIGHTobject.f0CandidatesMap = f0Structure.f0CandidatesMap;
STRAIGHTobject.f0CandidatesScoreMap = f0Structure.f0CandidatesScoreMap;
STRAIGHTobject.f0Structure = f0Structure;
STRAIGHTobject.f0candidatesPowerMap = f0Structure.f0candidatesPowerMap;
% 
% f0Structure=autoF0Tracking(STRAIGHTobject.F0Structure,x);
% f0Structure.vuv = refineVoicingDecision(x,f0Structure);

%%
GUIuserData=STRAIGHTobject;
refinedData = autoF0Tracking(GUIuserData,x);
refinedData.vuv = refineVoicingDecision(x,refinedData);
GUIuserData.f0 = refinedData.f0;
GUIuserData.vuv = refinedData.vuv;
GUIuserData.periodicityLevel = refinedData.periodicityLevel;

STRAIGHTobject.originalF0Structure = GUIuserData.f0Structure;
    STRAIGHTobject.refinedF0Structure = GUIuserData.f0Structure;
    STRAIGHTobject.refinedF0Structure.f0 = GUIuserData.f0;
    STRAIGHTobject.waveform = GUIuserData.waveform;
    if isfield(GUIuserData,'vuv')
        STRAIGHTobject.refinedF0Structure.vuv = GUIuserData.vuv;
    end;
    STRAIGHTobject.refinedF0Structure.periodicityLevel = ...
        GUIuserData.periodicityLevel;




%%
STRAIGHTobject.AperiodicityExtractionDate = datestr(now,30);

STRAIGHTobject.AperiodicityStructure = ...
    aperiodicityRatioSigmoid(STRAIGHTobject.waveform, ...
    STRAIGHTobject.refinedF0Structure,1,2,0); % originally 2

if isfield(STRAIGHTobject.refinedF0Structure,'vuv')
    STRAIGHTobject.AperiodicityStructure.vuv = ...
        STRAIGHTobject.refinedF0Structure.vuv;
   
end


STRAIGHTobject.SpectrumExtractionDate = datestr(now,30);
%prmIn.compensationCoefficient = [-0.2818;0.0397]; % this is experimental 15/Oct./2009
%prmIn.compensationCoefficient = -0.3804;
%prmIn.exponentControl = 1/10;
%prmIn.correctionForBlackman = 2.5;
%prmIn.debugperiodicityShaping = 2.5; % original 2.4
%prmIn.defaultF0 = 500;
SpectrumStructure = ...
    exSpectrumTSTRAIGHTGB(STRAIGHTobject.waveform, ...
    STRAIGHTobject.samplingFrequency, ...
    STRAIGHTobject.refinedF0Structure);%,prmIn);
%STRAIGHTobject.SpectrumStructure = ...
%    fixWindowEffects(SpectrumStructure,STRAIGHTobject.refinedF0Structure, ...
%    prmIn.correctionForBlackman,2,1/2); % 18/Oct./2009
STRAIGHTobject.SpectrumStructure = SpectrumStructure;
processedSpectrum = unvoicedProcessing(STRAIGHTobject);
STRAIGHTobject.SpectrumStructure.spectrogramSTRAIGHT = processedSpectrum;


if save==1
    STRAIGHTobject.lastUpdate = datestr(now);
    outFileName = ['StrObj' datestr(now,30) '.mat'];
    [file,path] = uiputfile(outFileName,'Save the STRAIGHT object');
    if length(file) == 1 && length(path) == 1
        if file == 0 || path == 0
            %okInd = 0;
            disp('Save is cancelled!');
            return;
        end;
    end;

    %pathReg = regexprep(path,'\s','\\ ');
    %eval(['save ' pathReg file ' STRAIGHTobject']);
    save('test.mat','STRAIGHTobject');
end


end

