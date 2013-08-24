function f_conversion( iFile,emotion,mixSp,workDir )
%F_CONVERSION Summary of this function goes here
%   Detailed explanation goes here

addpath('/home/lixian/speech/MatlabCode/STRAIGHT/TandemSTRAIGHTmonolithicPackage002Test');
prjdir=[pwd '/..'];
%  iFile=381;
% emotion='angry';
% method='tone_dall';
% durConv=1;
% f0Conv=1;
% spConv=1;
% ClusterSize=20;
% numMixOfPhrase=4;
% numMixOfTone=4;
% dctNumOfPhrase=4;
% workDir=[prjdir '/test/' method '/' emotion '/' num2str(numMixOfPhrase) 'mix' num2str(numMixOfTone) 'mix_' num2str(dctNumOfPhrase)];
% fprintf('generating %s/%d_%d%d%d.wav\n',workDir,iFile,durRate,f0Rate,spRate);

%%
ClusterSize=20;
%% read A
path=[prjdir '/neutral/wav/'];
file=[num2str(iFile) '.wav'];
SA=getSTRAIGHTobject(path,file,0);
AperiodicityStructure=SA.AperiodicityStructure;
SpectrumStructure =SA.SpectrumStructure;
path=[prjdir '/neutral/lab/phone/'];
file=[num2str(iFile) '.lab'];
labelA = readHTKLabel([path file]);

%%
durRatio=load([prjdir '/test/dur/' emotion '/F1/' num2str(iFile) '.predict_' num2str(ClusterSize)]);
durationModifier=ones(1,length(SpectrumStructure.temporalPositions));

for i=1:length(durRatio)
    start=labelA.segment(i,1)*1000/5+1;
    endx=labelA.segment(i,2)*1000/5+1;
    durationModifier(start:endx)=durRatio(i);
end

tempAp=AperiodicityStructure;
tempSp=SpectrumStructure;
for durConv=0:1
    for f0Conv=0:1
          for spConv=0:1
                fprintf('converting %s%d_c%d%d%d',workDir,iFile,durConv,f0Conv,spConv);
                AperiodicityStructure=tempAp;
                SpectrumStructure=tempSp;
                if(durConv)
                    temporalPositions=AperiodicityStructure.temporalPositions;
                    segmentDurations = diff(temporalPositions);
                    durationModifierReal = (durationModifier(1:end-1)+durationModifier(2:end))/2;
                    updatedTemporalPositions = ...
                    cumsum([temporalPositions(1); segmentDurations(:).*durationModifierReal(:)])';

                    AperiodicityStructure.temporalPositions = updatedTemporalPositions;
                    SpectrumStructure.temporalPositions = ...
                    AperiodicityStructure.temporalPositions;
                end

                %%
                if(f0Conv)
                     if exist([workDir '/' num2str(iFile) '.phrase.f0'])
                        temp=load([workDir '/' num2str(iFile) '.phrase.f0'])+load([workDir '/' num2str(iFile) '.tone.f0']);
                     
                     else
                         temp=load([workDir '/' num2str(iFile) '.tone.f0']);
                     end
                    
                % temp=load([directoryBase '/logf0/' num2str(iFile) '.f0' ]);
                    logf0=[temp ;0*ones(length(SpectrumStructure.temporalPositions)-length(temp),1)];
                    findvoiced=sign( diff([0; logf0(:)>0; 0]) );
                    voicedind=[find(findvoiced>0)'; find(findvoiced<0)'-1];
                   for k=1:size(voicedind,2)-1
                    logf0(voicedind(2,k):voicedind(1,k+1))=linspace(logf0(voicedind(2,k)),logf0(voicedind(1,k+1)),voicedind(1,k+1)-voicedind(2,k)+1);
                   end


                    first=find(logf0,1,'first');
                    last=find(logf0,1,'last');
                    logf0(1:first-1)=logf0(first);
                    logf0(last+1:end)=logf0(last);   
                    AperiodicityStructure.f0=exp(logf0);
                end

                %%
                if(spConv)
                 SpectrumStructure.spectrogramSTRAIGHT=load([prjdir '/vc/to' emotion '_' num2str(mixSp) 'mix' '/test/wav/neutral-' emotion  '_' num2str(mixSp) 'mix' '/' num2str(iFile) '.sp' ])';
                SpectrumStructure.spectrogramSTRAIGHT=load([prjdir '/vc/train/to' emotion '_ar/test/' num2str(128) '-mix' '/wav/convert-1' '/' num2str(iFile) '.sp' ])';

                end

                %%

                SynthesisStructure = ...
                    exTandemSTRAIGHTsynthNx(AperiodicityStructure,SpectrumStructure);
%                 soundsc(SynthesisStructure.synthesisOut,SynthesisStructure.samplingFrequency);
                wavwrite(SynthesisStructure.synthesisOut,SynthesisStructure.samplingFrequency,[workDir '/'  num2str(iFile) '_c' num2str(durConv) num2str(f0Conv) num2str(spConv) '.wav']);
                % SynthesisStructure = ...
                %     exGeneralSTRAIGHTsynthesisR2(AperiodicityStructure,SpectrumStructure);
                % soundsc(SynthesisStructure.synthesisOut,SynthesisStructure.samplingFrequency);
          end
    end
end
 end

