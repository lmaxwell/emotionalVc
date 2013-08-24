iFile=421;
method='tone_dall';
emotion='sad';
numMixOfPhrase=1;
numMixOfTone=1;
dctNumOfPhrase=3;

workDir=['../test/' method '/' emotion '/' num2str(numMixOfPhrase) 'mix' num2str(numMixOfTone) 'mix' '_' num2str(dctNumOfPhrase) '/'];

PhraseF0File=[workDir num2str(iFile) '.phrase.f0'];
ToneF0File=[workDir  num2str(iFile) '.tone.f0'];



Phrase=load(PhraseF0File);
Tone=load(ToneF0File);

F0=Phrase+Tone;
F0(Tone==0)=0;
time=0:0.005:(length(F0)-1)*0.005;
voicedTime=time(F0>0);
voicedF0=F0(F0>0);
voicedTone=Tone(Tone~=0);
phraseTime=time(Phrase>0);
Phrase=Phrase(Phrase>0);
figure;
subplot(2,1,1);
plot(phraseTime,Phrase,'.');
hold on;
axis([0 time(end)+0.2 4.9 6]);
plot(voicedTime,voicedF0,'r.');
subplot(2,1,2);

plot(voicedTime,voicedTone,'*');
axis([0 time(end)+0.2 4.9 6]);

%%
iFile=421;
method='gmm';
emotion='sad';
numMixOfPhrase=7;
numMixOfTone=8;
dctNumOfPhrase=4;

workDir=['../test/' method '/' emotion '/' num2str(numMixOfPhrase) 'mix' num2str(numMixOfTone) 'mix' '_' num2str(dctNumOfPhrase) '/'];


ToneF0File=[workDir  num2str(iFile) '.tone.f0'];



Tone=load(ToneF0File);

F0=Tone;
F0(Tone==0)=0;
time=0:0.005:(length(F0)-1)*0.005;
voicedTime=time(F0>0);
voicedF0=F0(F0>0);
voicedTone=Tone(Tone~=0);
figure;
subplot(2,1,1);


plot(voicedTime,voicedF0,'r.');
axis([0 time(end)+0.2 4.9 6]);
subplot(2,1,2);

plot(voicedTime,voicedTone,'*');
axis([0 time(end)+0.2 4.9 6]);


%%
%基频分解示意
iFile=244;
emotion='neutral';
dctNum=8;
 
phonelabfile=['../' emotion '/lab/phone/' num2str(iFile) '.lab' ];
phonelable=readHTKLabel2(phonelabfile);

wordlabfile=['../' emotion '/lab/word/' num2str(iFile) '.lab' ];
wordlable=readHTKLabel2(wordlabfile);


workDir=['../test/error_dctnum/' emotion  '/f0'  ];
originf0=load([workDir '/' num2str(iFile) '.f0']);
ophrase_f0=load([workDir '/' num2str(iFile) '.phrase_f0']);
otone_f0=load([workDir '/' num2str(iFile) '.tone_f0']);


ctone_f0=load([workDir '/../' num2str(dctNum) '/' num2str(iFile) '.tone_f0']);

reconstructf0=load([workDir '/../' num2str(dctNum) '/' num2str(iFile) '.ref0']);
time=0:0.005:(length(otone_f0)-1)*0.005;
voicedTime=time(otone_f0~=0);
phraseTime=time(ophrase_f0~=00);
ophrase_f0=ophrase_f0(ophrase_f0~=0);
otone_f0=otone_f0(otone_f0~=0);
ctone_f0=ctone_f0(ctone_f0~=0);
figure;

subplot(3,1,1);
hold on;
box off;
plot(voicedTime,originf0,'k');
hold on;
plot(voicedTime,reconstructf0,'k+');
hold  on;
legend('original f0','reconstruct f0');
for i=1:length(phonelable.segment(:,1))
plot([phonelable.segment(i,1) phonelable.segment(i,1)],[4.5 6],'k');
text((phonelable.segment(i,1)+phonelable.segment(i,2))/2,4.7,phonelable.phoneme(i,:));

hold on;
end
hold on;
plot([0,1.8], [6 ,6],'k');
hold on;
plot([1.8,1.8],[4.5,6],'k');
hold off;
title('F0');
set(gca,'xticklabel',sprintf('%1.1f|',get(gca,'xtick')));
set(gca,'yticklabel',sprintf('%1.1f|',get(gca,'ytick')));
subplot(3,1,2);
hold on;
box off;
plot(phraseTime,ophrase_f0,'k*');
hold on;
for i=1:length(wordlable.segment(:,1))
plot([wordlable.segment(i,1) wordlable.segment(i,1)],[-0.5 0.5],'k');
text((wordlable.segment(i,1)+wordlable.segment(i,2))/2,-0.3,wordlable.phoneme(i,:));

hold on;
end
hold on;
plot([0,1.8], [0.5,0.5],'k');
hold on;
plot([1.8,1.8],[-0.5,0.5],'k');

title('phrase level');
set(gca,'xticklabel',sprintf('%1.1f|',get(gca,'xtick')));
set(gca,'yticklabel',sprintf('%1.1f|',get(gca,'ytick')));
subplot(3,1,3);
hold on;
box off;
plot(voicedTime,ctone_f0,'k*')
hold on;
for i=1:length(phonelable.segment(:,1))
plot([phonelable.segment(i,1) phonelable.segment(i,1)],[-1 0.5],'k');

text((phonelable.segment(i,1)+phonelable.segment(i,2))/2,-0.8,phonelable.phoneme(i,:));
hold on;
end
hold on;
plot([0,1.8], [0.5,0.5],'k');
hold on;
plot([1.8,1.8],[-1,0.5],'k');

title('syllable level');
set(gca,'xticklabel',sprintf('%1.1f|',get(gca,'xtick')));
set(gca,'yticklabel',sprintf('%1.1f|',get(gca,'ytick')));

%% 
%重建误差图