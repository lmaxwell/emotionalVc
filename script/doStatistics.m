
emotion='angry';
method='twotier';
prjDir='../';
dir=[prjDir 'train/' method '/' emotion '/'];

neutralPhrase=load([dir 'neutral.phrase']);
emotionPhrase=load([dir emotion '.phrase']);

figure(1)
subplot(2,2,1)
plot(neutralPhrase(:,1),neutralPhrase(:,2),'.');
hold on;
plot(emotionPhrase(:,1),emotionPhrase(:,2),'r.');
hold off;
axis([5 6 -0.5 0.5]);
subplot(2,2,2)
plot(neutralPhrase(:,2),neutralPhrase(:,3),'.');
hold on;
plot(emotionPhrase(:,2),emotionPhrase(:,3),'r.');
hold off;
axis([-0.5 0.5 -0.5 0.5]);
subplot(2,2,3)
plot(neutralPhrase(:,1),neutralPhrase(:,3),'.');
hold on;
plot(emotionPhrase(:,1),emotionPhrase(:,3),'r.');
hold off;
axis([5 6  -0.5 0.5]);
figure(2)
subplot(3,2,1)
boxplot(neutralPhrase(:,1));
subplot(3,2,2)
boxplot(emotionPhrase(:,1));
subplot(3,2,3)
boxplot(neutralPhrase(:,2));
subplot(3,2,4)
boxplot(emotionPhrase(:,2));
subplot(3,2,5)
boxplot(neutralPhrase(:,3));
subplot(3,2,6)
boxplot(emotionPhrase(:,3));

%%
figure(3);
subplot(2,3,1);
hist(neutralPhrase(:,1));
subplot(2,3,4);
hist(emotionPhrase(:,1));
subplot(2,3,2);
hist(neutralPhrase(:,2));
subplot(2,3,5);
hist(emotionPhrase(:,2));
subplot(2,3,3);
hist(neutralPhrase(:,3));
subplot(2,3,6);
hist(emotionPhrase(:,3));
%%
figure(4)
subplot(2,2,1);
plot(neutralPhrase(:,1),emotionPhrase(:,1),'.');
subplot(2,2,2);
plot(neutralPhrase(:,2),emotionPhrase(:,2),'.');
subplot(2,2,3);
plot(neutralPhrase(:,3),emotionPhrase(:,3),'.');



%% 
neutralTone=load([dir 'neutral.tone']);
emotionTone=load([dir emotion '.tone']);

figure(1)
subplot(2,2,1)
plot(neutralTone(:,1),neutralTone(:,2),'.');
hold on;
plot(emotionTone(:,1),emotionTone(:,2),'r.');
hold off;
axis([-1 1 -1 1]);
subplot(2,2,2)
plot(neutralTone(:,2),neutralTone(:,3),'.');
hold on;
plot(emotionTone(:,2),emotionTone(:,3),'r.');
hold off;
axis([-1 1 -1 1]);
subplot(2,2,3)
plot(neutralTone(:,1),neutralTone(:,3),'.');
hold on;
plot(emotionTone(:,1),emotionTone(:,3),'r.');
hold off;
axis([-1 1 -1 1]);
figure(2)
subplot(3,2,1)
boxplot(neutralTone(:,1));
subplot(3,2,2)
boxplot(emotionTone(:,1));
subplot(3,2,3)
boxplot(neutralTone(:,2));
subplot(3,2,4)
boxplot(emotionTone(:,2));
subplot(3,2,5)
boxplot(neutralTone(:,3));
subplot(3,2,6)
boxplot(emotionTone(:,3));

%%
figure(3);
subplot(2,3,1);
hist(neutralTone(:,1));
subplot(2,3,4);
hist(emotionTone(:,1));
subplot(2,3,2);
hist(neutralTone(:,2));
subplot(2,3,5);
hist(emotionTone(:,2));
subplot(2,3,3);
hist(neutralTone(:,3));
subplot(2,3,6);
hist(emotionTone(:,3));
%%
figure(4)
subplot(2,2,1);
plot(neutralTone(:,1),emotionTone(:,1),'.');
subplot(2,2,2);
plot(neutralTone(:,2),emotionTone(:,2),'.');
subplot(2,2,3);
plot(neutralTone(:,3),emotionTone(:,3),'.');

%%
[idx c sumd d]=kmeans(neutralTone,5);