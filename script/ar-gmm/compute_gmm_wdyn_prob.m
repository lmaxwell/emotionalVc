clear;
trainDir='../../train/tone_dall/sad/3mix6mix_3';
gmmFile=[trainDir '/neutral_sad.tone.gmm'];
trainFile=[trainDir '/neutral_sad.tone'];

l=32;
mixNum=6;
command=['x2x +fa' num2str(l) ' '  trainFile '.f >'  trainFile];
system(command);
%% read gmm model
 gmmP=importdata(gmmFile);


    w=gmmP(1:mixNum);
    miux=zeros(l/2,mixNum);
    miuy=zeros(l/2,mixNum);
    miu=zeros(l/2,mixNum);
    sigma=zeros(l,l,mixNum);
    beginModel.w=w;
    for i=1:mixNum
        miux(:,i)=gmmP(mixNum+(i-1)*(l^2+l)+1:mixNum+(i-1)*(l^2+l)+l/2);
        miuy(:,i)=gmmP(mixNum+(i-1)*(l^2+l)+l/2+1:mixNum+(i-1)*(l^2+l)+l);
        miu(:,i)=[miux(1:l/4,i);miuy(1:l/4,i)];
        beginModel.miu(:,i)=[miux(:,i);miuy(:,i)];
        sigma(:,:,i)=reshape(gmmP(mixNum+(i-1)*(l^2+l)+l+1:mixNum+(i-1)*(l^2+l)+l+l^2),l,l);
        beginModel.sigma(:,:,i)=sigma(:,:,i);
    end

    sigmaxx=sigma(1:l/4,1:l/4,:);
    sigmayy=sigma(l/2+1:l/2+l/4,l/2+1:l/2+l/4,:);
    sigmaxy=sigma(1:l/4,l/2+1:l/2+l/4,:);
    sigmayx=sigma(l/2+1:l/2+l/4,1:l/4,:);
    sigma_s=[sigmaxx sigmayx;sigmaxy sigmayy];


    %%
    trainData=importdata(trainFile)';
numOfSeg=size(trainData,2);
beginZero=zeros(l,1);



postProb=zeros(mixNum,1);

liklihoodProb=0;


start=0;

    for iSeg=1:numOfSeg
        fprintf(1,'%dth segment\n',iSeg);
        %set preSegs
    
        
      sum_postProb=-1.0E10;
        maxProb=-1.0E10;
        maxMix=1;
        for iMix=1:mixNum
            sample=[trainData(1:l/4,start+iSeg);trainData(l/2+1:l/2+l/4,start+iSeg)];
            postProb(iMix)= log(w(iMix))-l/2*log(2*pi)-0.5*log(det(sigma_s(:,:,iMix)))-0.5.*(sample-miu(:,iMix))'/sigma_s(:,:,iMix)*(sample-miu(:,iMix))          ;
            %.*1./((2*pi)^(l/2).*sqrt(det(sigma(:,:,iMix)))).* ...
         %exp(-0.5.*(trainData(:,start+iSeg)-miu(:,iMix))'/sigma(:,:,iMix)*(trainData(:,start+iSeg)-miu(:,iMix)));
     %   log(w(i))-0.5*log(det(sigmaxx(:,:,i)))-4*log(2*pi)-0.5.*(x-miux(:,i))'/sigmaxx(:,:,i)*(x-miux(:,i));
            sum_postProb =log_add( sum_postProb ,postProb(iMix));
            if maxProb<postProb(iMix)
                maxMix=iMix;
                maxProb=postProb(iMix);
            end
        end
       
        postProb=exp(postProb-sum_postProb);
     
        liklihoodProb=liklihoodProb+sum_postProb;
      
      
        
    end

    

    


avgProb=liklihoodProb/numOfSeg;
fprintf(1,'%f\n',avgProb);
   