
function predict=gmm_ar_convert(x,l,order,mixNum,gmmFile,method)
% clear;
% pMix=3;
% sMix=6;
% trainDir=['../../train/tone_nodynamic/sad/'  num2str(pMix) 'mix'  num2str(sMix) 'mix_3'];
% gmmFile=[trainDir '/neutral_sad.tone.gmm_ar_full'];
% l=16;
% mixNum=sMix;
% order=1;
% method='full';

%% read model

if (strcmp(method,'full'))
 gmmP=importdata(gmmFile);


    w=gmmP(1:mixNum);
    miux=zeros(l/2,mixNum);
    miuy=zeros(l/2,mixNum);
    miu=zeros(l,mixNum);
    sigma=zeros(l,l,mixNum);
     ar=zeros(l,order,mixNum);
     arx=zeros(l/2,order,mixNum);
     ary=zeros(l/2,order,mixNum);
   % beginModel.w=w;
    for i=1:mixNum
        miux(:,i)=gmmP(mixNum+(i-1)*(l^2+l*(1+order))+1:mixNum+(i-1)*(l^2+ l*(1+order) )+l/2);
        miuy(:,i)=gmmP(mixNum+(i-1)*(l^2+l*(1+order) )+l/2+1:mixNum+(i-1)*(l^2+l*(1+order) )+l);
        for iO=1:order
            ar(:,iO,i)=gmmP(mixNum+(i-1)*(l^2+l*(1+order) )+l +l*(iO-1)+1 :mixNum+(i-1)*(l^2+l*(1+order) )+l*iO+l);
        end
        arx(:,iO,i)=ar(1:l/2,iO,i);
        ary(:,iO,i)=ar(l/2+1:l,iO,i);
        miu(:,i)=[miux(:,i);miuy(:,i)];
      %  beginModel.miu(:,i)=[miux(:,i);miuy(:,i)];
        sigma(:,:,i)=reshape(gmmP(mixNum+(i-1)*(l^2+l*(1+order) )+l*(1+order)+1: mixNum+(i-1)*(l^2+l*(1+order) )+l*(1+order)+l^2),l,l);
       % beginModel.sigma(:,:,i)=sigma(:,:,i);
    end

    sigmaxx=sigma(1:l/2,1:l/2,:);
    sigmayy=sigma(l/2+1:l,l/2+1:l,:);
    sigmaxy=sigma(1:l/2,l/2+1:l,:);
    sigmayx=sigma(l/2+1:l,1:l/2,:);
%% 
%x=[-0.0115096 0.085465 -0.0326061 0.00391863 -0.0131552 -0.00301318 0.00158419 0.00168599 0.0967152 -0.159517 0.029575 0.00595606 -0.016425 0.0107194 0.000881274 -0.00422446 -0.103581 0.177257 0.0704744 0.0171075 0.0113395 0.00681746 0.00332063 0.00574398 -0.0506981 -0.0123003 0.00256196 0.0156671 0.00075175 0.00924318 -0.00385466 0.00570197 -3.9615e-05 0.0414369 -0.00826747 0.0125788 -0.00141036 -0.000465152 0.00409928 -0.00359639 0.0422732 0.061477 -0.0160057 0.0195799 -0.0120588 0.00826337 -0.00393691 0.0013994 0.0288733 -0.0510451 -0.0121298 0.0153613 0.0179336 0.00302738 -0.00966935 -0.00851501]';
T=length(x)/l*2;
y=zeros(l/2,T);
postProb=zeros(mixNum,1);
liklihoodProb=0;
beginZero=zeros(l,1);
 preSegs=beginZero*ones(1,order);
 curMiux=zeros(l/2,mixNum);
for iT=1:T
    src=x((iT-1)*l/2+1:iT*l/2);
      sum_postProb=-1.0E10;
      maxMix=1;
       maxProb=-1.0E10;
    for iMix=1:mixNum
        curMiux(:,iMix)=miux(:,iMix);
        for iO=1:order
            curMiux(:,iMix)=curMiux(:,iMix)+arx(:,iO,iMix).*preSegs(1:l/2,iO);
        end
       postProb(iMix)= log(w(iMix))-l/2*log(2*pi)-0.5*log(det(sigmaxx(:,:,iMix)))-0.5.*(src-curMiux(:,iMix))'/sigmaxx(:,:,iMix)*(src-curMiux(:,iMix))          ;
        
            if maxProb<postProb(iMix)
                maxMix=iMix;
                maxProb=postProb(iMix);
            end
            sum_postProb =log_add( sum_postProb ,postProb(iMix));
       
             
    end
        postProb=exp(postProb-sum_postProb);
        
     %   liklihoodProb=liklihoodProb+sum_postProb;
    curMiuy=miuy(:,maxMix);
    for iO=1:order
        curMiuy=curMiuy+ary(:,iO,maxMix).*preSegs(l/2+1:l,iO);
    end
    
    y(:,iT)=curMiuy+sigmayx(:,:,maxMix)/sigmaxx(:,:,maxMix)*(src-curMiux(:,maxMix));
    preSegs(:,order)=[src;y(:,iT)];
    circshift(preSegs,1);
end

predict=reshape(y,T*l/2,1);
end

end