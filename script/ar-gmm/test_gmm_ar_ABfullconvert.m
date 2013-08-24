  
clear;
pMix=6;
sMix=8;
trainDir=['../../train/tone_nodynamic_diag/sad/'  num2str(pMix) 'mix'  num2str(sMix) 'mix_3'];
gmmFile=[trainDir '/neutral_sad.tone.gmm_ar_ABfull'];
l=16;
mixNum=sMix;
orderx=1;
ordery=1;
addpath('ar-gmm/lightspeed');
gmmP=importdata(gmmFile);


w=gmmP(1:mixNum);
    miux=zeros(l/2,mixNum);
    miuy=zeros(l/2,mixNum);
    miu=zeros(l,mixNum);
    sigma=zeros(l,l,mixNum);
    var=zeros(l,mixNum);
     ar=zeros(l,ordery,mixNum);
     arx=zeros(l/2,orderx,mixNum);
     ary=zeros(l/2,ordery,mixNum);
     bt=zeros(l/2,l/2,mixNum);
     gconst=zeros(mixNum,1);
  for i=1:mixNum
        miux(:,i)=gmmP(mixNum+(i-1)*(l^2+l/2*(1+orderx) +l/2*(1+ordery)+l/2*l/2 )+1:mixNum+(i-1)*(l^2+ l/2*(1+orderx) +l/2*(1+ordery)+l/2*l/2 )+l/2);
        miuy(:,i)=gmmP(mixNum+(i-1)*(l^2+l/2*(1+orderx) +l/2*(1+ordery) +l/2*l/2)+l/2+1:mixNum+(i-1)*(l^2+l/2*(1+orderx) +l/2*(1+ordery)+l/2*l/2 )+l);
        for iO=1:orderx
            arx(:,iO,i)=gmmP(mixNum+(i-1)*(l^2+l/2*(1+orderx) +l/2*(1+ordery) +l/2*l/2 )+l + l/2*(iO-1)+1 :mixNum+(i-1)*(l^2+l/2*(1+orderx) +l/2*(1+ordery) +l/2*l/2 )+l+ l/2*(iO-1)+l/2);
        end
         for iO=1:ordery
            ary(:,iO,i)=gmmP(mixNum+(i-1)*(l^2+l/2*(1+orderx) +l/2*(1+ordery) +l/2*l/2 )+l+l/2*orderx + l/2*(iO-1)+ 1 :mixNum+(i-1)*(l^2+l/2*(1+orderx) +l/2*(1+ordery) +l/2*l/2 )+ l+l/2*orderx+ l/2*(iO-1)+l/2);
            ar(:,iO,i)=[arx(:,iO,i); ary(:,iO,i)];
         end
        
       bt(:,:,i)=reshape(gmmP(mixNum+(i-1)*(l^2+l/2*(1+orderx) +l/2*(1+ordery)+l/2*l/2 )+ l +l/2*orderx + l/2*ordery+  1 :mixNum+(i-1)*(l^2+l/2*(1+orderx) +l/2*(1+ordery)+l/2*l/2 )+ l + l/2*orderx+ l/2*ordery+l/2*l/2 ),l/2,l/2);
        miu(:,i)=[miux(:,i);miuy(:,i)];
      %  beginModel.miu(:,i)=[miux(:,i);miuy(:,i)];
        sigma(:,:,i)=reshape(gmmP(mixNum+(i-1)*(l^2+l/2*(1+orderx) +l/2*(1+ordery)+l/2*l/2 )+ l/2*(1+orderx) +l/2*(1+ordery)+l/2*l/2 + 1: mixNum+(i-1)*(l^2+l/2*(1+orderx) +l/2*(1+ordery)+l/2*l/2 )+  l/2*(1+orderx) +l/2*(1+ordery) + l/2*l/2+l^2),l,l);
       % beginModel.sigma(:,:,i)=sigma(:,:,i);
     %  sigma(:,:,i)=[sigma(1:l/2,1:l/2,i),sigma(1:l/2,1:l/2,i)*bt(:,:,i)'; bt(:,:,i)*sigma(1:l/2,1:l/2,i),bt(:,:,i)*sigma(1:l/2,1:l/2,i)*bt(:,:,i)'+sigma(l/2+1:l,l/2+1:l,i) ];
        var(:,i)=diag(sigma(:,:,i));
         gconst(i)=log(w(i))-l/2*log(2*pi)-0.5*log(det( sigma(:,:,i) ));
  end

  %%
  x=[-0.0115096 0.085465 -0.0326061 0.00391863 -0.0131552 -0.00301318 0.00158419 0.00168599 0.0967152 -0.159517 0.029575 0.00595606 -0.016425 0.0107194 0.000881274 -0.00422446 -0.103581 0.177257 0.0704744 0.0171075 0.0113395 0.00681746 0.00332063 0.00574398 -0.0506981 -0.0123003 0.00256196 0.0156671 0.00075175 0.00924318 -0.00385466 0.00570197 -3.9615e-05 0.0414369 -0.00826747 0.0125788 -0.00141036 -0.000465152 0.00409928 -0.00359639 0.0422732 0.061477 -0.0160057 0.0195799 -0.0120588 0.00826337 -0.00393691 0.0013994 0.0288733 -0.0510451 -0.0121298 0.0153613 0.0179336 0.00302738 -0.00966935 -0.00851501]';
T=length(x)/l*2;
y=zeros(l/2,T);
postProb=zeros(mixNum,1);
liklihoodProb=0;
beginZero=zeros(l,1);
 preSegs=beginZero*ones(1,ordery);
 curMiux=zeros(l/2,mixNum);
for iT=1:T
    src=x((iT-1)*l/2+1:iT*l/2);

     
        curMiu=miu;
        for iO=1:ordery
            curMiu(1:l,:)=curMiu(1:l,:) + preSegs(1:l,iO)*ones(1,mixNum).*squeeze( ar(:,iO,:) );
            
        end
        
           for iMix=1:mixNum
            curMiu(l/2+1:l,iMix)=curMiu(l/2+1:l,iMix) + bt(:,:,iMix)*src;
        
      
   %   postProb(iMix)= log(w(iMix))-l/2*log(2*pi)-0.5*log(det(sigma(1:l,1:l,iMix)))-0.5.*([src;curMiu(l/2+1:l,iMix)]-curMiu(:,iMix))'/sigma(1:l,1:l,iMix)*([src;curMiu(l/2+1:l,iMix)]-curMiu(:,iMix))          ;
        
          postProb(iMix)= log(w(iMix))-l/2*log(2*pi)-0.5*log(det(sigma(1:l/2,1:l/2,iMix)))-0.5.*(src-curMiu(1:l/2,iMix))'/sigma(1:l/2,1:l/2,iMix)*(src-curMiu(1:l/2,iMix))          ;
       
             
    end
      %  postProb=gconst- ndsum(0.5.*([src*ones(1,mixNum);curMiu(l/2+1:l,:)]-curMiu).*([src*ones(1,mixNum);curMiu(l/2+1:l,:)]-curMiu)./var, 1) ;
          maxMix=find(postProb==max(postProb));
         y(:,iT)=curMiu(l/2+1:l, maxMix );
         
         
    if(ordery>0)     
    preSegs(:,ordery)=[src;y(:,iT)];
    
    circshift(preSegs,1);
    end
end

predict=reshape(y,T*l/2,1);