clear;
nMix=128;
emotion='sad';
trainDir=['../../vc/train/to' emotion '_ar/gmm_jde/' num2str(nMix) '-mix'];
gmmMean=[trainDir '/neutral-' emotion num2str(nMix) '.mean.txt' ];
gmmCov=[trainDir '/neutral-' emotion num2str(nMix) '.cov.txt' ];
gmmWght=[trainDir '/neutral-' emotion num2str(nMix) '.wght.txt' ];

trainFile=[trainDir '/../../dtw/' num2str(nMix) '-mix/neutral-' emotion '.mat-ar.txt'];
l=48;
%nMix=10;
mixNum=nMix;
orderx=0;
ordery=1;
vfloor=0.001;

logfd=fopen([trainDir '/log_xy_ABdiag_order' num2str(ordery)],'w');
%% read gmm model
 w=importdata(gmmWght);
 miu=reshape( importdata(gmmMean) ,l,mixNum);
sigma=reshape( importdata(gmmCov),l,l,mixNum);

miux=miu(1:l/2,:);
miuy=miu(l/2+1:l,:);
ar=zeros(l,ordery,mixNum);
ar2=zeros(l,ordery,mixNum);

preMiu=zeros(l,ordery,mixNum);
gconst=zeros(mixNum,1);
invSigma=zeros(l,l,mixNum);
for iMix=1:mixNum
    
   
    gconst(iMix)=log(w(iMix))-l/2*log(2*pi)-0.5*log(det(sigma(:,:,iMix)));
    
    invSigma(:,:,iMix)=inv( sigma(:,:,iMix) );
end 

   
    % set ar coefficient
    ar=zeros(l,ordery,mixNum);
   arx=zeros(l/2,orderx,mixNum);
ary=zeros(l/2,ordery,mixNum);
bt=zeros(l/2,mixNum);
    
    
%% train ar combined gmm

%first run  
trainData=importdata(trainFile)';

%%
% Data=trainData(:,trainData(1,:)~=0);
% %Data=Data(:,1:10000);
% nbStates=mixNum;
% [nbVar, nbData] = size(Data);
% 
% [Data_id, Centers] = kmeans(Data', nbStates);
% Mu = Centers';
% for i=1:nbStates
%   idtmp = find(Data_id==i);
%   Priors(i) = length(idtmp);
%   Sigma(:,:,i) = cov([Data(:,idtmp) Data(:,idtmp)]');
%   %% Add a tiny variance to avoid numerical instability
%  % Sigma(:,:,i) = Sigma(:,:,i) + 1E-5.*diag(ones(nbVar,1));
% end
% Priors = Priors ./ sum(Priors);
% 
% %%
% w=Priors';
% miu=Mu;
% sigma=Sigma;
%  gconst=zeros(mixNum,1);
% invSigma=zeros(l,l,mixNum);
% for iMix=1:mixNum
%     
%    
%     gconst(iMix)=log(w(iMix))-l/2*log(2*pi)-0.5*log(det(sigma(:,:,iMix)));
%     
%     invSigma(:,:,iMix)=inv( sigma(:,:,iMix) );
% end 
%%
indexOfSen=find(trainData(1,:)==0);
numOfSen=length(indexOfSen);
indexOfSen=[indexOfSen  size(trainData,2)+1];
beginZero=zeros(l,1);

A=zeros(l/2,orderx+1,orderx+1,mixNum);
C=zeros(l/2,orderx+1,mixNum);
B=zeros(l/2,ordery+2,ordery+2,mixNum);
D=zeros(l/2,ordery+2,mixNum);
mean=zeros(l,1);

postProb=zeros(mixNum,1);

liklihoodProb=0;
numOfSeg=0;
accW=zeros(mixNum,1);
accSig=zeros(l,mixNum);
accSig1=zeros(l,mixNum);
accSig2=zeros(l,ordery,mixNum);

accSigy3=zeros(l/2,mixNum);
vFloor=zeros(l,1);

for iSen=1:numOfSen
    start=indexOfSen(iSen)+1;
    xend=indexOfSen(iSen+1)-1;
    fprintf( 1 , '%dth sentence,start:%d,end:%d\n' , iSen , start , xend);
    for iSeg=0:xend-start
    %    fprintf(1,'%dth segment\n',iSeg);
        %set preSegs
        if iSeg==0
              preSegs=beginZero*ones(1,ordery-iSeg);
        elseif(iSeg<ordery)
           
            preSegs=[beginZero*ones(1,ordery-iSeg) trainData(:,start+iSeg-1:-1:start)];
           
        else
            preSegs=trainData(:,start+iSeg-1:-1:start+iSeg-ordery);
       end
        
        sum_postProb=-1.0E10;
        maxProb=-1.0E10;
        maxMix=1;
        for iMix=1:mixNum
            
          postProb(iMix)= gconst(iMix)-0.5.*(trainData(:,start+iSeg)-miu(:,iMix))'*invSigma(:,:,iMix)*(trainData(:,start+iSeg)-miu(:,iMix))          ;
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
         accW=accW+postProb;
        liklihoodProb=liklihoodProb+sum_postProb;
        for iMix=1:mixNum 
           
            for iOd=1:ordery  
                for jOd=1:ordery 
              
                
                 B(  :,iOd,jOd , iMix)= B(  :,iOd,jOd ,iMix)+ postProb(iMix)*preSegs(l/2+1:l,iOd).*preSegs(l/2+1:l,jOd);  %iOd th function , jOd th parameter
                end
              
              
                
                D(:,iOd,iMix)=D(:,iOd,iMix)+postProb(iMix)*trainData(l/2+1:l,start+iSeg).*preSegs(l/2+1:l,iOd);  %iOd th function
                B(:,ordery+1,iOd,iMix)=B(:,ordery+1,iOd,iMix)+postProb(iMix).*preSegs(l/2+1:l,iOd);  % order+1 th fuction, iOd th parameter
                B(:,ordery+2,iOd,iMix)=B(:,ordery+2,iOd,iMix)+postProb(iMix).*trainData(1:l/2,start+iSeg).*preSegs(l/2+1:l,iOd);
                B(:,iOd,ordery+1,iMix)=B(:,iOd,ordery+1,iMix)+ postProb(iMix).*preSegs(l/2+1:l,iOd);         % iOd th function , order+1 th parameter
               B(:,iOd,ordery+2,iMix)=B(:,iOd,ordery+2,iMix)+ postProb(iMix).*trainData(1:l/2,start+iSeg).*preSegs(l/2+1:l,iOd);  
               
               accSig2(:,iOd,iMix)=accSig2(:,iOd,iMix)+postProb(iMix)*trainData(1:l,start+iSeg).*preSegs(1:l,iOd);
            end
               A(:,orderx+1,orderx+1,iMix)=   A(:,orderx+1,orderx+1,iMix)+postProb(iMix);   %order+1 th function ,order+1 th parameter
               C(:,orderx+1,iMix)=C(:,orderx+1,iMix)+trainData(1:l/2,start+iSeg).*postProb(iMix); % order+1th function
                
              B(:,ordery+1,ordery+1,iMix)=   B(:,ordery+1,ordery+1,iMix)+postProb(iMix);   %order+1 th function ,order+1 th parameter
               D(:,ordery+1,iMix)=D(:,ordery+1,iMix)+trainData(l/2+1:l,start+iSeg).*postProb(iMix); % order+1th function
                B(:,ordery+1,ordery+2,iMix)=   B(:,ordery+1,ordery+2,iMix)+postProb(iMix).*trainData(1:l/2,start+iSeg);  
                  B(:,ordery+2,ordery+1,iMix)=   B(:,ordery+2,ordery+1,iMix)+postProb(iMix).*trainData(1:l/2,start+iSeg);
                    B(:,ordery+2,ordery+2,iMix)=   B(:,ordery+2,ordery+2,iMix)+postProb(iMix).*trainData(1:l/2,start+iSeg).*trainData(1:l/2,start+iSeg);
               D(:,ordery+2,iMix)=D(:,ordery+2,iMix)+trainData(l/2+1:l,start+iSeg).*trainData(1:l/2,start+iSeg).*postProb(iMix); 
               
               
%                     accSigx(:,iMix)=accSigx(:,iMix)+postProb(iMix)*trainData(1:l/2,start+iSeg).*trainData(1:l/2,start+iSeg);
%                
%                    accSigx1(:,iMix)=accSigx1(:,iMix)+postProb(iMix)*trainData(1:l/2,start+iSeg);
%                    
%                   accSigy(:,iMix)=accSigy(:,iMix)+postProb(iMix)*trainData(l/2+1:l,start+iSeg).*trainData(l/2+1:l,start+iSeg);
%                   accSigy1(:,iMix)=accSigy1(:,iMix)+postProb(iMix)*trainData(l/2+1:l,start+iSeg);
%          for iO=1:order
%                 accSigx2(:,iO,iMix)=accSigx2(:,iO,iMix)+postProb(iMix)*trainData(1:l/2,start+iSeg).*preSegs(1:l/2,iO);
%                 accSigy2(:,iO,iMix)=accSigy2(:,iO,iMix)+postProb(iMix)*trainData(l/2+1:l,start+iSeg).*preSegs(l/2+1:l,iO);
%          end
%                 
%          accSigy3(:,iMix)=accSigy3(:,iMix)+postProb(iMix).*trainData(l/2+1:l,start+iSeg).*trainData(1:l/2,start+iSeg);
   
                  accSig(:,iMix)=accSig(:,iMix)+postProb(iMix)*trainData(1:l,start+iSeg).*trainData(1:l,start+iSeg);
               
                   accSig1(:,iMix)=accSig1(:,iMix)+postProb(iMix)*trainData(1:l,start+iSeg);
                   
               
        
                
          
         
                
         accSigy3(:,iMix)=accSigy3(:,iMix)+postProb(iMix).*trainData(l/2+1:l,start+iSeg).*trainData(1:l/2,start+iSeg);
         
         
        end
        
        
         mean(:)=mean(:)+trainData(:,start+iSeg);  %for Sigma update
        
        vFloor=vFloor+(trainData(:,start+iSeg)-miu(:,maxMix)).* (trainData(:,start+iSeg)-miu(:,maxMix));
       
        numOfSeg=numOfSeg+1;
        
    end
    mean=mean/numOfSeg;
    

    
    
end
vFloor=vFloor*vfloor/numOfSeg;
    for iMix=1:mixNum
        for iL=1:l/2
              tempA=squeeze(A(iL,:,:,iMix));
            if min(svd(tempA))<1e-7
                temp=pinv(tempA)* C(iL,:,iMix)';
            else
                temp=tempA \ C(iL,:,iMix)';
            end
           
            arx(iL,:,iMix)=temp(1:orderx);
            miux(iL,iMix)=temp(orderx+1);
        end
        for iL=1:l/2
            tempB=squeeze( B(iL,:,:,iMix) );
             if min(svd(tempB))<1e-7
                temp=pinv(tempB)* D(iL,:,iMix)';
            else
                temp=tempB \ D(iL,:,iMix)';
             end
            ary(iL,:,iMix)=temp(1:ordery);
            miuy(iL,iMix)=temp(ordery+1);
            bt(iL,iMix)=temp(ordery+2);
               %bt(iL,iMix)=0;
        end
        
       w(iMix)=accW(iMix)/numOfSeg;
         sigmaxx(:,:,iMix)=diag(accSig(1:l/2,iMix)-accSig1(1:l/2,iMix).*miux(:,iMix))./accW(iMix);
         
          sigmayy(:,:,iMix)=diag(accSig(l/2+1:l,iMix)-accSig1(l/2+1:l,iMix).*miuy(:,iMix)  -accSigy3(:,iMix).*bt(:,iMix)        )./accW(iMix);
          
     for iO=1:ordery  
           
          %  sigmaxx(:,:,iMix)=sigmaxx(:,:,iMix)-diag(accSig2(1:l/2,iO,iMix).*arx(:,iO,iMix) )./accW(iMix);
            sigmayy(:,:,iMix)=sigmayy(:,:,iMix)-diag( accSig2(l/2+1:l,iO,iMix).*ary(:,iO,iMix)  )./accW(iMix);

     end
     sigma(:,:,iMix)=blkdiag(sigmaxx(:,:,iMix),sigmayy(:,:,iMix));
      [D_,p_]=chol(sigma(:,:,iMix));
      if p_>0
        sigma(:,:,iMix)=sigma(:,:,iMix)+ diag(vFloor);
      end
      
        gconst(iMix)=log(w(iMix))-l/2*log(2*pi)-0.5*log(det(sigma(:,:,iMix)));
    invSigma(:,:,iMix)=inv( sigma(:,:,iMix) );
      
      
    end

ogmm=fopen([trainDir '/gmm_ar_diag_1thit_order' num2str(ordery) ],'w');
for iMix=1:mixNum
    fprintf(ogmm,'%f\n',w(iMix)); 
end
for iMix=1:mixNum
    fprintf(ogmm,'%f\n',miu(:,iMix));
    for iO=1:ordery
        fprintf(ogmm,'%f\n',ar(:,ordery,iMix));
    end
    fprintf(ogmm,'%f\n',sigma(:,:,iMix));
end
%%
preProb=liklihoodProb/numOfSeg;
curProb=0;
 fprintf(1,'initialProb %f\n',preProb);
%iteration
while(abs((curProb-preProb)/preProb)>1e-5)
A=zeros(l/2,orderx+1,orderx+1,mixNum);
C=zeros(l/2,orderx+1,mixNum);
B=zeros(l/2,ordery+2,ordery+2,mixNum);
D=zeros(l/2,ordery+2,mixNum);
mean=zeros(l,1);

postProb=zeros(mixNum,1);

liklihoodProb=0;
numOfSeg=0;
accW=zeros(mixNum,1);
accSig=zeros(l,mixNum);
accSig1=zeros(l,mixNum);
accSig2=zeros(l,ordery,mixNum);

accSigy3=zeros(l/2,mixNum);

for iSen=1:numOfSen
    start=indexOfSen(iSen)+1;
    xend=indexOfSen(iSen+1)-1;
  %  fprintf( 1 , '%dth sentence,start:%d,end:%d\n' , iSen , start , xend);
    for iSeg=0:xend-start
 %       fprintf(1,'%dth segment\n',iSeg);
        %set preSegs
         if iSeg==0
              preSegs=beginZero*ones(1,ordery-iSeg);
        elseif(iSeg<ordery)
           
            preSegs=[beginZero*ones(1,ordery-iSeg) trainData(:,start+iSeg-1:-1:start)];
           
        else
            preSegs=trainData(:,start+iSeg-1:-1:start+iSeg-ordery);
        end
        
        sum_postProb=-1.0E10;
       
        for iMix=1:mixNum

            if iSeg<0
                curMiu=beginModel.miu(:,iMix);
                
                 postProb(iMix)= w(iMix).*1./((2*pi)^(l/2).*sqrt(det(beginModel.sigma(:,:,iMix)))).* ...
            exp(-0.5.*(trainData(:,start+iSeg)-curMiu)'/beginModel.sigma(:,:,iMix)*(trainData(:,start+iSeg)-curMiu));
            else
                curMiu=[miux(:,iMix);miuy(:,iMix)];
                for iO=1:ordery
                %    curMiu(1:l/2)=curMiu(1:l/2)+arx(:,iO,iMix).*preSegs(1:l/2,iO);
                    curMiu(l/2+1:l)=curMiu(l/2+1:l)+ary(:,iO,iMix).*preSegs(l/2+1:l,iO);
                  
                end
                   curMiu(l/2+1:l)=curMiu(l/2+1:l)+bt(:,iMix).*trainData(1:l/2,start+iSeg);
                
               postProb(iMix)=gconst(iMix)-0.5.*(trainData(:,start+iSeg)-curMiu)'*invSigma(:,:,iMix)*(trainData(:,start+iSeg)-curMiu)          ;
            %.*1./((2*pi)^(l/2).*sqrt(det(sigma(:,:,iMix)))).* ...
         %exp(-0.5.*(trainData(:,start+iSeg)-miu(:,iMix))'/sigma(:,:,iMix)*(trainData(:,start+iSeg)-miu(:,iMix)));
     %   log(w(i))-0.5*log(det(sigmaxx(:,:,i)))-4*log(2*pi)-0.5.*(x-miux(:,i))'/sigmaxx(:,:,i)*(x-miux(:,i));
            end
            sum_postProb =log_add( sum_postProb ,postProb(iMix));
         %       postProb(iMix)= w(iMix).*1./((2*pi)^(l/2).*sqrt(det(sigma(:,:,iMix)))).* ...
           % exp(-0.5.*(trainData(:,start+iSeg)-curMiu)'/sigma(:,:,iMix)*(trainData(:,start+iSeg)-curMiu));
            
             
        end
        postProb=exp(postProb-sum_postProb);
        accW=accW+postProb;
        liklihoodProb=liklihoodProb+sum_postProb;
        
         for iMix=1:mixNum 
           
            for iOd=1:ordery  
                for jOd=1:ordery 
              
                
                 B(  :,iOd,jOd , iMix)= B(  :,iOd,jOd ,iMix)+ postProb(iMix)*preSegs(l/2+1:l,iOd).*preSegs(l/2+1:l,jOd);  %iOd th function , jOd th parameter
                end
              
              
                
                D(:,iOd,iMix)=D(:,iOd,iMix)+postProb(iMix)*trainData(l/2+1:l,start+iSeg).*preSegs(l/2+1:l,iOd);  %iOd th function
                B(:,ordery+1,iOd,iMix)=B(:,ordery+1,iOd,iMix)+postProb(iMix).*preSegs(l/2+1:l,iOd);  % order+1 th fuction, iOd th parameter
                B(:,ordery+2,iOd,iMix)=B(:,ordery+2,iOd,iMix)+postProb(iMix).*trainData(1:l/2,start+iSeg).*preSegs(l/2+1:l,iOd);
                B(:,iOd,ordery+1,iMix)=B(:,iOd,ordery+1,iMix)+ postProb(iMix).*preSegs(l/2+1:l,iOd);         % iOd th function , order+1 th parameter
               B(:,iOd,ordery+2,iMix)=B(:,iOd,ordery+2,iMix)+ postProb(iMix).*trainData(1:l/2,start+iSeg).*preSegs(l/2+1:l,iOd);  
               
               accSig2(:,iOd,iMix)=accSig2(:,iOd,iMix)+postProb(iMix)*trainData(1:l,start+iSeg).*preSegs(1:l,iOd);
            end
               A(:,orderx+1,orderx+1,iMix)=   A(:,orderx+1,orderx+1,iMix)+postProb(iMix);   %order+1 th function ,order+1 th parameter
               C(:,orderx+1,iMix)=C(:,orderx+1,iMix)+trainData(1:l/2,start+iSeg).*postProb(iMix); % order+1th function
                
              B(:,ordery+1,ordery+1,iMix)=   B(:,ordery+1,ordery+1,iMix)+postProb(iMix);   %order+1 th function ,order+1 th parameter
               D(:,ordery+1,iMix)=D(:,ordery+1,iMix)+trainData(l/2+1:l,start+iSeg).*postProb(iMix); % order+1th function
                B(:,ordery+1,ordery+2,iMix)=   B(:,ordery+1,ordery+2,iMix)+postProb(iMix).*trainData(1:l/2,start+iSeg);  
                  B(:,ordery+2,ordery+1,iMix)=   B(:,ordery+2,ordery+1,iMix)+postProb(iMix).*trainData(1:l/2,start+iSeg);
                    B(:,ordery+2,ordery+2,iMix)=   B(:,ordery+2,ordery+2,iMix)+postProb(iMix).*trainData(1:l/2,start+iSeg).*trainData(1:l/2,start+iSeg);
               D(:,ordery+2,iMix)=D(:,ordery+2,iMix)+trainData(l/2+1:l,start+iSeg).*trainData(1:l/2,start+iSeg).*postProb(iMix); 
               
               
%                     accSigx(:,iMix)=accSigx(:,iMix)+postProb(iMix)*trainData(1:l/2,start+iSeg).*trainData(1:l/2,start+iSeg);
%                
%                    accSigx1(:,iMix)=accSigx1(:,iMix)+postProb(iMix)*trainData(1:l/2,start+iSeg);
%                    
%                   accSigy(:,iMix)=accSigy(:,iMix)+postProb(iMix)*trainData(l/2+1:l,start+iSeg).*trainData(l/2+1:l,start+iSeg);
%                   accSigy1(:,iMix)=accSigy1(:,iMix)+postProb(iMix)*trainData(l/2+1:l,start+iSeg);
%          for iO=1:order
%                 accSigx2(:,iO,iMix)=accSigx2(:,iO,iMix)+postProb(iMix)*trainData(1:l/2,start+iSeg).*preSegs(1:l/2,iO);
%                 accSigy2(:,iO,iMix)=accSigy2(:,iO,iMix)+postProb(iMix)*trainData(l/2+1:l,start+iSeg).*preSegs(l/2+1:l,iO);
%          end
%                 
%          accSigy3(:,iMix)=accSigy3(:,iMix)+postProb(iMix).*trainData(l/2+1:l,start+iSeg).*trainData(1:l/2,start+iSeg);
   
                  accSig(:,iMix)=accSig(:,iMix)+postProb(iMix)*trainData(1:l,start+iSeg).*trainData(1:l,start+iSeg);
               
                   accSig1(:,iMix)=accSig1(:,iMix)+postProb(iMix)*trainData(1:l,start+iSeg);
                   
               
        
                
          
         
                
         accSigy3(:,iMix)=accSigy3(:,iMix)+postProb(iMix).*trainData(l/2+1:l,start+iSeg).*trainData(1:l/2,start+iSeg);
         
         
        end
        
         
        mean(:)=mean(:)+trainData(:,start+iSeg);  %for Sigma update
        numOfSeg=numOfSeg+1;
    end
    
   
    
    
end
preProb=curProb;

 curProb=liklihoodProb/numOfSeg;
mean=mean/numOfSeg;
 fprintf(1,'curProb %f\n',curProb);
  for iMix=1:mixNum
        for iL=1:l/2
              tempA=squeeze(A(iL,:,:,iMix));
            if min(svd(tempA))<1e-7
                temp=pinv(tempA)* C(iL,:,iMix)';
            else
                temp=tempA \ C(iL,:,iMix)';
            end
           
            arx(iL,:,iMix)=temp(1:orderx);
            miux(iL,iMix)=temp(orderx+1);
        end
        for iL=1:l/2
            tempB=squeeze( B(iL,:,:,iMix) );
             if min(svd(tempB))<1e-7
                temp=pinv(tempB)* D(iL,:,iMix)';
            else
                temp=tempB \ D(iL,:,iMix)';
             end
            ary(iL,:,iMix)=temp(1:ordery);
            miuy(iL,iMix)=temp(ordery+1);
            bt(iL,iMix)=temp(ordery+2);
           % bt(iL,iMix)=0;
        end
        
       w(iMix)=accW(iMix)/numOfSeg;
         sigmaxx(:,:,iMix)=diag(accSig(1:l/2,iMix)-accSig1(1:l/2,iMix).*miux(:,iMix))./accW(iMix);
         
          sigmayy(:,:,iMix)=diag(accSig(l/2+1:l,iMix)-accSig1(l/2+1:l,iMix).*miuy(:,iMix)  -accSigy3(:,iMix).*bt(:,iMix)        )./accW(iMix);
          
     for iO=1:ordery  
           
          %  sigmaxx(:,:,iMix)=sigmaxx(:,:,iMix)-diag(accSig2(1:l/2,iO,iMix).*arx(:,iO,iMix) )./accW(iMix);
            sigmayy(:,:,iMix)=sigmayy(:,:,iMix)-diag( accSig2(l/2+1:l,iO,iMix).*ary(:,iO,iMix)  )./accW(iMix);

     end
     sigma(:,:,iMix)=blkdiag(sigmaxx(:,:,iMix),sigmayy(:,:,iMix));
     
         [D_,p_]=chol(sigma(:,:,iMix));
      if p_>0
        sigma(:,:,iMix)=sigma(:,:,iMix)+ diag(vFloor);
      end
        gconst(iMix)=log(w(iMix))-l/2*log(2*pi)-0.5*log(det(sigma(:,:,iMix)));
    invSigma(:,:,iMix)=inv( sigma(:,:,iMix) );
   end
    %
end

