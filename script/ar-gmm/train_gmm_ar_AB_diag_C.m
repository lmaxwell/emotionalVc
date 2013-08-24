clear;
emotion='angry';
pMix=6;
sMix=12;
trainDir=['../../train/tone_nodynamic_diag/' emotion '/' num2str(pMix) 'mix' num2str(sMix) 'mix_3'];
gmmFile=[trainDir '/neutral_' emotion '.tone.gmm'];
trainFile=[trainDir '/neutral_' emotion '.tone_ar'];
l=16;
mixNum=12;
order=1;
vfloor=0.001;
%% read gmm model
 gmmP=importdata(gmmFile);

 
   w=gmmP(1:mixNum);
    miux=zeros(l/2,mixNum);
    miuy=zeros(l/2,mixNum);
       miu=zeros(l,mixNum);
    sigma=zeros(l,l,mixNum);
       beginModel.w=w;
    for i=1:mixNum
        miux(:,i)=gmmP(mixNum+(i-1)*2*l+1:mixNum+(i-1)*2*l+l/2);
        miuy(:,i)=gmmP(mixNum+(i-1)*2*l+l/2+1:mixNum+(i-1)*2*l+l);
         miu(:,i)=[miux(:,i);miuy(:,i)];
        sigma(:,:,i)=diag(gmmP(mixNum+(i-1)*2*l+l+1:mixNum+(i-1)*2*l+2*l));
         beginModel.miu(:,i)=[miux(:,i);miuy(:,i)];
       
        beginModel.sigma(:,:,i)=sigma(:,:,i);
    end
    
    sigmaxx=sigma(1:l/2,1:l/2,:);
    sigmayy=sigma(l/2+1:l,l/2+1:l,:);
    sigmaxy=sigma(1:l/2,l/2+1:l,:);
    sigmayx=sigma(l/2+1:l,1:l/2,:);
 

%    
    % set ar coefficient
    ar=zeros(l,order,mixNum);
   arx=zeros(l/2,order,mixNum);
ary=zeros(l/2,order,mixNum);
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

%%
indexOfSen=find(trainData(1,:)==0);
numOfSen=length(indexOfSen);
indexOfSen=[indexOfSen  size(trainData,2)+1];
beginZero=zeros(l,1);

A=zeros(l/2,order+1,order+1,mixNum);
C=zeros(l/2,order+1,mixNum);
B=zeros(l/2,order+2,order+2,mixNum);
D=zeros(l/2,order+2,mixNum);
mean=zeros(l,1);

postProb=zeros(mixNum,1);

liklihoodProb=0;
numOfSeg=0;
accW=zeros(mixNum,1);
accSigx=zeros(l/2,mixNum);
accSigx1=zeros(l/2,mixNum);
accSigx2=zeros(l/2,order,mixNum);
accSigy=zeros(l/2,mixNum);
accSigy1=zeros(l/2,mixNum);
accSigy2=zeros(l/2,order,mixNum);
accSigy3=zeros(l/2,mixNum);
vFloor=zeros(l,1);

for iSen=1:numOfSen
    start=indexOfSen(iSen)+1;
    xend=indexOfSen(iSen+1)-1;
    fprintf( 1 , '%dth sentence,start:%d,end:%d\n' , iSen , start , xend);
    for iSeg=0:xend-start
        fprintf(1,'%dth segment\n',iSeg);
        %set preSegs
        if iSeg==0
              preSegs=beginZero*ones(1,order-iSeg);
        elseif(iSeg<order)
           
            preSegs=[beginZero*ones(1,order-iSeg) trainData(:,start+iSeg-1:-1:start)];
           
        else
            preSegs=trainData(:,start+iSeg-1:-1:start+iSeg-order);
       end
        
        sum_postProb=-1.0E10;
        maxProb=-1.0E10;
        maxMix=1;
        for iMix=1:mixNum
            
            postProb(iMix)= log(w(iMix))-l/2*log(2*pi)-0.5*log(det(sigma(:,:,iMix)))-0.5.*(trainData(:,start+iSeg)-miu(:,iMix))'/sigma(:,:,iMix)*(trainData(:,start+iSeg)-miu(:,iMix))          ;
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
           
            for iOd=1:order  
                for jOd=1:order 
                A(  :,iOd,jOd , iMix)= A(  :,iOd,jOd ,iMix)+ postProb(iMix)*preSegs(1:l/2,iOd).*preSegs(1:l/2,jOd);  %iOd th function , jOd th parameter
                
                 B(  :,iOd,jOd , iMix)= B(  :,iOd,jOd ,iMix)+ postProb(iMix)*preSegs(l/2+1:l,iOd).*preSegs(l/2+1:l,jOd);  %iOd th function , jOd th parameter
                end
                C(:,iOd,iMix)=C(:,iOd,iMix)+postProb(iMix)*trainData(1:l/2,start+iSeg).*preSegs(1:l/2,iOd);  %iOd th function
                A(:,order+1,iOd,iMix)=A(:,order+1,iOd,iMix)+postProb(iMix).*preSegs(1:l/2,iOd);  % order+1 th fuction, iOd th parameter
                A(:,iOd,order+1,iMix)=A(:,iOd,order+1,iMix)+ postProb(iMix).*preSegs(1:l/2,iOd);         % iOd th function , order+1 th parameter
                
                D(:,iOd,iMix)=D(:,iOd,iMix)+postProb(iMix)*trainData(l/2+1:l,start+iSeg).*preSegs(l/2+1:l,iOd);  %iOd th function
                B(:,order+1,iOd,iMix)=B(:,order+1,iOd,iMix)+postProb(iMix).*preSegs(l/2+1:l,iOd);  % order+1 th fuction, iOd th parameter
                B(:,order+2,iOd,iMix)=B(:,order+2,iOd,iMix)+postProb(iMix).*trainData(1:l/2,start+iSeg).*preSegs(l/2+1:l,iOd);
                B(:,iOd,order+1,iMix)=B(:,iOd,order+1,iMix)+ postProb(iMix).*preSegs(l/2+1:l,iOd);         % iOd th function , order+1 th parameter
               B(:,iOd,order+2,iMix)=B(:,iOd,order+2,iMix)+ postProb(iMix).*trainData(1:l/2,start+iSeg).*preSegs(l/2+1:l,iOd);  
            end
               A(:,order+1,order+1,iMix)=   A(:,order+1,order+1,iMix)+postProb(iMix);   %order+1 th function ,order+1 th parameter
               C(:,order+1,iMix)=C(:,order+1,iMix)+trainData(1:l/2,start+iSeg).*postProb(iMix); % order+1th function
                
              B(:,order+1,order+1,iMix)=   B(:,order+1,order+1,iMix)+postProb(iMix);   %order+1 th function ,order+1 th parameter
               D(:,order+1,iMix)=D(:,order+1,iMix)+trainData(l/2+1:l,start+iSeg).*postProb(iMix); % order+1th function
                B(:,order+1,order+2,iMix)=   B(:,order+1,order+2,iMix)+postProb(iMix).*trainData(1:l/2,start+iSeg);  
                  B(:,order+2,order+1,iMix)=   B(:,order+2,order+1,iMix)+postProb(iMix).*trainData(1:l/2,start+iSeg);
                    B(:,order+2,order+2,iMix)=   B(:,order+2,order+2,iMix)+postProb(iMix).*trainData(1:l/2,start+iSeg).*trainData(1:l/2,start+iSeg);
               D(:,order+2,iMix)=D(:,order+2,iMix)+trainData(l/2+1:l,start+iSeg).*trainData(1:l/2,start+iSeg).*postProb(iMix); 
               
               
                    accSigx(:,iMix)=accSigx(:,iMix)+postProb(iMix)*trainData(1:l/2,start+iSeg).*trainData(1:l/2,start+iSeg);
               
                   accSigx1(:,iMix)=accSigx1(:,iMix)+postProb(iMix)*trainData(1:l/2,start+iSeg);
                   
                  accSigy(:,iMix)=accSigy(:,iMix)+postProb(iMix)*trainData(l/2+1:l,start+iSeg).*trainData(l/2+1:l,start+iSeg);
                  accSigy1(:,iMix)=accSigy1(:,iMix)+postProb(iMix)*trainData(l/2+1:l,start+iSeg);
         for iO=1:order
                accSigx2(:,iO,iMix)=accSigx2(:,iO,iMix)+postProb(iMix)*trainData(1:l/2,start+iSeg).*preSegs(1:l/2,iO);
                accSigy2(:,iO,iMix)=accSigy2(:,iO,iMix)+postProb(iMix)*trainData(l/2+1:l,start+iSeg).*preSegs(l/2+1:l,iO);
         end
                
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
           
            arx(iL,:,iMix)=temp(1:order);
            miux(iL,iMix)=temp(order+1);
        end
        for iL=1:l/2
            tempB=squeeze( B(iL,:,:,iMix) );
             if min(svd(tempB))<1e-7
                temp=pinv(tempB)* D(iL,:,iMix)';
            else
                temp=tempB \ D(iL,:,iMix)';
             end
            ary(iL,:,iMix)=temp(1:order);
            miuy(iL,iMix)=temp(order+1);
            bt(iL,iMix)=temp(order+2);
               %bt(iL,iMix)=0;
        end
        
       w(iMix)=accW(iMix)/numOfSeg;
         sigmaxx(:,:,iMix)=diag(accSigx(:,iMix)-accSigx1(:,iMix).*miux(:,iMix))./accW(iMix);
         
          sigmayy(:,:,iMix)=diag(accSigy(:,iMix)-accSigy1(:,iMix).*miuy(:,iMix)  -accSigy3(:,iMix).*bt(:,iMix)        )./accW(iMix);
          
     for iO=1:order  
           
            sigmaxx(:,:,iMix)=sigmaxx(:,:,iMix)-diag(accSigx2(:,iO,iMix).*arx(:,iO,iMix) )./accW(iMix);
            sigmayy(:,:,iMix)=sigmayy(:,:,iMix)-diag( accSigy2(:,iO,iMix).*ary(:,iO,iMix)  )./accW(iMix);

     end
     sigma(:,:,iMix)=blkdiag(sigmaxx(:,:,iMix),sigmayy(:,:,iMix));
      [D_,p_]=chol(sigma(:,:,iMix));
      if p_>0
        sigma(:,:,iMix)=sigma(:,:,iMix)+ diag(vFloor);
      end
      
      
    end

ogmm=fopen([trainDir '/gmm_ar_diag_1thit_order' num2str(order) ],'w');
for iMix=1:mixNum
    fprintf(ogmm,'%f\n',w(iMix)); 
end
for iMix=1:mixNum
    fprintf(ogmm,'%f\n',miu(:,iMix));
    for iO=1:order
        fprintf(ogmm,'%f\n',ar(:,order,iMix));
    end
    fprintf(ogmm,'%f\n',sigma(:,:,iMix));
end
%%
preProb=liklihoodProb/numOfSeg;
curProb=0;
 fprintf(1,'initialProb %f\n',preProb);
%iteration
ite=1;
while(abs((curProb-preProb)/preProb)>1e-5 && ite<20)
A=zeros(l/2,order+1,order+1,mixNum);
C=zeros(l/2,order+1,mixNum);
B=zeros(l/2,order+2,order+2,mixNum);
D=zeros(l/2,order+2,mixNum);
mean=zeros(l,1);

postProb=zeros(mixNum,1);

liklihoodProb=0;
numOfSeg=0;
accW=zeros(mixNum,1);
accSigx=zeros(l/2,mixNum);
accSigx1=zeros(l/2,mixNum);
accSigx2=zeros(l/2,order,mixNum);
accSigy=zeros(l/2,mixNum);
accSigy1=zeros(l/2,mixNum);
accSigy2=zeros(l/2,order,mixNum);
accSigy3=zeros(l/2,mixNum);

for iSen=1:numOfSen
    start=indexOfSen(iSen)+1;
    xend=indexOfSen(iSen+1)-1;
  %  fprintf( 1 , '%dth sentence,start:%d,end:%d\n' , iSen , start , xend);
    for iSeg=0:xend-start
 %       fprintf(1,'%dth segment\n',iSeg);
        %set preSegs
         if iSeg==0
              preSegs=beginZero*ones(1,order-iSeg);
        elseif(iSeg<order)
           
            preSegs=[beginZero*ones(1,order-iSeg) trainData(:,start+iSeg-1:-1:start)];
           
        else
            preSegs=trainData(:,start+iSeg-1:-1:start+iSeg-order);
        end
        
        sum_postProb=-1.0E10;
       
        for iMix=1:mixNum

            if iSeg<0
                curMiu=beginModel.miu(:,iMix);
                
                 postProb(iMix)= w(iMix).*1./((2*pi)^(l/2).*sqrt(det(beginModel.sigma(:,:,iMix)))).* ...
            exp(-0.5.*(trainData(:,start+iSeg)-curMiu)'/beginModel.sigma(:,:,iMix)*(trainData(:,start+iSeg)-curMiu));
            else
                curMiu=[miux(:,iMix);miuy(:,iMix)];
                for iO=1:order
                    curMiu(1:l/2)=curMiu(1:l/2)+arx(:,iO,iMix).*preSegs(1:l/2,iO);
                    curMiu(l/2+1:l)=curMiu(l/2+1:l)+ary(:,iO,iMix).*preSegs(l/2+1:l,iO);
                  
                end
                   curMiu(l/2+1:l)=curMiu(l/2+1:l)+bt(:,iMix).*trainData(1:l/2,start+iSeg);
                
                 postProb(iMix)= log(w(iMix))-l/2*log(2*pi)-0.5*log(det(sigma(:,:,iMix)))-0.5.*(trainData(:,start+iSeg)-curMiu)'/sigma(:,:,iMix)*(trainData(:,start+iSeg)-curMiu)          ;
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
           
            for iOd=1:order  
                for jOd=1:order 
                A(  :,iOd,jOd , iMix)= A(  :,iOd,jOd ,iMix)+ postProb(iMix)*preSegs(1:l/2,iOd).*preSegs(1:l/2,jOd);  %iOd th function , jOd th parameter
                
                 B(  :,iOd,jOd , iMix)= B(  :,iOd,jOd ,iMix)+ postProb(iMix)*preSegs(l/2+1:l,iOd).*preSegs(l/2+1:l,jOd);  %iOd th function , jOd th parameter
                end
                C(:,iOd,iMix)=C(:,iOd,iMix)+postProb(iMix)*trainData(1:l/2,start+iSeg).*preSegs(1:l/2,iOd);  %iOd th function
                A(:,order+1,iOd,iMix)=A(:,order+1,iOd,iMix)+postProb(iMix).*preSegs(1:l/2,iOd);  % order+1 th fuction, iOd th parameter
                A(:,iOd,order+1,iMix)=A(:,iOd,order+1,iMix)+ postProb(iMix).*preSegs(1:l/2,iOd);         % iOd th function , order+1 th parameter
                
                D(:,iOd,iMix)=D(:,iOd,iMix)+postProb(iMix)*trainData(l/2+1:l,start+iSeg).*preSegs(l/2+1:l,iOd);  %iOd th function
                B(:,order+1,iOd,iMix)=B(:,order+1,iOd,iMix)+postProb(iMix).*preSegs(l/2+1:l,iOd);  % order+1 th fuction, iOd th parameter
                B(:,order+2,iOd,iMix)=B(:,order+2,iOd,iMix)+postProb(iMix).*trainData(1:l/2,start+iSeg).*preSegs(l/2+1:l,iOd);
                B(:,iOd,order+1,iMix)=B(:,iOd,order+1,iMix)+ postProb(iMix).*preSegs(l/2+1:l,iOd);         % iOd th function , order+1 th parameter
               B(:,iOd,order+2,iMix)=B(:,iOd,order+2,iMix)+ postProb(iMix).*trainData(1:l/2,start+iSeg).*preSegs(l/2+1:l,iOd);  
            end
               A(:,order+1,order+1,iMix)=   A(:,order+1,order+1,iMix)+postProb(iMix);   %order+1 th function ,order+1 th parameter
               C(:,order+1,iMix)=C(:,order+1,iMix)+trainData(1:l/2,start+iSeg).*postProb(iMix); % order+1th function
                
              B(:,order+1,order+1,iMix)=   B(:,order+1,order+1,iMix)+postProb(iMix);   %order+1 th function ,order+1 th parameter
               D(:,order+1,iMix)=D(:,order+1,iMix)+trainData(l/2+1:l,start+iSeg).*postProb(iMix); % order+1th function
              B(:,order+1,order+2,iMix)=   B(:,order+1,order+2,iMix)+postProb(iMix).*trainData(1:l/2,start+iSeg);  
                  B(:,order+2,order+1,iMix)=   B(:,order+2,order+1,iMix)+postProb(iMix).*trainData(1:l/2,start+iSeg);
                    B(:,order+2,order+2,iMix)=   B(:,order+2,order+2,iMix)+postProb(iMix).*trainData(1:l/2,start+iSeg).*trainData(1:l/2,start+iSeg);
               D(:,order+2,iMix)=D(:,order+2,iMix)+trainData(l/2+1:l,start+iSeg).*trainData(1:l/2,start+iSeg).*postProb(iMix); 
               
               
                    accSigx(:,iMix)=accSigx(:,iMix)+postProb(iMix)*trainData(1:l/2,start+iSeg).*trainData(1:l/2,start+iSeg);
               
                   accSigx1(:,iMix)=accSigx1(:,iMix)+postProb(iMix)*trainData(1:l/2,start+iSeg);
                   
                  accSigy(:,iMix)=accSigy(:,iMix)+postProb(iMix)*trainData(l/2+1:l,start+iSeg).*trainData(l/2+1:l,start+iSeg);
                  accSigy1(:,iMix)=accSigy1(:,iMix)+postProb(iMix)*trainData(l/2+1:l,start+iSeg);
         for iO=1:order
                accSigx2(:,iO,iMix)=accSigx2(:,iO,iMix)+postProb(iMix)*trainData(1:l/2,start+iSeg).*preSegs(1:l/2,iO);
                accSigy2(:,iO,iMix)=accSigy2(:,iO,iMix)+postProb(iMix)*trainData(l/2+1:l,start+iSeg).*preSegs(l/2+1:l,iO);
         end
                
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
           
            arx(iL,:,iMix)=temp(1:order);
            miux(iL,iMix)=temp(order+1);
        end
        for iL=1:l/2
            tempB=squeeze( B(iL,:,:,iMix) );
             if min(svd(tempB))<1e-7
                temp=pinv(tempB)* D(iL,:,iMix)';
            else
                temp=tempB \ D(iL,:,iMix)';
             end
            ary(iL,:,iMix)=temp(1:order);
            miuy(iL,iMix)=temp(order+1);
            bt(iL,iMix)=temp(order+2);
           % bt(iL,iMix)=0;
        end
        
       w(iMix)=accW(iMix)/numOfSeg;
         sigmaxx(:,:,iMix)=diag(accSigx(:,iMix)-accSigx1(:,iMix).*miux(:,iMix))./accW(iMix);
         
          sigmayy(:,:,iMix)=diag(accSigy(:,iMix)-accSigy1(:,iMix).*miuy(:,iMix)  -accSigy3(:,iMix).*bt(:,iMix)        )./accW(iMix);
          
     for iO=1:order  
           
            sigmaxx(:,:,iMix)=sigmaxx(:,:,iMix)-diag(accSigx2(:,iO,iMix).*arx(:,iO,iMix))./accW(iMix);
            sigmayy(:,:,iMix)=sigmayy(:,:,iMix)-diag( accSigy2(:,iO,iMix).*ary(:,iO,iMix)  )./accW(iMix);

     end
     sigma(:,:,iMix)=blkdiag(sigmaxx(:,:,iMix),sigmayy(:,:,iMix));
         [D_,p_]=chol(sigma(:,:,iMix));
      if p_>0
           fprintf(1,'add variance floor:%d\n',nFloor);
        sigma(:,:,iMix)=sigma(:,:,iMix)+ diag(vFloor);
      end
   end
    %
    ite=ite+1;
end

%%
ogmm=fopen([trainDir '/neutral_' emotion '.tone.gmm_ar_ABdiag'],'w');
for iMix=1:mixNum
    fprintf(ogmm,'%f\n',w(iMix));
end
for iMix=1:mixNum
    fprintf(ogmm,'%f\n',miux(:,iMix));
      fprintf(ogmm,'%f\n',miuy(:,iMix));
    for iO=1:order
        fprintf(ogmm,'%f\n',arx(:,iO,iMix));
      
    end
    for iO=1:order
        fprintf(ogmm,'%f\n',ary(:,iO,iMix));
    end
    
    fprintf(ogmm,'%f\n',bt(:,iMix));
    fprintf(ogmm,'%f\n',sigma(:,:,iMix));
end
