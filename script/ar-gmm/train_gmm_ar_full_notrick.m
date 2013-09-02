clear;
pMix=6;
sMix=12;
emotion='sad';
trainDir=['../../train/tone_nodynamic/' emotion '/' num2str(pMix) 'mix'  num2str(sMix) 'mix_3'];
gmmFile=[trainDir '/neutral_' emotion '.tone.gmm'];
trainFile=[trainDir '/neutral_' emotion '.tone_ar'];
l=16;
mixNum=sMix;
%mixNum=3;
order=0  ;
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
        miux(:,i)=gmmP(mixNum+(i-1)*(l^2+l)+1:mixNum+(i-1)*(l^2+l)+l/2);
        miuy(:,i)=gmmP(mixNum+(i-1)*(l^2+l)+l/2+1:mixNum+(i-1)*(l^2+l)+l);
        miu(:,i)=[miux(:,i);miuy(:,i)];
        beginModel.miu(:,i)=[miux(:,i);miuy(:,i)];
        sigma(:,:,i)=reshape(gmmP(mixNum+(i-1)*(l^2+l)+l+1:mixNum+(i-1)*(l^2+l)+l+l^2),l,l);
        beginModel.sigma(:,:,i)=sigma(:,:,i);
    end

    sigmaxx=sigma(1:l/2,1:l/2,:);
    sigmayy=sigma(l/2+1:l,l/2+1:l,:);
    sigmaxy=sigma(1:l/2,l/2+1:l,:);
    sigmayx=sigma(l/2+1:l,1:l/2,:);
   
    % set ar coefficient
    ar=zeros(l,order,mixNum);
   

    
    
%% train ar combined gmm

%first run  
trainData=importdata(trainFile)';
%
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
vFloor=zeros(l,1);
ite=1;
curProb=-1E10;
preProb=0;
while(ite<29 && abs(curProb-preProb)>1e-4)

A=zeros(order*l+1,order*l+1,l,mixNum);
B=zeros(order*l+1,l,mixNum);
mean=zeros(l,1);

postProb=zeros(mixNum,1);

liklihoodProb=0;
numOfSeg=0;
accW=zeros(mixNum,1);
accSig=zeros(l,l,mixNum);
accSig1=zeros(l,mixNum);
accSig2=zeros(l,l,order,mixNum);
accSig3=zeros(l,l,order+1,order+1,mixNum);



for iSen=1:numOfSen
    start=indexOfSen(iSen)+1;
    xend=indexOfSen(iSen+1)-1;
   % fprintf( 1 , '%dth sentence,start:%d,end:%d\n' , iSen , start , xend);
    for iSeg=0:xend-start
   %    fprintf(1,'%dth segment\n',iSeg);
     %   %set preSegs
        if iSeg==0
              preSegs=beginZero*ones(1,order-iSeg);
        elseif(iSeg<order)
           
            preSegs=[ trainData(:,start+iSeg-1:-1:start) beginZero*ones(1,order-iSeg)];
           
        else
            preSegs=trainData(:,start+iSeg-1:-1:start+iSeg-order);
       end
        
      sum_postProb=-1.0E10;
        maxProb=-1.0E10;
        maxMix=1;
        for iMix=1:mixNum
            curMiu=miu(:,iMix);
                for iO=1:order
                    curMiu=curMiu+ar(:,iO,iMix).*preSegs(:,iO);
                end
               postProb(iMix)= log(w(iMix))-l/2*log(2*pi)-0.5*log(det(sigma(:,:,iMix)))-0.5.*(trainData(:,start+iSeg)-curMiu)'/sigma(:,:,iMix)*(trainData(:,start+iSeg)-curMiu) ; 
           % postProb(iMix)= log(w(iMix))-l/2*log(2*pi)-0.5*log(det(sigma(:,:,iMix)))-0.5.*(trainData(:,start+iSeg)-miu(:,iMix))'/sigma(:,:,iMix)*(trainData(:,start+iSeg)-miu(:,iMix))          ;
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
          for iL=1:l
        
            s=[ reshape(preSegs(:,:),1,l*order) 1  ]';
          
            for iMix=1:mixNum
                A(:,:,iL,iMix )= A(:,:,iL,iMix) + postProb(iMix).*s*s';
                 B(:,iL,iMix)=B(:,iL,iMix) + postProb(iMix).*s*trainData(iL,start+iSeg);
            end
            
       end
        
        
        
        for iMix=1:mixNum 
           
            for iOd=1:order  
                for jOd=1:order 
          %      A(  :,iOd,jOd , iMix)= A(  :,iOd,jOd ,iMix)+ postProb(iMix)*preSegs(:,iOd).*preSegs(:,jOd);  %iOd th function , jOd th parameter
                
                   accSig3(:,:,iOd,jOd,iMix)=accSig3(:,:,iOd,jOd,iMix) + postProb(iMix) .* preSegs(:,iOd) *preSegs(:,jOd)';
                end
            %    B(:,iOd,iMix)=B(:,iOd,iMix)+postProb(iMix)*trainData(:,start+iSeg).*preSegs(:,iOd);  %iOd th function
              %  A(:,order+1,iOd,iMix)=A(:,order+1,iOd,iMix)+postProb(iMix).*preSegs(:,iOd);  % order+1 th fuction, iOd th parameter
                %A(:,iOd,order+1,iMix)=A(:,iOd,order+1,iMix)+ postProb(iMix).*preSegs(:,iOd);         % iOd th function , order+1 th parameter
                
             
                accSig3(:,:,iOd,order+1)=accSig3(:,:,iOd,order+1) + postProb(iMix) .* preSegs(:,iOd)*ones(1,l);
               accSig3(:,:,order+1,iOd)=accSig3(:,:,order+1,iOd) + postProb(iMix) .* ones(l,1)*preSegs(:,iOd)';
            end
               %A(:,order+1,order+1,iMix)=   A(:,order+1,order+1,iMix)+postProb(iMix);   %order+1 th function ,order+1 th parameter
               %B(:,order+1,iMix)=B(:,order+1,iMix)+trainData(:,start+iSeg).*postProb(iMix); % order+1th function
                
               accSig(:,:,iMix)=accSig(:,:,iMix)+postProb(iMix)*trainData(:,start+iSeg)*trainData(:,start+iSeg)';
               
          accSig3(:,:,order+1,order+1)=accSig3(:,:,order+1,order+1) + postProb(iMix) ;
                 
                
                       
         accSig1(:,iMix)=accSig1(:,iMix)+postProb(iMix)*trainData(:,start+iSeg);
              
                       
                       for iO=1:order
                            accSig2(:,:,iO,iMix)=accSig2(:,:,iO,iMix)+postProb(iMix)*trainData(:,start+iSeg)*preSegs(:,iO)';
                       end
            
%                 accSig1(:,:,iMix)=accSig1(:,:,iMix)+postProb(iMix)*preSegs;
%          for iO=1:order
%              for jO=1:order
%                 accSig2(:,:,iO,jO,iMix)=accSig2(:,:,iO,jO,iMix)+postProb(iMix)*preSegs(:,iO)*preSegs(:,jO)';
%              end
%          end
               
        end
         mean(:)=mean(:)+trainData(:,start+iSeg);  %for Sigma update
        if (ite==1)
        vFloor=vFloor+(trainData(:,start+iSeg)-miu(:,maxMix)).* (trainData(:,start+iSeg)-miu(:,maxMix));
        end
        numOfSeg=numOfSeg+1;
        
    end
    mean=mean/numOfSeg;
    

    
    
end
preProb=curProb;
curProb=liklihoodProb/numOfSeg;
 fprintf(1,'curProb %f\n',curProb);

if(ite==1)
vFloor=vFloor*vfloor/numOfSeg;
end
%%
artemp=zeros(l,l,order);
    for iMix=1:mixNum
        for iL=1:l
              tempA=squeeze(A(:,:,iL,iMix));
            if min(svd(tempA))<1e-4
                temp=pinv(tempA)* B(:,iL,iMix);
          %      fprintf(1,'%dmix,%ddim,pinv\n',iMix,iL);
            else
                temp=tempA \ B(:,iL,iMix);
            end
            
           temptemp=reshape(temp(1:order*l),l,order)';
           for iO=1:order
            artemp(iL,:,iO)=temptemp(iO,:);
             ar(iL,iO,iMix)=artemp(iL,iL,iO);
           end
           
            miu(iL,iMix)=temp(order*l+1);
        end
       w(iMix)=accW(iMix)/numOfSeg;
       
       sigma(:,:,iMix)=accSig(:,:,iMix)./accW(iMix) - accSig1(:,iMix)*miu(:,iMix)'./accW(iMix);
       for iO=1:order
        sigma(:,:,iMix) = sigma(:,:,iMix) - accSig2(:,:,iO,iMix)*artemp(:,:,iO)'./accW(iMix);
       end
%       sigma(:,:,iMix)=( accSig(:,:,iMix)-miu(:,iMix)*accSig1(:,iMix)'-accSig1(:,iMix)*miu(:,iMix)' )/accW(iMix) +miu(:,iMix)*miu(:,iMix)';
%          for iO=1:order
%            sigma(:,:,iMix)=sigma(:,:,iMix)-accSig2(:,:,iO,iMix)*artemp(:,:,iO)'./accW(iMix) - artemp(:,:,iO)*accSig2(:,:,iO,iMix)'./accW(iMix);
%           sigma(:,:,iMix)=sigma(:,:,iMix)+ diag(miu(:,iMix))*accSig3(:,:,order+1,iO)*artemp(:,:,iO)'./accW(iMix);
%               sigma(:,:,iMix)=sigma(:,:,iMix)+ artemp(:,:,iO)*accSig3(:,:,iO,order+1)*diag(miu(:,iMix))./accW(iMix);
%            for jO=1:order
%                sigma(:,:,iMix)=sigma(:,:,iMix) +  artemp(:,:,iO) * accSig3(:,:,iO,jO,iMix) *artemp(:,:,jO)'./accW(iMix);
%            end
%          end
       varx=diag(sigma(1:l/2,1:l/2,iMix));
       vary=diag(sigma(l/2+1:l,l/2+1:l,iMix));
       varxy=diag(sigma(1:l/2,l/2+1:l,iMix));
       sigma(:,:,iMix)=[diag(varx) diag(varxy);diag(varxy) diag(vary) ];
%      sigma(:,:,iMix)=accSig(:,:,iMix)/accW(iMix)-miu(:,iMix)*miu(:,iMix)';
%      for iO=1:order
%             tempSig1=miu(:,iMix)*accSig1(:,iO,iMix)'*diag(ar(:,iO,iMix))/accW(iMix);
%             sigma(:,:,iMix)=sigma(:,:,iMix)-tempSig1-tempSig1';
%             for jO=1:order
%                  sigma(:,:,iMix)=sigma(:,:,iMix)-diag(ar(:,iO,iMix))*accSig2(:,:,iO,jO,iMix)*diag(ar(:,jO,iMix))'/accW(iMix);
%             end
%      end
    nFloor=1;
      [D_,dp]=chol(sigma(:,:,iMix));
      if dp>0
         fprintf(1,'add variance floor:%d\n',nFloor);
        sigma(:,:,iMix)= diag(vFloor)+sigma(:,:,iMix);
        [D_,dp]=chol(sigma(:,:,iMix));
        nFloor=nFloor+1;
      end
    end
end


%%
preProb=liklihoodProb/numOfSeg;
curProb=0;
 fprintf(1,'initialProb %f\n',preProb);
%iteration
while(abs((curProb-preProb)/preProb)>1e-5)
liklihoodProb=0;
 numOfSeg=0;
 accW=zeros(mixNum,1);
 A=zeros(l,order+1,order+1,mixNum);
B=zeros(l,order+1,mixNum);
mean=zeros(l,1);
accSig=zeros(l,l,mixNum);
accSig1=zeros(l,mixNum);
accSig2=zeros(l,l,order,mixNum);
accSig3=zeros(l,l,order+1,order+1,mixNum);

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
           
                  preSegs=[ trainData(:,start+iSeg-1:-1:start) beginZero*ones(1,order-iSeg)];
           
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
                curMiu=miu(:,iMix);
                for iO=1:order
                    curMiu=curMiu+ar(:,iO,iMix).*preSegs(:,iO);
                end
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
                A(  :,iOd,jOd , iMix)= A(  :,iOd,jOd ,iMix)+ postProb(iMix)*preSegs(:,iOd).*preSegs(:,jOd);  %iOd th function , jOd th parameter
                
                   accSig3(:,:,iOd,jOd,iMix)=accSig3(:,:,iOd,jOd,iMix) + postProb(iMix) .* preSegs(:,iOd) *preSegs(:,jOd)';
                end
                B(:,iOd,iMix)=B(:,iOd,iMix)+postProb(iMix)*trainData(:,start+iSeg).*preSegs(:,iOd);  %iOd th function
                A(:,order+1,iOd,iMix)=A(:,order+1,iOd,iMix)+postProb(iMix).*preSegs(:,iOd);  % order+1 th fuction, iOd th parameter
                A(:,iOd,order+1,iMix)=A(:,iOd,order+1,iMix)+ postProb(iMix).*preSegs(:,iOd);         % iOd th function , order+1 th parameter
                
             
                accSig3(:,:,iOd,order+1)=accSig3(:,:,iOd,order+1) + postProb(iMix) .* preSegs(:,iOd)*ones(1,l);
               accSig3(:,:,order+1,iOd)=accSig3(:,:,order+1,iOd) + postProb(iMix) .* ones(l,1)*preSegs(:,iOd)';
            end
               A(:,order+1,order+1,iMix)=   A(:,order+1,order+1,iMix)+postProb(iMix);   %order+1 th function ,order+1 th parameter
               B(:,order+1,iMix)=B(:,order+1,iMix)+trainData(:,start+iSeg).*postProb(iMix); % order+1th function
                
               accSig(:,:,iMix)=accSig(:,:,iMix)+postProb(iMix)*trainData(:,start+iSeg)*trainData(:,start+iSeg)';
               
          accSig3(:,:,order+1,order+1)=accSig3(:,:,order+1,order+1) + postProb(iMix) ;
                 
                
                       
         accSig1(:,iMix)=accSig1(:,iMix)+postProb(iMix)*trainData(:,start+iSeg);
              
                       
                       for iO=1:order
                            accSig2(:,:,iO,iMix)=accSig2(:,:,iO,iMix)+postProb(iMix)*trainData(:,start+iSeg)*preSegs(:,iO)';
                       end
%                
%                    accSig1(:,:,iMix)=accSig1(:,:,iMix)+postProb(iMix)*preSegs(:,:);
%          for iO=1:order
%              for jO=1:order
%                 accSig2(:,:,iO,jO,iMix)=accSig2(:,:,iO,jO,iMix)+postProb(iMix)*preSegs(:,iO)*preSegs(:,jO)';
%              end
%          end
%                
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
        for iL=1:l
            tempA=squeeze(A(iL,:,:,iMix));
            if min(svd(tempA))<1e-5
                temp=pinv(tempA)* B(iL,:,iMix)';
            else
                temp=tempA \ B(iL,:,iMix)';
            end
            
            ar(iL,:,iMix)=temp(1:order);
            miu(iL,iMix)=temp(order+1);
        end
         w(iMix)=accW(iMix)/numOfSeg;
     sigma(:,:,iMix)=( accSig(:,:,iMix)-miu(:,iMix)*accSig1(:,iMix)'-accSig1(:,iMix)*miu(:,iMix)' )/accW(iMix) +miu(:,iMix)*miu(:,iMix)';
       for iO=1:order
           sigma(:,:,iMix)=sigma(:,:,iMix)-accSig2(:,:,iO,iMix)*diag(ar(:,iO,iMix))./accW(iMix) - diag(ar(:,iO,iMix))*accSig2(:,:,iO,iMix)'./accW(iMix);
         %  sigma(:,:,iMix)=sigma(:,:,iMix)+ diag(miu(:,iMix))*accSig3(:,:,order+1,iO)*diag(ar(:,iO,iMix))./accW(iMix);
           %   sigma(:,:,iMix)=sigma(:,:,iMix)+ diag(ar(:,iO,iMix))*accSig3(:,:,iO,order+1)*diag(miu(:,iMix))./accW(iMix);
           for jO=1:order
         %      sigma(:,:,iMix)=sigma(:,:,iMix) +  diag(ar(:,iO,iMix)) * accSig3(:,:,iO,jO,iMix) *diag(ar(:,jO,iMix))./accW(iMix);
           end
       end
        varx=diag(sigma(1:l/2,1:l/2,iMix));
       vary=diag(sigma(l/2+1:l,l/2+1:l,iMix));
       varxy=diag(sigma(1:l/2,l/2+1:l,iMix));
       sigma(:,:,iMix)=[diag(varx) diag(varxy);diag(varxy) diag(vary) ];
       
%         sigma(:,:,iMix)=accSig(:,:,iMix)/accW(iMix)-miu(:,iMix)*miu(:,iMix)';
%      for iO=1:order
%             tempSig1=miu(:,iMix)*accSig1(:,iO,iMix)'*diag(ar(:,iO,iMix))/accW(iMix);
%             sigma(:,:,iMix)=sigma(:,:,iMix)-tempSig1-tempSig1';
%             for jO=1:order
%                  sigma(:,:,iMix)=sigma(:,:,iMix)-diag(ar(:,iO,iMix))*accSig2(:,:,iO,jO,iMix)*diag(ar(:,jO,iMix))'/accW(iMix);
%             end
%      end
      [D_,dp]=chol(sigma(:,:,iMix));
      nFloor=1;
     if dp>0
         fprintf(1,'add variance floor:%d\n',nFloor);
        sigma(:,:,iMix)= diag(vFloor)+sigma(:,:,iMix);
        [D_,dp]=chol(sigma(:,:,iMix));
        nFloor=nFloor+1;
      end
    end
    
end


%% write gmm
ogmm=fopen([trainDir '/neutral_' emotion '.tone.gmm_ar_full'],'w');
for iMix=1:mixNum
    fprintf(ogmm,'%f\n',w(iMix));
end
for iMix=1:mixNum
    fprintf(ogmm,'%f\n',miu(:,iMix));
    for iO=1:order
        fprintf(ogmm,'%f\n',ar(:,iO,iMix));
    end
    fprintf(ogmm,'%f\n',sigma(:,:,iMix));
end
