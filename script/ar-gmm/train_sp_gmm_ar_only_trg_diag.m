clear;
nMix=128;
emotion='sad';
trainDir=['../../vc/train/to' emotion '_ar/gmm_jde/' num2str(nMix) '-mix'];
gmmMean=[trainDir '/neutral-' emotion num2str(nMix) '.mean.txt' ];
gmmCov=[trainDir '/neutral-' emotion num2str(nMix) '.cov.txt' ];
gmmWght=[trainDir '/neutral-' emotion num2str(nMix) '.wght.txt' ];

trainFile=[trainDir '/../../dtw/' num2str(nMix) '-mix/' emotion '.mat.txt'];
l=24;
nMix=25;
mixNum=nMix;
order=1;
vfloor=0.001;

logfd=fopen([trainDir '/log_diag_trg_order' num2str(order)],'w');
%% read gmm model
%  w=importdata(gmmWght);
%  miu=reshape( importdata(gmmMean) ,l,mixNum);
% sigma=reshape( importdata(gmmCov),l,l,mixNum);
% 
% miux=miu(1:l/2,:);
% miuy=miu(l/2+1:l,:);
% ar=zeros(l,order,mixNum);
% ar2=zeros(l,order,mixNum);
% 
% preMiu=zeros(l,order,mixNum);
% gconst=zeros(mixNum,1);
% invSigma=zeros(l,l,mixNum);
% for iMix=1:mixNum
%     
%    
%     gconst(iMix)=log(w(iMix))-l/2*log(2*pi)-0.5*log(det(sigma(:,:,iMix)));
%     
%     invSigma(:,:,iMix)=inv( sigma(:,:,iMix) );
% end
 %%
%     w=gmmP(1:mixNum);
%     miux=zeros(l/2,mixNum);
%     miuy=zeros(l/2,mixNum);
%     miu=zeros(l,mixNum);
%     sigma=zeros(l,l,mixNum);
%     beginModel.w=w;
%     for i=1:mixNum
%         miux(:,i)=gmmP(mixNum+(i-1)*(l^2+l)+1:mixNum+(i-1)*(l^2+l)+l/2);
%         miuy(:,i)=gmmP(mixNum+(i-1)*(l^2+l)+l/2+1:mixNum+(i-1)*(l^2+l)+l);
%         miu(:,i)=[miux(:,i);miuy(:,i)];
%         beginModel.miu(:,i)=[miux(:,i);miuy(:,i)];
%         sigma(:,:,i)=reshape(gmmP(mixNum+(i-1)*(l^2+l)+l+1:mixNum+(i-1)*(l^2+l)+l+l^2),l,l);
%         beginModel.sigma(:,:,i)=sigma(:,:,i);
%     end
% 
%     sigmaxx=sigma(1:l/2,1:l/2,:);
%     sigmayy=sigma(l/2+1:l,l/2+1:l,:);
%     sigmaxy=sigma(1:l/2,l/2+1:l,:);
%     sigmayx=sigma(l/2+1:l,1:l/2,:);
%    
%     % set ar coefficient
%     ar=zeros(l,order,mixNum);
%    

    
    
%% train ar combined gmm

%first run  
trainData=importdata(trainFile)';


%%
Data=trainData(:,trainData(1,:)~=0);
Data=Data(:,1:10000);
nbStates=nMix;
[nbVar, nbData] = size(Data);

[Data_id, Centers] = kmeans(Data', nbStates);
Mu = Centers';
for i=1:nbStates
  idtmp = find(Data_id==i);
  Priors(i) = length(idtmp);
  Sigma(:,:,i) = cov([Data(:,idtmp) Data(:,idtmp)]');
  %% Add a tiny variance to avoid numerical instability
  Sigma(:,:,i) = Sigma(:,:,i) + 1E-5.*diag(ones(nbVar,1));
end
Priors = Priors ./ sum(Priors);

%%
w=Priors';
miu=Mu;
sigma=Sigma;
for iMix=1:mixNum
    
   
    gconst(iMix)=log(w(iMix))-l/2*log(2*pi)-0.5*log(det(sigma(:,:,iMix)));
    
    invSigma(:,:,iMix)=inv( sigma(:,:,iMix) );
end
%%
indexOfSen=find(trainData(1,:)==0);
numOfSen=length(indexOfSen);
indexOfSen=[indexOfSen  size(trainData,2)+1];

%%
beginZero=zeros(l,1);

A=zeros(l,order,order,mixNum);
B=zeros(l,order,mixNum);
mean=zeros(l,1);

accPre=zeros(l,order,mixNum);
accCur=zeros(l,mixNum);

postProb=zeros(mixNum,1);

liklihoodProb=0;
numOfSeg=0;
accW=zeros(mixNum,1);
accSig=zeros(l,mixNum);
accSig1=zeros(l,mixNum);
accSig2=zeros(l,order,mixNum);
vFloor=zeros(l,1);

for iSen=1:numOfSen/20
    start=indexOfSen(iSen)+1;
    xend=indexOfSen(iSen+1)-1;
    fprintf( 1 , '%dth sentence,start:%d,end:%d\n' , iSen , start , xend);
     fprintf( logfd , '%dth sentence,start:%d,end:%d\n' , iSen , start , xend);
    for iSeg=0:xend-start
  %      fprintf(1,'%dth segment\n',iSeg);
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
           
            for iOd=1:order  
                for jOd=1:order 
                A(  :,iOd,jOd , iMix)= A(  :,iOd,jOd ,iMix)+ postProb(iMix)*preSegs(:,iOd).*preSegs(:,jOd);  %iOd th function , jOd th parameter
                end
                B(:,iOd,iMix)=B(:,iOd,iMix)+postProb(iMix)*trainData(:,start+iSeg).*preSegs(:,iOd);  %iOd th function
               accPre(:,iOd,iMix)=accPre(:,iOd,iMix)+postProb(iMix)*preSegs(:,iOd);
            end
                accCur(:,iMix)=accCur(:,iMix)+postProb(iMix)*trainData(:,start+iSeg);          
           
                
                 accSig(:,iMix)=accSig(:,iMix)+postProb(iMix)*trainData(:,start+iSeg).*trainData(:,start+iSeg);
               
               %    accSig1(:,iMix)=accSig1(:,iMix)+postProb(iMix)*trainData(:,start+iSeg);
         for iO=1:order
                accSig2(:,iO,iMix)=accSig2(:,iO,iMix)+postProb(iMix)*trainData(:,start+iSeg).*preSegs(:,iO);
           
         end
                  
              

               
        end     
%             for iOd=1:order  
%                 for jOd=1:order 
%                 A(  :,iOd,jOd , :)= A(  :,iOd,jOd ,:)+ preSegs(:,iOd).*preSegs(:,jOd)*postProb';  %iOd th function , jOd th parameter
%                 end
%                 B(:,iOd,:)=B(:,iOd,:)+trainData(:,start+iSeg).*preSegs(:,iOd)*postProb';  %iOd th function
%                 A(:,order+1,iOd,:)=A(:,order+1,iOd,:)+preSegs(:,iOd)*postProb';  % order+1 th fuction, iOd th parameter
%                 A(:,iOd,order+1,:)=A(:,iOd,order+1,:)+ preSegs(:,iOd)*postProb;         % iOd th function , order+1 th parameter
%             end
%                A(:,order+1,order+1,:)=   A(:,order+1,order+1,:)+ones(l,1)*postProb';   %order+1 th function ,order+1 th parameter
%                B(:,order+1,:)=B(:,order+1,:)+trainData(:,start+iSeg)*postProb'; % order+1th function
%                 
%                for iMix=1:mixNum
%                accSig(:,:,iMix)=accSig(:,:,iMix)+postProb(iMix)*trainData(:,start+iSeg)*trainData(:,start+iSeg)';
%                
%                for iL=1:l
%                    accSig1(iL,iMix)=accSig1(iL,iMix)+postProb(iMix)*trainData(iL,start+iSeg);
%                    for jL=1:l
%                        
%                        for iO=1:order
%                             accSig2(iL,jL,iO,iMix)=accSig2(iL,jL,iO,iMix)+postProb(iMix)*trainData(iL,start+iSeg)*preSegs(jL,iO);
%                        end
%                    end
%                end
% 
%                end
        
         mean(:)=mean(:)+trainData(:,start+iSeg);  %for Sigma update
        
        vFloor=vFloor+(trainData(:,start+iSeg)-miu(:,maxMix)).* (trainData(:,start+iSeg)-miu(:,maxMix));
       
        numOfSeg=numOfSeg+1;
        
    end
    mean=mean/numOfSeg;
    

    fprintf(1,'total liklihood:%f , avg liklihood:%f\n',liklihoodProb,liklihoodProb/numOfSeg);
        fprintf(logfd,'total liklihood:%f , avg liklihood:%f\n',liklihoodProb,liklihoodProb/numOfSeg);
end

vFloor=vFloor*vfloor/numOfSeg;


%%
for iMix=1:mixNum
    for iO=1:order
        for jO=1:order
            A(:,iO,jO,iMix)=A(:,iO,jO,iMix)/accW(iMix)-accPre(:,iO,iMix)./accW(iMix).*accPre(:,jO,iMix)./accW(iMix);
        end
        
        B(:,iO,iMix)=B(:,iO,iMix)./accW(iMix)-accCur(:,iMix)./accW(iMix).*accPre(:,iO,iMix)./accW(iMix);
    end
end

    for iMix=1:mixNum
        for iL=1:l
              tempA=squeeze(A(iL,:,:,iMix));
            if min(svd(tempA))<1e-6
                fprintf(1,'%d,%d svd<1e-6\n',iMix,iL);
                temp=pinv(tempA)* B(iL,:,iMix)';
            else
                temp=tempA \ B(iL,:,iMix)';
            end
           for iO=1:order
               if abs(temp(iO))<1
                    ar(iL,iO,iMix)=temp(iO);
               else
                    ar(iL,iO,iMix)=temp(iO)/(abs(temp(iO))+1);
               end
           end
           ar2(iL,:,iMix)=temp(1:order);
        end
        
        miu(:,iMix)=accCur(:,iMix)/accW(iMix);
        preMiu(:,:,iMix)=accPre(:,:,iMix)/accW(iMix);
          
        
        
       w(iMix)=accW(iMix)/numOfSeg;
             sigma(:,:,iMix)=diag(accSig(:,iMix)./accW(iMix)-    accCur(:,iMix)./accW(iMix).*accCur(:,iMix)./accW(iMix)  );
     for iO=1:order  
           
            sigma(:,:,iMix)=sigma(:,:,iMix)-diag(B(:,iO,iMix).*ar(:,iO,iMix));
     end
%      sigma(:,:,iMix)=accSig(:,:,iMix)/accW(iMix)-miu(:,iMix)*miu(:,iMix)';
%      for iO=1:order
%             tempSig1=miu(:,iMix)*accSig1(:,iO,iMix)'*diag(ar(:,iO,iMix))/accW(iMix);
%             sigma(:,:,iMix)=sigma(:,:,iMix)-tempSig1-tempSig1';
%             for jO=1:order
%                  sigma(:,:,iMix)=sigma(:,:,iMix)-diag(ar(:,iO,iMix))*accSig2(:,:,iO,jO,iMix)*diag(ar(:,jO,iMix))'/accW(iMix);
%             end
%      end
       sigma(:,:,iMix)= (sigma(:,:,iMix)+sigma(:,:,iMix)')/2;
      [D_,p_]=chol(sigma(:,:,iMix));
      if p_>0
        sigma(:,:,iMix)=sigma(:,:,iMix)+ diag(vFloor);
      end
      
    gconst(iMix)=log(w(iMix))-l/2*log(2*pi)-0.5*log(det(sigma(:,:,iMix)));
    invSigma(:,:,iMix)=inv( sigma(:,:,iMix) );
    end
  
 %%
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
  fprintf(logfd,'initialProb %f\n',preProb);
%iteration
ite=1;
while(abs((curProb-preProb)/preProb)>1e-5 && ite<20)
liklihoodProb=0;
 numOfSeg=0;
A=zeros(l,order,order,mixNum);
B=zeros(l,order,mixNum);
mean=zeros(l,1);

accPre=zeros(l,order,mixNum);
accCur=zeros(l,mixNum);

postProb=zeros(mixNum,1);

accW=zeros(mixNum,1);
accSig=zeros(l,mixNum);
accSig1=zeros(l,mixNum);
accSig2=zeros(l,order,mixNum);
for iSen=1:numOfSen/20
    start=indexOfSen(iSen)+1;
    xend=indexOfSen(iSen+1)-1;
      fprintf( 1 , '%dth sentence,start:%d,end:%d\n' , iSen , start , xend);
     fprintf( logfd , '%dth sentence,start:%d,end:%d\n' , iSen , start , xend);
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
         postProb=zeros(mixNum,1);
        for iMix=1:mixNum

            if iSeg<0
                curMiu=beginModel.miu(:,iMix);
                 postProb(iMix)= w(iMix).*1./((2*pi)^(l/2).*sqrt(det(beginModel.sigma(:,:,iMix)))).* ...
            exp(-0.5.*(trainData(:,start+iSeg)-curMiu)'/beginModel.sigma(:,:,iMix)*(trainData(:,start+iSeg)-curMiu));
            else
                curMiu=miu(:,iMix);
                for iO=1:order
                    curMiu=curMiu+ar(:,iO,iMix).*(preSegs(:,iO)-preMiu(:,iO,iMix));
                end
              
            %.*1./((2*pi)^(l/2).*sqrt(det(sigma(:,:,iMix)))).* ...
         %exp(-0.5.*(trainData(:,start+iSeg)-miu(:,iMix))'/sigma(:,:,iMix)*(trainData(:,start+iSeg)-miu(:,iMix)));
     %   log(w(i))-0.5*log(det(sigmaxx(:,:,i)))-4*log(2*pi)-0.5.*(x-miux(:,i))'/sigmaxx(:,:,i)*(x-miux(:,i));
            end
            postProb(iMix)=log(w(iMix))-l/2*log(2*pi)-0.5*log(det(sigma(:,:,iMix)))-0.5.*(trainData(:,start+iSeg)-curMiu)'/sigma(:,:,iMix)*(trainData(:,start+iSeg)-curMiu)          ;
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
                end
                B(:,iOd,iMix)=B(:,iOd,iMix)+postProb(iMix)*trainData(:,start+iSeg).*preSegs(:,iOd);  %iOd th function
               accPre(:,iOd,iMix)=accPre(:,iOd,iMix)+postProb(iMix)*preSegs(:,iOd);
            end
                accCur(:,iMix)=accCur(:,iMix)+postProb(iMix)*trainData(:,start+iSeg);          
           
                
                 accSig(:,iMix)=accSig(:,iMix)+postProb(iMix)*trainData(:,start+iSeg).*trainData(:,start+iSeg);
               
               %    accSig1(:,iMix)=accSig1(:,iMix)+postProb(iMix)*trainData(:,start+iSeg);
         for iO=1:order
                accSig2(:,iO,iMix)=accSig2(:,iO,iMix)+postProb(iMix)*trainData(:,start+iSeg).*preSegs(:,iO);
           
         end
                  
              

               
        end  
        
         
        mean(:)=mean(:)+trainData(:,start+iSeg);  %for Sigma update
        numOfSeg=numOfSeg+1;
        
         
    end
    
      fprintf(1,'total liklihood:%f , avg liklihood:%f\n',liklihoodProb,liklihoodProb/numOfSeg);
        fprintf(logfd,'total liklihood:%f , avg liklihood:%f\n',liklihoodProb,liklihoodProb/numOfSeg);
    
    
end
preProb=curProb;

 curProb=liklihoodProb/numOfSeg;
mean=mean/numOfSeg;
 fprintf(1,'curProb %f\n',curProb);
  fprintf(logfd,'curProb %f\n',curProb);
for iMix=1:mixNum
    for iO=1:order
        for jO=1:order
            A(:,iO,jO,iMix)=A(:,iO,jO,iMix)/accW(iMix)-accPre(:,iO,iMix)./accW(iMix).*accPre(:,jO,iMix)./accW(iMix);
        end
        
        B(:,iO,iMix)=B(:,iO,iMix)./accW(iMix)-accCur(:,iMix)./accW(iMix).*accPre(:,iO,iMix)./accW(iMix);
    end
end

    for iMix=1:mixNum
        for iL=1:l
              tempA=squeeze(A(iL,:,:,iMix));
            if min(svd(tempA))<1e-6
                fprintf(1,'%d,%d svd<1e-6\n',iMix,iL);
                temp=pinv(tempA)* B(iL,:,iMix)';
            else
                temp=tempA \ B(iL,:,iMix)';
            end
           for iO=1:order
               if abs(temp(iO))<1
                    ar(iL,iO,iMix)=temp(iO);
               else
                    ar(iL,iO,iMix)=temp(iO)/(abs(temp(iO))+1);
               end
           end
           ar2(iL,:,iMix)=temp(1:order);
        end
        
        miu(:,iMix)=accCur(:,iMix)/accW(iMix);
        preMiu(:,:,iMix)=accPre(:,:,iMix)/accW(iMix);
          
        
        
       w(iMix)=accW(iMix)/numOfSeg;
             sigma(:,:,iMix)=diag(accSig(:,iMix)./accW(iMix)-    accCur(:,iMix)./accW(iMix).*accCur(:,iMix)./accW(iMix)  );
     for iO=1:order  
           
            sigma(:,:,iMix)=sigma(:,:,iMix)-diag(B(:,iO,iMix).*ar(:,iO,iMix));
     end
%      sigma(:,:,iMix)=accSig(:,:,iMix)/accW(iMix)-miu(:,iMix)*miu(:,iMix)';
%      for iO=1:order
%             tempSig1=miu(:,iMix)*accSig1(:,iO,iMix)'*diag(ar(:,iO,iMix))/accW(iMix);
%             sigma(:,:,iMix)=sigma(:,:,iMix)-tempSig1-tempSig1';
%             for jO=1:order
%                  sigma(:,:,iMix)=sigma(:,:,iMix)-diag(ar(:,iO,iMix))*accSig2(:,:,iO,jO,iMix)*diag(ar(:,jO,iMix))'/accW(iMix);
%             end
%      end
       sigma(:,:,iMix)= (sigma(:,:,iMix)+sigma(:,:,iMix)')/2;
      [D_,p_]=chol(sigma(:,:,iMix));
      if p_>0
        sigma(:,:,iMix)=sigma(:,:,iMix)+ diag(vFloor);
      end
      
    gconst(iMix)=log(w(iMix))-l/2*log(2*pi)-0.5*log(det(sigma(:,:,iMix)));
    invSigma(:,:,iMix)=inv( sigma(:,:,iMix) );
    end
    ite=ite+1;
end


%% write gmm
% ogmm=fopen([trainDir '/neutral_' emotion '.tone.gmm_ar_full'],'w');
% for iMix=1:mixNum
%     fprintf(ogmm,'%f\n',w(iMix));
% end
% for iMix=1:mixNum
%     fprintf(ogmm,'%f\n',miu(:,iMix));
%     for iO=1:order
%         fprintf(ogmm,'%f\n',ar(:,order,iMix));
%     end
%     fprintf(ogmm,'%f\n',sigma(:,:,iMix));
% end
