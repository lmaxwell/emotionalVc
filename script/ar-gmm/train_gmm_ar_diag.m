clear;
trainDir='../../train/tone_nodynamic_diag/sad/3mix6mix_3';
gmmFile=[trainDir '/neutral_sad.tone.gmm'];
trainFile=[trainDir '/neutral_sad.tone_ar'];
l=16;
mixNum=6;
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
 

   
    % set ar coefficient
    ar=zeros(l,order,mixNum);
   

    
    
%% train ar combined gmm

%first run  
trainData=importdata(trainFile)';
indexOfSen=find(trainData(1,:)==0);
numOfSen=length(indexOfSen);
indexOfSen=[indexOfSen  size(trainData,2)+1];
beginZero=zeros(l,1);

A=zeros(l,order+1,order+1,mixNum);
B=zeros(l,order+1,mixNum);
mean=zeros(l,1);

postProb=zeros(mixNum,1);

liklihoodProb=0;
numOfSeg=0;
accW=zeros(mixNum,1);
accSig=zeros(l,mixNum);
accSig1=zeros(l,mixNum);
accSig2=zeros(l,order,mixNum);
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
                A(  :,iOd,jOd , iMix)= A(  :,iOd,jOd ,iMix)+ postProb(iMix)*preSegs(:,iOd).*preSegs(:,jOd);  %iOd th function , jOd th parameter
                end
                B(:,iOd,iMix)=B(:,iOd,iMix)+postProb(iMix)*trainData(:,start+iSeg).*preSegs(:,iOd);  %iOd th function
                A(:,order+1,iOd,iMix)=A(:,order+1,iOd,iMix)+postProb(iMix).*preSegs(:,iOd);  % order+1 th fuction, iOd th parameter
                A(:,iOd,order+1,iMix)=A(:,iOd,order+1,iMix)+ postProb(iMix).*preSegs(:,iOd);         % iOd th function , order+1 th parameter
            end
               A(:,order+1,order+1,iMix)=   A(:,order+1,order+1,iMix)+postProb(iMix);   %order+1 th function ,order+1 th parameter
               B(:,order+1,iMix)=B(:,order+1,iMix)+trainData(:,start+iSeg).*postProb(iMix); % order+1th function
                
                    accSig(:,iMix)=accSig(:,iMix)+postProb(iMix)*trainData(:,start+iSeg).*trainData(:,start+iSeg);
               
                   accSig1(:,iMix)=accSig1(:,iMix)+postProb(iMix)*trainData(:,start+iSeg);
         for iO=1:order
                accSig2(:,iO,iMix)=accSig2(:,iO,iMix)+postProb(iMix)*trainData(:,start+iSeg).*preSegs(:,iO);
           
         end
               
        end
         mean(:)=mean(:)+trainData(:,start+iSeg);  %for Sigma update
        
        vFloor=vFloor+(trainData(:,start+iSeg)-miu(:,maxMix)).* (trainData(:,start+iSeg)-miu(:,maxMix));
       
        numOfSeg=numOfSeg+1;
        
    end
    mean=mean/numOfSeg;
    

    
    
end
vFloor=vFloor*vfloor/numOfSeg;
    for iMix=1:mixNum
        for iL=1:l
              tempA=squeeze(A(iL,:,:,iMix));
            if min(svd(tempA))<1e-7
                temp=pinv(tempA)* B(iL,:,iMix)';
            else
                temp=tempA \ B(iL,:,iMix)';
            end
           
            ar(iL,:,iMix)=temp(1:order);
            miu(iL,iMix)=temp(order+1);
        end
       w(iMix)=accW(iMix)/numOfSeg;
         sigma(:,:,iMix)=diag(accSig(:,iMix)-accSig1(:,iMix).*miu(:,iMix))./accW(iMix);
     for iO=1:order  
           
            sigma(:,:,iMix)=sigma(:,:,iMix)-diag(accSig2(:,iO,iMix).*ar(:,iO,iMix))./accW(iMix);
            
%             for dim=1:l
%               
%                     if(sigma(dim,dim,iMix))<vFloor(dim)
%                         sigma(dim,dim,iMix)=vFloor(dim);
%                     end
%            
%             end
     end
     
%       [D_,p_]=chol(sigma(:,:,iMix));
%       if p_>0
%         sigma(:,:,iMix)=sigma(:,:,iMix)+ diag(vFloor);
%       end
%       
      
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
while(abs((curProb-preProb)/preProb)>1e-5)
liklihoodProb=0;
 numOfSeg=0;
 accW=zeros(mixNum,1);
 A=zeros(l,order+1,order+1,mixNum);
B=zeros(l,order+1,mixNum);
mean=zeros(l,1);
accSig=zeros(l,mixNum);
accSig1=zeros(l,mixNum);
accSig2=zeros(l,order,mixNum);
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
                end
                B(:,iOd,iMix)=B(:,iOd,iMix)+postProb(iMix)*trainData(:,start+iSeg).*preSegs(:,iOd);  %iOd th function
                A(:,order+1,iOd,iMix)=A(:,order+1,iOd,iMix)+postProb(iMix).*preSegs(:,iOd);  % order+1 th fuction, iOd th parameter
                A(:,iOd,order+1,iMix)=A(:,iOd,order+1,iMix)+ postProb(iMix).*preSegs(:,iOd);         % iOd th function , order+1 th parameter
            end
               A(:,order+1,order+1,iMix)=   A(:,order+1,order+1,iMix)+postProb(iMix);   %order+1 th function ,order+1 th parameter
               B(:,order+1,iMix)=B(:,order+1,iMix)+trainData(:,start+iSeg).*postProb(iMix); % order+1th function
               
              
               accSig(:,iMix)=accSig(:,iMix)+postProb(iMix)*trainData(:,start+iSeg).*trainData(:,start+iSeg);
               
                   accSig1(:,iMix)=accSig1(:,iMix)+postProb(iMix)*trainData(:,start+iSeg);
         for iO=1:order
                accSig2(:,iO,iMix)=accSig2(:,iO,iMix)+postProb(iMix)*trainData(:,start+iSeg).*preSegs(:,iO);
           
         end
               
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
       if accW(iMix)>0
        for iL=1:l
            tempA=squeeze(A(iL,:,:,iMix));
            if min(svd(tempA))<1e-7
                temp=pinv(tempA)* B(iL,:,iMix)';
            else
                temp=tempA \ B(iL,:,iMix)';
            end
            
            ar(iL,:,iMix)=temp(1:order);
            miu(iL,iMix)=temp(order+1);
        end
         w(iMix)=accW(iMix)/numOfSeg;
         sigma(:,:,iMix)=diag(accSig(:,iMix)-accSig1(:,iMix).*miu(:,iMix))./accW(iMix);
     for iO=1:order  
           
            sigma(:,:,iMix)=sigma(:,:,iMix)-diag(accSig2(:,iO,iMix).*ar(:,iO,iMix))./accW(iMix);
            
%             for dim=1:l
%               
%                     if(sigma(dim,dim,iMix))<vFloor(dim)
%                         sigma(dim,dim,iMix)=vFloor(dim);
%                     end
%            
%             end
     end
     
%       [D_,dp]=chol(sigma(:,:,iMix));
%       nFloor=1;
%        if dp>0
%          fprintf(1,'add variance floor:%d\n',nFloor);
%         sigma(:,:,iMix)= diag(vFloor)+sigma(:,:,iMix);
%         [D_,dp]=chol(sigma(:,:,iMix));
%         nFloor=nFloor+1;
%       end
    end
   end
    %
end

