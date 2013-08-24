oldSigma=sigma;
for iMix=1:mixNum
    sigma(:,:,iMix)=(sigma(:,:,iMix)+sigma(:,:,iMix)')/2;
end
for iMix=1:mixNum
    gconst(iMix)=log(w(iMix))-l/2*log(2*pi)-0.5*log(det(sigma(:,:,iMix)));
    invSigma(:,:,iMix)=inv( sigma(:,:,iMix) );
end
preProb=liklihoodProb/numOfSeg;
curProb=0;
 fprintf(1,'initialProb %f\n',preProb);
  fprintf(logfd,'initialProb %f\n',preProb);
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
for iSen=1:numOfSen
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
               
              
               accSig(:,:,iMix)=accSig(:,:,iMix)+postProb(iMix)*trainData(:,start+iSeg)*trainData(:,start+iSeg)';
              
                accSig1(:,iMix)=accSig1(:,iMix)+postProb(iMix)*trainData(:,start+iSeg);
              
                       
                       for iO=1:order
                            accSig2(:,:,iO,iMix)=accSig2(:,:,iO,iMix)+postProb(iMix)*preSegs(:,iO)*trainData(:,start+iSeg)';
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
    
      fprintf(1,'total liklihood:%f , avg liklihood:%f\n',liklihoodProb,liklihoodProb/numOfSeg);
        fprintf(logfd,'total liklihood:%f , avg liklihood:%f\n',liklihoodProb,liklihoodProb/numOfSeg);
    
    
end
preProb=curProb;

 curProb=liklihoodProb/numOfSeg;
mean=mean/numOfSeg;
 fprintf(1,'curProb %f\n',curProb);
  fprintf(logfd,'curProb %f\n',curProb);
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
         
        sigma(:,:,iMix)=( accSig(:,:,iMix)-miu(:,iMix)*accSig1(:,iMix)')/accW(iMix);
       for iO=1:order
           sigma(:,:,iMix)=sigma(:,:,iMix)-diag(ar(:,iO,iMix))*accSig2(:,:,iO,iMix)./accW(iMix);
       end
         
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