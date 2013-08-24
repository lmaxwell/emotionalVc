% clear;
% emotion='sad';
% pMix=6;
% sMix=12;
% trainDir=['../../train/tone_nodynamic_diag/' emotion '/' num2str(pMix) 'mix' num2str(sMix) 'mix_3'];
% gmmFile=[trainDir '/neutral_' emotion '.phrase.gmm_full'];
% trainFile=[trainDir '/neutral_' emotion '.phrase_ar'];
% l=6;
% mixNum=pMix;
% order=1;
% vfloor=0.001;
% addpath('lightspeed');
% %% read gmm model
%  gmmP=importdata(gmmFile);
% 
% 
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
%    gconst=zeros(1,mixNum)';
%     % set ar coefficient
%     ar=zeros(l,order,mixNum);
%     
%     sigmaxx=sigma(1:l/2,1:l/2,:);
%     sigmayy=sigma(l/2+1:l,l/2+1:l,:);
%     sigmaxy=sigma(1:l/2,l/2+1:l,:);
%     sigmayx=sigma(l/2+1:l,1:l/2,:);
%  
%   for iMix=1:mixNum
%     
%    var(:,iMix)=diag(sigma(:,:,iMix));
%     gconst(iMix)=log(w(iMix))-l/2*log(2*pi)-0.5*log(prod(var(:,iMix)));
%     
%    
% end 
% %    
%     % set ar coefficient
%     ar=zeros(l,order,mixNum);
%    arx=zeros(l/2,order,mixNum);
% ary=zeros(l/2,order,mixNum);
% bt=zeros(l/2,l/2,mixNum);
%     
%     
% %% train ar combined gmm
% 
% %first run  
% trainData=importdata(trainFile)';
% 
% %%
% % Data=trainData(:,trainData(1,:)~=0);
% % %Data=Data(:,1:10000);
% % nbStates=mixNum;
% % [nbVar, nbData] = size(Data);
% % 
% % [Data_id, Centers] = kmeans(Data', nbStates);
% % Mu = Centers';
% % for i=1:nbStates
% %   idtmp = find(Data_id==i);
% %   Priors(i) = length(idtmp);
% %   Sigma(:,:,i) = cov([Data(:,idtmp) Data(:,idtmp)]');
% %   %% Add a tiny variance to avoid numerical instability
% %  % Sigma(:,:,i) = Sigma(:,:,i) + 1E-5.*diag(ones(nbVar,1));
% % end
% % Priors = Priors ./ sum(Priors);
% % 
% % %%
% % w=Priors';
% % miu=Mu;
% % sigma=Sigma;
% 
% %%
% indexOfSen=find(trainData(1,:)==0);
% numOfSen=length(indexOfSen);
% indexOfSen=[indexOfSen  size(trainData,2)+1];
% beginZero=zeros(l,1);
% 
% A=zeros(order,order,l/2,mixNum);
% C=zeros(mixNum,order,l/2);
% B=zeros(order+l/2,order+l/2,l/2,mixNum);
% D=zeros(mixNum,order+l/2,l/2);
% mean=zeros(l,1);
% 
% postProb=zeros(mixNum,1);
% 
% liklihoodProb=0;
% numOfSeg=0;
% accW=zeros(mixNum,1);
% accSigx=zeros(l/2,mixNum);
% accSigx1=zeros(l/2,mixNum);
% accSigx2=zeros(l/2,mixNum,order);
% accSigy=zeros(l/2,mixNum);
% accSigy1=zeros(l/2,mixNum);
% accSigy2=zeros(l/2,mixNum,order);
% accSigy3=zeros(l/2,l/2,mixNum);
% accPre=zeros(mixNum,l,order);
% vFloor=zeros(l,1);
% 
% for iSen=1:numOfSen
%     start=indexOfSen(iSen)+1;
%     xend=indexOfSen(iSen+1)-1;
%     fprintf( 1 , '%dth sentence,start:%d,end:%d\n' , iSen , start , xend);
%     for iSeg=0:xend-start
%         fprintf(1,'%dth segment\n',iSeg);
%         %set preSegs
%         if iSeg==0
%               preSegs=beginZero*ones(1,order-iSeg);
%         elseif(iSeg<order)
%            
%             preSegs=[beginZero*ones(1,order-iSeg) trainData(:,start+iSeg-1:-1:start)];
%            
%         else
%             preSegs=trainData(:,start+iSeg-1:-1:start+iSeg-order);
%        end
%         
%         sum_postProb=-1.0E10;
%         maxProb=-1.0E10;
%      
%         postProb=gconst- ndsum(0.5.*(trainData(:,start+iSeg)*ones(1,mixNum)-miu).*(trainData(:,start+iSeg)*ones(1,mixNum)-miu)./var, 1) ;
%           maxMix=find(postProb==max(postProb));
%         sum_postProb=logsumexp(postProb);
%         postProb=exp(postProb-sum_postProb);
%          accW=accW+postProb;
%         liklihoodProb=liklihoodProb+sum_postProb;
%         
%         
%         for iL=1:l/2
%             s=[preSegs(iL,:) ]';
%             r=[preSegs(l/2+iL,:)  trainData(1:l/2,start+iSeg)' ]';
%             C(:,:,iL)=C(:,:,iL)+ postProb*s'.*trainData(iL,start+iSeg);
%             D(:,:,iL)=D(:,:,iL)+postProb*r'.*trainData(l/2+iL,start+iSeg);
%             for iMix=1:mixNum
%                 A(:,:,iL,iMix )= A(:,:,iL,iMix) + postProb(iMix).*s*s';
%                 B(:,:,iL,iMix)=B(:,:,iL,iMix) + postProb(iMix).*r*r';
%             end
%         end
%         
%         
%         
%         accSigx=accSigx+trainData(1:l/2,start+iSeg).*trainData(1:l/2,start+iSeg)*postProb';
%         accSigx1=accSigx1+trainData(1:l/2,start+iSeg)*postProb';
%           accSigy=accSigy+trainData(l/2+1:l,start+iSeg).*trainData(l/2+1:l,start+iSeg)*postProb';
%           accSigy1=accSigy1+trainData(l/2+1:l,start+iSeg)*postProb';
%            %accSigy3=accSigy3+trainData(l/2+1:l,start+iSeg).*trainData(1:l/2,start+iSeg)*postProb';
%        for iMix=1:mixNum
%             accSigy3(:,:,iMix)= accSigy3(:,:,iMix)+postProb(iMix).*trainData(l/2+1:l,start+iSeg)*trainData(1:l/2,start+iSeg)';
%        end
%          for iO=1:order
%                 accSigx2(:,:,iO)=accSigx2(:,:,iO)+trainData(1:l/2,start+iSeg).*preSegs(1:l/2,iO)*postProb';
%                 accSigy2(:,:,iO)=accSigy2(:,:,iO)+trainData(l/2+1:l,start+iSeg).*preSegs(l/2+1:l,iO)*postProb';
%                 accPre(:,:,iO)=accPre(:,:,iO)+postProb*preSegs(:,iO)';
%          end
%          
%         
%         
%         vFloor=vFloor+(trainData(:,start+iSeg)-miu(:,maxMix)).* (trainData(:,start+iSeg)-miu(:,maxMix));
%        
%         numOfSeg=numOfSeg+1;
%         
%     end
%     mean=mean/numOfSeg;
%     
% 
%     
%     
% end
% vFloor=vFloor*vfloor/numOfSeg;
% 
% %%
% 
%    
% for iL=1:l/2
%     for iMix=1:mixNum
%         preMiu(iMix,iL,1:order)=accPre(iMix,iL,1:order)./accW(iMix);
%          preMiu(iMix,iL+l/2,1:order)=accPre(iMix,iL+l/2,1:order)./accW(iMix);
%         s=accPre(iMix,iL,1:order)';
%         r=[accPre(iMix,iL+l/2,1:order) accSigx1(:,iMix)']';
%         A(:,:,iL,iMix)=A(:,:,iL,iMix)./accW(iMix)-s*s'./accW(iMix)./accW(iMix);
%          B(:,:,iL,iMix)=B(:,:,iL,iMix)./accW(iMix)-r*r'./accW(iMix)./accW(iMix);
%          C(iMix,:,iL)=C(iMix,:,iL)./accW(iMix)-accSigx1(iL,iMix).*s'./accW(iMix)./accW(iMix);
%          D(iMix,:,iL)=D(iMix,:,iL)./accW(iMix)-accSigy1(iL,iMix).*r'./accW(iMix)./accW(iMix);
%     end
% 
% end
% 
% 
% 
% 
% 
%     for iMix=1:mixNum
%         for iL=1:l/2
%               tempA=squeeze(A(:,:,iL,iMix));
%             if min(svd(tempA))<1e-7
%                 temp=pinv(tempA)* C(iMix,:,iL)';
%             else
%                 temp=tempA \ C(iMix,:,iL)';
%             end
%            
%             arx(iL,:,iMix)=temp(1:order);
%       
%             ar(iL,:,iMix)=temp(1:order);
%         end
%         for iL=1:l/2
%             tempB=squeeze( B(:,:,iL,iMix) );
%              if min(svd(tempB))<1e-7
%                 temp=pinv(tempB)* D(iMix,:,iL)';
%             else
%                 temp=tempB \ D(iMix,:,iL)';
%              end
%             ary(iL,:,iMix)=temp(1:order);
%               ar(l/2+iL,:,iMix)=temp(1:order);
%   
%             bt(iL,:,iMix)=temp(order+1:order+l/2);
%                %bt(iL,iMix)=0;
%         end
%         miux(:,iMix)=accSigx1(:,iMix)./accW(iMix);
%         miuy(:,iMix)=accSigy1(:,iMix)./accW(iMix);
%        
%        
%        w(iMix)=accW(iMix)/numOfSeg;
%          sigmaxx(:,:,iMix)=diag(accSigx(:,iMix)-accSigx1(:,iMix).*miux(:,iMix))./accW(iMix);
%          
%           sigmayy(:,:,iMix)=diag(accSigy(:,iMix)-accSigy1(:,iMix).*miuy(:,iMix))./accW(iMix) -squeeze(D(iMix,order+1:order+l/2,:))*bt(:,:,iMix)' ;
%           
%      for iO=1:order  
%            
%           sigmaxx(:,:,iMix)=sigmaxx(:,:,iMix)-diag(squeeze(C(iMix,iO,:)).*arx(:,iO,iMix) );
%             sigmayy(:,:,iMix)=sigmayy(:,:,iMix)-diag( squeeze(D(iMix,iO,:)).*ary(:,iO,iMix)  );
% 
%      end
%      sigma(:,:,iMix)=blkdiag(sigmaxx(:,:,iMix),sigmayy(:,:,iMix));
%       [D_,p_]=chol(sigma(:,:,iMix));
%       if p_>0
%         sigma(:,:,iMix)=sigma(:,:,iMix)+ diag(vFloor);
%       end
%          var(:,iMix)=diag(sigma(:,:,iMix));
%             gconst(iMix)=log(w(iMix))-l/2*log(2*pi)-0.5*log(prod( var(1:l/2,iMix) )) -0.5*log(prod(var(l/2+1:l,iMix)));
%       % gconst2(iMix)=log(w(iMix))-l/2*log(2*pi)-0.5*log(prod( var(:,iMix) ));
%     end



%%
preProb=liklihoodProb/numOfSeg;
curProb=0;
preProb=1;
 fprintf(1,'initialProb %f\n',preProb);
%iteration
ite=1;
while(abs((curProb-preProb)/preProb)>1e-5 && ite <30)
A=zeros(order,order,l/2,mixNum);
C=zeros(mixNum,order,l/2);
B=zeros(order+l/2,order+l/2,l/2,mixNum);
D=zeros(mixNum,order+l/2,l/2);
mean=zeros(l,1);

postProb=zeros(mixNum,1);

liklihoodProb=0;
numOfSeg=0;
accW=zeros(mixNum,1);
accSigx=zeros(l/2,mixNum);
accSigx1=zeros(l/2,mixNum);
accSigx2=zeros(l/2,mixNum,order);
accSigy=zeros(l/2,mixNum);
accSigy1=zeros(l/2,mixNum);
accSigy2=zeros(l/2,mixNum,order);
accSigy3=zeros(l/2,l/2,mixNum);
accPre=zeros(mixNum,l,order);
vFloor=zeros(l,1);

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
       
     curMiu=[miux ;miuy];
        for iO=1:order
            curMiu(1:l,:)=curMiu(1:l,:) + (preSegs(1:l,iO)*ones(1,mixNum)- squeeze(preMiu(:,:,iO))' ).*squeeze( ar(:,iO,:) );
            
        end
        for iMix=1:mixNum
            curMiu(l/2+1:l,iMix)=curMiu(l/2+1:l,iMix) + bt(:,:,iMix)*(trainData(1:l/2,start+iSeg)-miux(:,iMix));
        end
        postProb=gconst- ndsum(0.5.*(trainData(:,start+iSeg)*ones(1,mixNum)-curMiu).*(trainData(:,start+iSeg)*ones(1,mixNum)-curMiu)./var, 1) ;
        sum_postProb=logsumexp(postProb);
        postProb=exp(postProb-sum_postProb);
        accW=accW+postProb;
        liklihoodProb=liklihoodProb+sum_postProb;
        
         
        for iL=1:l/2
            s=[preSegs(iL,:) ]';
            r=[preSegs(l/2+iL,:)  trainData(1:l/2,start+iSeg)' ]';
            C(:,:,iL)=C(:,:,iL)+ postProb*s'.*trainData(iL,start+iSeg);
            D(:,:,iL)=D(:,:,iL)+postProb*r'.*trainData(l/2+iL,start+iSeg);
            for iMix=1:mixNum
                A(:,:,iL,iMix )= A(:,:,iL,iMix) + postProb(iMix).*s*s';
                B(:,:,iL,iMix)=B(:,:,iL,iMix) + postProb(iMix).*r*r';
            end
        end
        
        
        
        accSigx=accSigx+trainData(1:l/2,start+iSeg).*trainData(1:l/2,start+iSeg)*postProb';
        accSigx1=accSigx1+trainData(1:l/2,start+iSeg)*postProb';
          accSigy=accSigy+trainData(l/2+1:l,start+iSeg).*trainData(l/2+1:l,start+iSeg)*postProb';
          accSigy1=accSigy1+trainData(l/2+1:l,start+iSeg)*postProb';
           %accSigy3=accSigy3+trainData(l/2+1:l,start+iSeg).*trainData(1:l/2,start+iSeg)*postProb';
       for iMix=1:mixNum
            accSigy3(:,:,iMix)= accSigy3(:,:,iMix)+postProb(iMix).*trainData(l/2+1:l,start+iSeg)*trainData(1:l/2,start+iSeg)';
       end
         for iO=1:order
                accSigx2(:,:,iO)=accSigx2(:,:,iO)+trainData(1:l/2,start+iSeg).*preSegs(1:l/2,iO)*postProb';
                accSigy2(:,:,iO)=accSigy2(:,:,iO)+trainData(l/2+1:l,start+iSeg).*preSegs(l/2+1:l,iO)*postProb';
                accPre(:,:,iO)=accPre(:,:,iO)+postProb*preSegs(:,iO)';
         end
         
     
        numOfSeg=numOfSeg+1;
    end
    
   
    
    
end
preProb=curProb;

 curProb=liklihoodProb/numOfSeg;
mean=mean/numOfSeg;
 fprintf(1,'curProb %f\n',curProb);
 
 
  for iL=1:l/2
    for iMix=1:mixNum
   preMiu(iMix,iL,1:order)=accPre(iMix,iL,1:order)./accW(iMix);
         preMiu(iMix,iL+l/2,1:order)=accPre(iMix,iL+l/2,1:order)./accW(iMix);
        s=accPre(iMix,iL,1:order)';
        r=[accPre(iMix,iL+l/2,1:order) accSigx1(:,iMix)']';
        A(:,:,iL,iMix)=A(:,:,iL,iMix)./accW(iMix)-s*s'./accW(iMix)./accW(iMix);
         B(:,:,iL,iMix)=B(:,:,iL,iMix)./accW(iMix)-r*r'./accW(iMix)./accW(iMix);
         C(iMix,:,iL)=C(iMix,:,iL)./accW(iMix)-accSigx1(iL,iMix).*s'./accW(iMix)./accW(iMix);
         D(iMix,:,iL)=D(iMix,:,iL)./accW(iMix)-accSigy1(iL,iMix).*r'./accW(iMix)./accW(iMix);
    end

end





    for iMix=1:mixNum
        for iL=1:l/2
              tempA=squeeze(A(:,:,iL,iMix));
            if min(svd(tempA))<1e-7
                temp=pinv(tempA)* C(iMix,:,iL)';
            else
                temp=tempA \ C(iMix,:,iL)';
            end
           
            arx(iL,:,iMix)=temp(1:order);
      
            ar(iL,:,iMix)=temp(1:order);
        end
        for iL=1:l/2
            tempB=squeeze( B(:,:,iL,iMix) );
             if min(svd(tempB))<1e-7
                temp=pinv(tempB)* D(iMix,:,iL)';
            else
                temp=tempB \ D(iMix,:,iL)';
             end
            ary(iL,:,iMix)=temp(1:order);
              ar(l/2+iL,:,iMix)=temp(1:order);
  
            bt(iL,:,iMix)=temp(order+1:order+l/2);
               %bt(iL,iMix)=0;
        end
        miux(:,iMix)=accSigx1(:,iMix)./accW(iMix);
        miuy(:,iMix)=accSigy1(:,iMix)./accW(iMix);
        
       w(iMix)=accW(iMix)/numOfSeg;
         sigmaxx(:,:,iMix)=diag(accSigx(:,iMix)-accSigx1(:,iMix).*miux(:,iMix))./accW(iMix);
         
          sigmayy(:,:,iMix)=diag(accSigy(:,iMix)-accSigy1(:,iMix).*miuy(:,iMix))./accW(iMix) -squeeze(D(iMix,order+1:order+l/2,:))*bt(:,:,iMix)'./accW(iMix) ;
          
     for iO=1:order  
           
            sigmaxx(:,:,iMix)=sigmaxx(:,:,iMix)-diag(squeeze(C(iMix,iO,:)).*arx(:,iO,iMix) );
            sigmayy(:,:,iMix)=sigmayy(:,:,iMix)-diag( squeeze(D(iMix,iO,:)).*ary(:,iO,iMix)  );

     end
     sigma(:,:,iMix)=blkdiag(sigmaxx(:,:,iMix),sigmayy(:,:,iMix));
      [D_,p_]=chol(sigma(:,:,iMix));
      if p_>0
        sigma(:,:,iMix)=sigma(:,:,iMix)+ diag(vFloor);
      end
         var(:,iMix)=diag(sigma(:,:,iMix));
            gconst(iMix)=log(w(iMix))-l/2*log(2*pi)-0.5*log(prod( var(1:l/2,iMix) )) -0.5*log(prod(var(l/2+1:l,iMix)));
      % gconst2(iMix)=log(w(iMix))-l/2*log(2*pi)-0.5*log(prod( var(:,iMix) ));
    end
    ite=ite+1;
end

%%
ogmm=fopen([trainDir '/neutral_' emotion '.phrase.gmm_ar_ABfull'],'w');
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
    
    fprintf(ogmm,'%f\n',bt(:,:,iMix));
    fprintf(ogmm,'%f\n',sigma(:,:,iMix));
end
