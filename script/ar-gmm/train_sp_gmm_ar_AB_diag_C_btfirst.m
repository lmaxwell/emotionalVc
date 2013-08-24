clear;
nMix=64;
emotion='sad';
trainDir=['../../vc/train/to' emotion '_ar/gmm_jde/' num2str(nMix) '-mix'];
gmmMean=[trainDir '/neutral-' emotion num2str(nMix) '.mean.txt' ];
gmmCov=[trainDir '/neutral-' emotion num2str(nMix) '.cov.txt' ];
gmmWght=[trainDir '/neutral-' emotion num2str(nMix) '.wght.txt' ];

gmmFile=[trainDir '/../../gmm_para/'  num2str(nMix) '-mix' '/neutral-' emotion num2str(nMix) '_ordery0'  '_btzero' num2str(1) '.para.txt' ];
trainFile=[trainDir '/../../dtw/' num2str(nMix) '-mix/neutral-' emotion '.mat-ar.txt'];
l=48;
%nMix=64;
mixNum=64;
orderx=1;
ordery=1;
vfloor=0.001;
btzero=1;
order=1;
addpath('lightspeed');
logfd=fopen([trainDir '/log_xy_ABdiag_order' num2str(ordery) '_btfirst' num2str(btzero)],'w');
%% read gmm model
   miux=zeros(l/2,mixNum);
    miuy=zeros(l/2,mixNum);
    miu=zeros(l,mixNum);
    sigma=zeros(l,l,mixNum);
    var=zeros(l,mixNum);
     ar=zeros(l,ordery,mixNum);
     arx=zeros(l/2,ordery,mixNum);
     ary=zeros(l/2,ordery,mixNum);
     bt=zeros(l/2,mixNum);
     gconst=zeros(mixNum,1);
gmmP=importdata(gmmFile);

orderx=0;
ordery=0;
   w=gmmP(1:mixNum);
   
  for i=1:mixNum
        miux(:,i)=gmmP(mixNum+(i-1)*(l+l/2*(1+orderx) +l/2*(1+ordery)+l/2 )+1:mixNum+(i-1)*(l+ l/2*(1+orderx) +l/2*(1+ordery)+l/2 )+l/2);
        miuy(:,i)=gmmP(mixNum+(i-1)*(l+l/2*(1+orderx) +l/2*(1+ordery) +l/2)+l/2+1:mixNum+(i-1)*(l+l/2*(1+orderx) +l/2*(1+ordery)+l/2 )+l);
        for iO=1:orderx
            arx(:,iO,i)=gmmP(mixNum+(i-1)*(l+l/2*(1+orderx) +l/2*(1+ordery) +l/2 )+l + l/2*(iO-1)+1 :mixNum+(i-1)*(l+l/2*(1+orderx) +l/2*(1+ordery) +l/2 )+l+ l/2*(iO-1)+l/2);
        end
         for iO=1:ordery
            ary(:,iO,i)=gmmP(mixNum+(i-1)*(l+l/2*(1+orderx) +l/2*(1+ordery) +l/2 )+l+l/2*orderx + l/2*(iO-1)+ 1 :mixNum+(i-1)*(l+l/2*(1+orderx) +l/2*(1+ordery) +l/2 )+ l+l/2*orderx+ l/2*(iO-1)+l/2);
            ar(:,iO,i)=[arx(:,iO,i); ary(:,iO,i)];
         end
        
       bt(:,i)=gmmP(mixNum+(i-1)*(l+l/2*(1+orderx) +l/2*(1+ordery)+l/2 )+ l +l/2*orderx + l/2*ordery+  1 :mixNum+(i-1)*(l+l/2*(1+orderx) +l/2*(1+ordery)+l/2 )+ l + l/2*orderx+ l/2*ordery+l/2 );
        miu(:,i)=[miux(:,i);miuy(:,i)];
      %  beginModel.miu(:,i)=[miux(:,i);miuy(:,i)];
       % sigma(:,:,i)=reshape(gmmP(mixNum+(i-1)*(l^2+l/2*(1+orderx) +l/2*(1+ordery)+l/2 )+ l/2*(1+orderx) +l/2*(1+ordery)+l/2 + 1: mixNum+(i-1)*(l^2+l/2*(1+orderx) +l/2*(1+ordery)+l/2 )+  l/2*(1+orderx) +l/2*(1+ordery) + l/2+l^2),l,l);
       % beginModel.sigma(:,:,i)=sigma(:,:,i);
        var(:,i)=gmmP(mixNum+(i-1)*(l+l/2*(1+orderx) +l/2*(1+ordery)+l/2 )+ l/2*(1+orderx) +l/2*(1+ordery)+l/2 + 1: mixNum+(i-1)*(l+l/2*(1+orderx) +l/2*(1+ordery)+l/2 )+  l/2*(1+orderx) +l/2*(1+ordery) + l/2+l);
         gconst(i)=log(w(i))-l/2*log(2*pi)-0.5*log(prod( var(:,i) ));
  end
    
  orderx=1;ordery=1;
%% train ar combined gmm

%first run  
trainData=importdata([trainFile '-0' ])';
preData=importdata([trainFile '-' num2str(1)])';
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
% %
% w=Priors';
% miu=Mu;
% sigma=Sigma;
%  gconst=zeros(mixNum,1);
% invSigma=zeros(l,l,mixNum);
% var=zeros(l,mixNum);
% for iMix=1:mixNum
%     
%    var(:,iMix)=diag( sigma(:,:,iMix) );
%   gconst(iMix)=log(w(iMix))-l/2*log(2*pi)-0.5*log(prod( var(:,iMix) ));
%     
%     invSigma(:,:,iMix)=inv( sigma(:,:,iMix) );
% end 
clear sigma Sigma Mu Prios Data nbStates Centers Data_id;
%%
indexOfSen=find(trainData(1,:)==0);
numOfSen=length(indexOfSen);
indexOfSen=[indexOfSen  size(trainData,2)+1];
beginZero=zeros(l,1);

A=zeros(mixNum,l/2,orderx+1,orderx+1);
C=zeros(mixNum,l/2,orderx+1);
B=zeros(mixNum,l/2,ordery+2,ordery+2);
D=zeros(mixNum,l/2,ordery+2);


postProb=zeros(mixNum,1);

liklihoodProb=0;
numOfSeg=0;
accW=zeros(mixNum,1);
accSig=zeros(l,mixNum);
accSig1=zeros(l,mixNum);
accSig2=zeros(mixNum,l,ordery);

accSigy3=zeros(l/2,mixNum);
accMiu=zeros(l,mixNum);

vFloor=zeros(l,1);

for iSen=1:numOfSen
    tic;
    start=indexOfSen(iSen)+1;
    xend=indexOfSen(iSen+1)-1;
    fprintf( 1 , '%dth sentence,start:%d,end:%d\n' , iSen , start , xend);
    for iSeg=0:xend-start
    %    fprintf(1,'%dth segment\n',iSeg);
        %set preSegs
%         if iSeg==0
%               preSegs=beginZero*ones(1,ordery-iSeg);
%         elseif(iSeg<ordery)
%            
%             preSegs=[beginZero*ones(1,ordery-iSeg) trainData(:,start+iSeg-1:-1:start)];
%            
%         else
%             preSegs=trainData(:,start+iSeg-1:-1:start+iSeg-ordery);
%        end
%     if(all( trainData(1:l/2,start+iSeg)==trainData(1:l/2,start+iSeg-1)) || all( trainData(l/2+1:l,start+iSeg)==trainData(l/2+1:l,start+iSeg-1))  )
%         continue;
%     end
     preSegs(:,1)=preData(:,start+iSeg);

      %  sum_postProb=-1.0E10;
        maxProb=-1.0E10;
    %    maxMix=1;
%        for iMix=1:mixNum
%             
%               
%           postProb(iMix)= gconst(iMix)-sum(0.5.*(trainData(:,start+iSeg)-miu(:,iMix)).*(trainData(:,start+iSeg)-miu(:,iMix))./var(:,iMix) );
%         
%          
% %   
% %            sum_postProb =log_add( sum_postProb ,postProb(iMix));
% %  
% %             if maxProb<postProb(iMix)
% %                 maxMix=iMix;
% %                 maxProb=postProb(iMix);
% %             end
%        end

      curMiu=[miux ;miuy];
        for iO=1:ordery
            curMiu(1:l/2,:)=curMiu(1:l/2,:)+preSegs(1:l/2,iO)*ones(1,mixNum).*squeeze(arx(:,iO,:));
            curMiu(l/2+1:l,:)=curMiu(l/2+1:l,:) + preSegs(l/2+1:l,iO)*ones(1,mixNum).*squeeze( ary(:,iO,:) );
        end
        
        curMiu(l/2+1:l,:)=curMiu(l/2+1:l,:) + trainData(1:l/2,start+iSeg)*ones(1,mixNum).*bt;
        
        postProb=gconst- ndsum(0.5.*(trainData(:,start+iSeg)*ones(1,mixNum)-curMiu).*(trainData(:,start+iSeg)*ones(1,mixNum)-curMiu)./var, 1) ;
        maxMix=find( postProb==max(postProb));
       sum_postProb=logsumexp(postProb);
       postProb=exp(postProb-sum_postProb);
         accW=accW+postProb;
        liklihoodProb=liklihoodProb+sum_postProb;
        
         for iOd=1:ordery  
                for jOd=1:ordery 
                    A(:,:,iOd,jOd)=A(:,:,iOd,jOd)+postProb* ( preSegs(1:l/2,iOd).*preSegs(1:l/2,jOd) )'; 
                    B(:,:,iOd,jOd)=B(:,:,iOd,jOd)+postProb* ( preSegs(l/2+1:l,iOd).*preSegs(l/2+1:l,jOd) )'; 
                end
                C(:,:,iOd)=C(:,:,iOd)+postProb*( trainData(1:l/2,start+iSeg).*preSegs(1:l/2,iOd) )';
                D(:,:,iOd)=D(:,:,iOd)+postProb*( trainData(l/2+1:l,start+iSeg).*preSegs(l/2+1:l,iOd) )';
                
                 A(:,:,ordery+1,iOd)=A(:,:,ordery+1,iOd)+postProb*preSegs(1:l/2,iOd)' ;
                  A(:,:,iOd,ordery+1)=A(:,:,iOd,ordery+1)+ postProb*preSegs(1:l/2,iOd)'; 
                B(:,:,ordery+1,iOd)=B(:,:,ordery+1,iOd)+postProb*preSegs(l/2+1:l,iOd)' ;
                 B(:,:,ordery+2,iOd)=B(:,:,ordery+2,iOd)+postProb*( trainData(1:l/2,start+iSeg).*preSegs(l/2+1:l,iOd) )';
                  B(:,:,iOd,ordery+1)=B(:,:,iOd,ordery+1)+ postProb*preSegs(l/2+1:l,iOd)'; 
                   B(:,:,iOd,ordery+2)=B(:,:,iOd,ordery+2)+ postProb*( trainData(1:l/2,start+iSeg).*preSegs(l/2+1:l,iOd) )';  
                   
                    accSig2(:,:,iOd)=accSig2(:,:,iOd)+postProb*( trainData(1:l,start+iSeg).*preSegs(1:l,iOd) )';
                    %  accSig2(:,l/2+1:l,iOd)=accSig2(:,l/2+1:l,iOd)+postProb*( trainData(l/2+1:l,start+iSeg).*preSegs(l/2+1:l,iOd) )';
         end
         
      %    A(:,:,orderx+1,orderx+1)=   A(:,:,orderx+1,orderx+1)+postProb*ones(1,l/2);   
      %    C(:,:,orderx+1)=C(:,:,orderx+1)+postProb*trainData(1:l/2,start+iSeg)'; 
        
      
      % B(:,:,ordery+1,ordery+1)=   B(:,:,ordery+1,ordery+1)+postProb*ones(1,l/2);  
          
          % D(:,:,ordery+1)=D(:,:,ordery+1)+postProb*trainData(l/2+1:l,start+iSeg)'  
       %   B(:,:,ordery+1,ordery+2)=   B(:,:,ordery+1,ordery+2)+postProb*trainData(1:l/2,start+iSeg)';  
         %    B(:,:,ordery+2,ordery+1)=   B(:,:,ordery+2,ordery+1)+postProb*trainData(1:l/2,start+iSeg)';
         %  B(:,:,ordery+2,ordery+2)=   B(:,:,ordery+2,ordery+2)+postProb* ( trainData(1:l/2,start+iSeg).*trainData(1:l/2,start+iSeg) )';
         %     D(:,:,ordery+2)=D(:,:,ordery+2)+postProb*( trainData(l/2+1:l,start+iSeg).*trainData(1:l/2,start+iSeg) )'; 
             
                   accSig(:,:)=accSig(:,:)+trainData(1:l,start+iSeg).*trainData(1:l,start+iSeg)*postProb';
               
              %     accSig1(l/2+1:l,:)=accSig1(l/2+1:l,:)+trainData(l/2+1:l,start+iSeg)*postProb'; %=D(:,:,ordery+1)
 
         accSigy3(:,:)=accSigy3(:,:)+trainData(l/2+1:l,start+iSeg).*trainData(1:l/2,start+iSeg)*postProb' ; %=D(:,:,ordery+2)
              
         accMiu(:,:)=accMiu(:,:)+trainData(1:l,start+iSeg)*postProb' ;     %   =B(:,:,ordery+1,ordery+2) = B(:,:,ordery+2,ordery+1)
         
               

        
        
       
        
        vFloor=vFloor+(trainData(:,start+iSeg)-miu(:,maxMix)).* (trainData(:,start+iSeg)-miu(:,maxMix));
       
        numOfSeg=numOfSeg+1;
        
    end
 
    

    
    toc;
end
vFloor=vFloor*vfloor/numOfSeg;
C(:,:,ordery+1)=accMiu(1:l/2,:)';
D(:,:,ordery+1)=accMiu(l/2+1:l,:)';
A(:,:,ordery+1,ordery+1)=accW*ones(1,l/2);
B(:,:,ordery+1,ordery+1)=accW*ones(1,l/2);
B(:,:,ordery+1,ordery+2)=  accMiu(1:l/2,:)' ;
B(:,:,ordery+2,ordery+1)=  accMiu(1:l/2,:)';
  B(:,:,ordery+2,ordery+2)=accSig(1:l/2,:)';
   D(:,:,ordery+2)=accSigy3';
      accSig1(1:l/2,:)=accMiu(1:l/2,:);
   accSig1(l/2+1:l,:)=accMiu(l/2+1:l,:);
%%
    for iMix=1:mixNum
               for iL=1:l/2
              tempA=squeeze(A(iMix,iL,:,:));
              tempC=squeeze(C(iMix,iL,:,:));
            if min(svd(tempA))<1e-7
                temp=pinv(tempA)*tempC;
            else
                temp=tempA \tempC;
            end
           
            arx(iL,:,iMix)=temp(1:orderx);
            miux(iL,iMix)=temp(orderx+1);
            
            eig=roots([1 -arx(iL,:,iMix)]);
           
           while(find(abs(eig)>1))
             ii=find(abs(eig)>1);
             fprintf(1,'%d mix ; %d dimension ;arx=%f\n',iMix,ii,eig(ii));
             eig(ii)=eig(ii)/abs(eig(ii));
             
           end
           tpoly=poly(eig);
           arx(iL,:,iMix)=-tpoly(2:orderx+1);
            
        end
       % miux(:,iMix)=accMiu(1:l/2,iMix)./accW(iMix);
           
        for iL=1:l/2
            tempB=squeeze( B(iMix,iL,:,:) );
          tempD=squeeze( D(iMix,iL,:) );
             if min(svd(tempB))<1e-7
                temp=pinv(tempB)* tempD ;
            else
                temp=tempB \ tempD;
             end
            ary(iL,:,iMix)=temp(1:ordery);
            miuy(iL,iMix)=temp(ordery+1);
         
       %     bt(iL,iMix)=temp(ordery+2);
            eig=roots([1 -ary(iL,:,iMix)]);
           
           while(find(abs(eig)>1))
             ii=find(abs(eig)>1);
                  fprintf(1,'%d mix ; %d dimension ;ary=%f\n',iMix,ii,eig(ii));
             eig(ii)=eig(ii)/abs(eig(ii));
         
           end
           tpoly=poly(eig);
           ary(iL,:,iMix)=-tpoly(2:ordery+1);
            %bt(iL,iMix)=0;
        end
    
       w(iMix)=accW(iMix)/numOfSeg;
         %sigmaxx(:,:,iMix)=diag(accSig(1:l/2,iMix)./accW(iMix)-miux(:,iMix).*miux(:,iMix));
         var(1:l/2,iMix)=(accSig(1:l/2,iMix)-accSig1(1:l/2,iMix).*miux(:,iMix))./accW(iMix);
        %  sigmayy(:,:,iMix)=diag(accSig(l/2+1:l,iMix)-accSig1(l/2+1:l,iMix).*miuy(:,iMix)  -accSigy3(:,iMix).*bt(:,iMix)        )./accW(iMix);
          var(l/2+1:l,iMix)= ( accSig(l/2+1:l,iMix)-accSig1(l/2+1:l,iMix).*miuy(:,iMix)  -accSigy3(:,iMix).*bt(:,iMix)        )./accW(iMix) ;
     for iO=1:ordery  
           var(1:l/2,iMix)=var(1:l/2,iMix)-accSig2(iMix,1:l/2,iO)'.*arx(:,iO,iMix)./accW(iMix);
          %  sigmaxx(:,:,iMix)=sigmaxx(:,:,iMix)-diag(accSig2(1:l/2,iO,iMix).*arx(:,iO,iMix) )./accW(iMix);
           % sigmayy(:,:,iMix)=sigmayy(:,:,iMix)-diag( accSig2(iMix,l/2+1:l,iO)'.*ary(:,iO,iMix)  )./accW(iMix);
            var(l/2+1:l,iMix)=var(l/2+1:l,iMix)-accSig2(iMix,l/2+1:l,iO)'.*ary(:,iO,iMix) ./accW(iMix);
     end
   %  sigma(:,:,iMix)=blkdiag(sigmaxx(:,:,iMix),sigmayy(:,:,iMix));
%       [D_,p_]=chol(sigma(:,:,iMix));
%       if p_>0
%         sigma(:,:,iMix)=sigma(:,:,iMix)+ diag(vFloor);
%       end
      
        gconst(iMix)=log(w(iMix))-l/2*log(2*pi)-0.5*log(prod( var(:,iMix) ));
  %  invSigma(:,:,iMix)=inv( sigma(:,:,iMix) );
      
      
    end

% ogmm=fopen([trainDir '/gmm_ar_diag_1thit_order' num2str(ordery) ],'w');
% for iMix=1:mixNum
%     fprintf(ogmm,'%f\n',w(iMix)); 
% end
% for iMix=1:mixNum
%     fprintf(ogmm,'%f\n',miu(:,iMix));
%     for iO=1:ordery
%         fprintf(ogmm,'%f\n',ar(:,ordery,iMix));
%     end
%     fprintf(ogmm,'%f\n',sigma(:,:,iMix));
% end
%%
;
preProb=liklihoodProb/numOfSeg;
curProb=0;
 fprintf(1,'initialProb %f\n',preProb);
  fprintf(logfd,'initialProb %f\n',preProb);
%iteration
%%
ite=1;
while(abs((curProb-preProb)/preProb)>1e-5 && ite<20)
A=zeros(mixNum,l/2,orderx+1,orderx+1);
C=zeros(mixNum,l/2,orderx+1);
B=zeros(mixNum,l/2,ordery+2,ordery+2);
D=zeros(mixNum,l/2,ordery+2);


postProb=zeros(mixNum,1);

liklihoodProb=0;
numOfSeg=0;
accW=zeros(mixNum,1);
accSig=zeros(l,mixNum);
accSig1=zeros(l,mixNum);
accSig2=zeros(mixNum,l,ordery);

accSigy3=zeros(l/2,mixNum);
accMiu=zeros(l,mixNum);




for iSen=1:numOfSen
    tic;
    start=indexOfSen(iSen)+1;
    xend=indexOfSen(iSen+1)-1;
  %  fprintf( 1 , '%dth sentence,start:%d,end:%d\n' , iSen , start , xend);
    for iSeg=0:xend-start
 %       fprintf(1,'%dth segment\n',iSeg);
        %set preSegs
%          if iSeg==0
%               preSegs=beginZero*ones(1,ordery-iSeg);
%         elseif(iSeg<ordery)
%            
%             preSegs=[beginZero*ones(1,ordery-iSeg) trainData(:,start+iSeg-1:-1:start)];
%            
%         else
%             preSegs=trainData(:,start+iSeg-1:-1:start+iSeg-ordery);
%         end
%   if(all( trainData(1:l/2,start+iSeg)==trainData(1:l/2,start+iSeg-1)) || all( trainData(l/2+1:l,start+iSeg)==trainData(l/2+1:l,start+iSeg-1))  )
%         continue;
%     end
        preSegs(:,1)=preData(:,start+iSeg);
        sum_postProb=-1.0E10;
       
%         for iMix=1:mixNum
% 
%            
%                 curMiu=[miux(:,iMix);miuy(:,iMix)];
%                 for iO=1:ordery
%                 %    curMiu(1:l/2)=curMiu(1:l/2)+arx(:,iO,iMix).*preSegs(1:l/2,iO);
%                     curMiu(l/2+1:l)=curMiu(l/2+1:l)+ary(:,iO,iMix).*preSegs(l/2+1:l,iO);
%                   
%                 end
%                    curMiu(l/2+1:l)=curMiu(l/2+1:l)+bt(:,iMix).*trainData(1:l/2,start+iSeg);
%     
%                   postProb(iMix)= gconst(iMix)-sum(0.5.*(trainData(:,start+iSeg)-curMiu).*(trainData(:,start+iSeg)-curMiu)./var(:,iMix) );
%      
%              
%         end
        
        curMiu=[miux ;miuy];
        for iO=1:ordery
            curMiu(1:l/2,:)=curMiu(1:l/2,:)+preSegs(1:l/2,iO)*ones(1,mixNum).*squeeze(arx(:,iO,:));
            curMiu(l/2+1:l,:)=curMiu(l/2+1:l,:) + preSegs(l/2+1:l,iO)*ones(1,mixNum).*squeeze( ary(:,iO,:) );
        end
        
        curMiu(l/2+1:l,:)=curMiu(l/2+1:l,:) + trainData(1:l/2,start+iSeg)*ones(1,mixNum).*bt;
        
        postProb=gconst- ndsum(0.5.*(trainData(:,start+iSeg)*ones(1,mixNum)-curMiu).*(trainData(:,start+iSeg)*ones(1,mixNum)-curMiu)./var, 1) ;
         sum_postProb=logsumexp(postProb);
        postProb=exp(postProb-sum_postProb);
        accW=accW+postProb;
        liklihoodProb=liklihoodProb+sum_postProb;
        
          for iOd=1:ordery  
                for jOd=1:ordery 
                    A(:,:,iOd,jOd)=A(:,:,iOd,jOd)+postProb* ( preSegs(1:l/2,iOd).*preSegs(1:l/2,jOd) )'; 
                    B(:,:,iOd,jOd)=B(:,:,iOd,jOd)+postProb* ( preSegs(l/2+1:l,iOd).*preSegs(l/2+1:l,jOd) )'; 
                end
                C(:,:,iOd)=C(:,:,iOd)+postProb*( trainData(1:l/2,start+iSeg).*preSegs(1:l/2,iOd) )';
                D(:,:,iOd)=D(:,:,iOd)+postProb*( trainData(l/2+1:l,start+iSeg).*preSegs(l/2+1:l,iOd) )';
                
                 A(:,:,ordery+1,iOd)=A(:,:,ordery+1,iOd)+postProb*preSegs(1:l/2,iOd)' ;
                  A(:,:,iOd,ordery+1)=A(:,:,iOd,ordery+1)+ postProb*preSegs(1:l/2,iOd)'; 
                B(:,:,ordery+1,iOd)=B(:,:,ordery+1,iOd)+postProb*preSegs(l/2+1:l,iOd)' ;
                 B(:,:,ordery+2,iOd)=B(:,:,ordery+2,iOd)+postProb*( trainData(1:l/2,start+iSeg).*preSegs(l/2+1:l,iOd) )';
                  B(:,:,iOd,ordery+1)=B(:,:,iOd,ordery+1)+ postProb*preSegs(l/2+1:l,iOd)'; 
                   B(:,:,iOd,ordery+2)=B(:,:,iOd,ordery+2)+ postProb*( trainData(1:l/2,start+iSeg).*preSegs(l/2+1:l,iOd) )';  
                   
                    accSig2(:,:,iOd)=accSig2(:,:,iOd)+postProb*( trainData(1:l,start+iSeg).*preSegs(1:l,iOd) )';
                    %  accSig2(:,l/2+1:l,iOd)=accSig2(:,l/2+1:l,iOd)+postProb*( trainData(l/2+1:l,start+iSeg).*preSegs(l/2+1:l,iOd) )';
         end
         
      %    A(:,:,orderx+1,orderx+1)=   A(:,:,orderx+1,orderx+1)+postProb*ones(1,l/2);   
      %    C(:,:,orderx+1)=C(:,:,orderx+1)+postProb*trainData(1:l/2,start+iSeg)'; 
        
      
      % B(:,:,ordery+1,ordery+1)=   B(:,:,ordery+1,ordery+1)+postProb*ones(1,l/2);  
          
          % D(:,:,ordery+1)=D(:,:,ordery+1)+postProb*trainData(l/2+1:l,start+iSeg)'  
       %   B(:,:,ordery+1,ordery+2)=   B(:,:,ordery+1,ordery+2)+postProb*trainData(1:l/2,start+iSeg)';  
         %    B(:,:,ordery+2,ordery+1)=   B(:,:,ordery+2,ordery+1)+postProb*trainData(1:l/2,start+iSeg)';
         %  B(:,:,ordery+2,ordery+2)=   B(:,:,ordery+2,ordery+2)+postProb* ( trainData(1:l/2,start+iSeg).*trainData(1:l/2,start+iSeg) )';
         %     D(:,:,ordery+2)=D(:,:,ordery+2)+postProb*( trainData(l/2+1:l,start+iSeg).*trainData(1:l/2,start+iSeg) )'; 
             
                   accSig(:,:)=accSig(:,:)+trainData(1:l,start+iSeg).*trainData(1:l,start+iSeg)*postProb';
               
              %     accSig1(l/2+1:l,:)=accSig1(l/2+1:l,:)+trainData(l/2+1:l,start+iSeg)*postProb'; %=D(:,:,ordery+1)
 
         accSigy3(:,:)=accSigy3(:,:)+trainData(l/2+1:l,start+iSeg).*trainData(1:l/2,start+iSeg)*postProb' ; %=D(:,:,ordery+2)
              
         accMiu(:,:)=accMiu(:,:)+trainData(1:l,start+iSeg)*postProb' ;     %   =B(:,:,ordery+1,ordery+2) = B(:,:,ordery+2,ordery+1)
         
         
         
       
        numOfSeg=numOfSeg+1;
    end
    
   
    
    toc;
end
preProb=curProb;

 curProb=liklihoodProb/numOfSeg;

 fprintf(1,'curProb %f\n',curProb);
  fprintf(logfd,'curProb %f\n',curProb);

C(:,:,ordery+1)=accMiu(1:l/2,:)';
D(:,:,ordery+1)=accMiu(l/2+1:l,:)';
A(:,:,ordery+1,ordery+1)=accW*ones(1,l/2);
B(:,:,ordery+1,ordery+1)=accW*ones(1,l/2);
B(:,:,ordery+1,ordery+2)=  accMiu(1:l/2,:)' ;
B(:,:,ordery+2,ordery+1)=  accMiu(1:l/2,:)';
  B(:,:,ordery+2,ordery+2)=accSig(1:l/2,:)';
   D(:,:,ordery+2)=accSigy3';
      accSig1(1:l/2,:)=accMiu(1:l/2,:);
   accSig1(l/2+1:l,:)=accMiu(l/2+1:l,:);
  
  for iMix=1:mixNum
        for iL=1:l/2
              tempA=squeeze(A(iMix,iL,:,:));
              tempC=squeeze(C(iMix,iL,:,:));
            if min(svd(tempA))<1e-7
                temp=pinv(tempA)*tempC;
            else
                temp=tempA \tempC;
            end
           
            arx(iL,:,iMix)=temp(1:orderx);
            miux(iL,iMix)=temp(orderx+1);
            
             eig=roots([1 -arx(iL,:,iMix)]);
           
           while(find(abs(eig)>1))
             ii=find(abs(eig)>1);
             eig(ii)=eig(ii)/abs(eig(ii));
             
           end
           tpoly=poly(eig);
           arx(iL,:,iMix)=-tpoly(2:orderx+1);
            
        end
       % miux(:,iMix)=accMiu(1:l/2,iMix)./accW(iMix);
           
        for iL=1:l/2
            tempB=squeeze( B(iMix,iL,:,:) );
          tempD=squeeze( D(iMix,iL,:) );
             if min(svd(tempB))<1e-7
                temp=pinv(tempB)* tempD ;
            else
                temp=tempB \ tempD;
             end
            ary(iL,:,iMix)=temp(1:ordery);
            miuy(iL,iMix)=temp(ordery+1);
         
       %     bt(iL,iMix)=temp(ordery+2);
            eig=roots([1 -ary(iL,:,iMix)]);
           
           while(find(abs(eig)>1))
             ii=find(abs(eig)>1);
             eig(ii)=eig(ii)/abs(eig(ii));
             
           end
           tpoly=poly(eig);
           ary(iL,:,iMix)=-tpoly(2:ordery+1);
            %bt(iL,iMix)=0;
        end
       
        
       w(iMix)=accW(iMix)/numOfSeg;
         var(1:l/2,iMix)=(accSig(1:l/2,iMix)-accSig1(1:l/2,iMix).*miux(:,iMix))./accW(iMix);
        %  sigmayy(:,:,iMix)=diag(accSig(l/2+1:l,iMix)-accSig1(l/2+1:l,iMix).*miuy(:,iMix)  -accSigy3(:,iMix).*bt(:,iMix)        )./accW(iMix);
          var(l/2+1:l,iMix)= ( accSig(l/2+1:l,iMix)-accSig1(l/2+1:l,iMix).*miuy(:,iMix)  -accSigy3(:,iMix).*bt(:,iMix)        )./accW(iMix) ;
     for iO=1:ordery  
           var(1:l/2,iMix)=var(1:l/2,iMix)-accSig2(iMix,1:l/2,iO)'.*arx(:,iO,iMix)./accW(iMix);
          %  sigmaxx(:,:,iMix)=sigmaxx(:,:,iMix)-diag(accSig2(1:l/2,iO,iMix).*arx(:,iO,iMix) )./accW(iMix);
           % sigmayy(:,:,iMix)=sigmayy(:,:,iMix)-diag( accSig2(iMix,l/2+1:l,iO)'.*ary(:,iO,iMix)  )./accW(iMix);
            var(l/2+1:l,iMix)=var(l/2+1:l,iMix)-accSig2(iMix,l/2+1:l,iO)'.*ary(:,iO,iMix) ./accW(iMix);
     end
   %  sigma(:,:,iMix)=blkdiag(sigmaxx(:,:,iMix),sigmayy(:,:,iMix));
%       [D_,p_]=chol(sigma(:,:,iMix));
%       if p_>0
%         sigma(:,:,iMix)=sigma(:,:,iMix)+ diag(vFloor);
%       end
      
        gconst(iMix)=log(w(iMix))-l/2*log(2*pi)-0.5*log(prod( var(:,iMix) ));
  %  invSigma(:,:,iMix)=inv( sigma(:,:,iMix) );
   end
    %
    ite=ite+1;
    ;
end

%%

% ogmm=fopen([trainDir '/../../gmm_para/'  num2str(nMix) '-mix' '/neutral-' emotion num2str(nMix) '_ordery' num2str(ordery)  '_btzero' num2str(btzero) '.para.txt' ],'w');
% for iMix=1:mixNum
%     fprintf(ogmm,'%f\n',w(iMix));
% end
% for iMix=1:mixNum
%     fprintf(ogmm,'%f\n',miux(:,iMix));
%       fprintf(ogmm,'%f\n',miuy(:,iMix));
%     for iO=1:ordery
%         fprintf(ogmm,'%f\n',arx(:,iO,iMix));
%         fprintf(ogmm,'%f\n',ary(:,iO,iMix));
%     end
%     fprintf(ogmm,'%f\n',bt(:,iMix));
%     fprintf(ogmm,'%f\n',var(:,iMix));
% end
% fclose(ogmm);