
emotion=strvcat('neutral','sad','happy','angry');
prjDir='..';

workDir=[prjDir '/test/error_dctnum'];

figure;
taxis=[4:16];
%rmse=zeros(9,1);
for iE=1:4
for dctNum=4:16
    
    rmse(dctNum-3)=importdata([workDir '/'  deblank(emotion(iE,:))  '/'  num2str(dctNum)  '/reconstructRmse']);
end
if(iE==1)
    s='*k';
elseif(iE==2)
        s='ok';
elseif(iE==3)
    s='^k';
else
    s='+k';
end
plot(taxis,rmse,s);
hold on;
axis([4 16 0 9]);
end
legend('neutral','sad','happy','angry');
ylabel('HZ');

