function sol=RMSE(MS,F)
MS=double(MS);
F=double(F);
[n,m,d]=size(F);
D=(MS(:,:,1:d)-F).^2;
sol=sqrt(sum(sum(sum(D)))/(n*m*d));
