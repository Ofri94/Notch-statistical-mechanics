function [z,T,Sa,Sb,phase] = nth_csl_sites(N)

%this script is meant to solve numerically the partition function and sensitivity
%function of nth sites notch gene without cooperativity and with a repressor.
%alpha and beta ranges logarithmically from 0 to 10



% a=[0, logspace(-3,0.7,1001)];%alpha vec - spanned logarithmically to improve memory
% b=[0, (logspace(-3,0.7,1001))]';%beta vec - spanned logarithmically to improve memory
 a=0:0.001:1;
 b=(0:0.001:1)';
A=meshgrid(a);%creating a metrix of alpha in order to test for each beta

%combining into a partition function
z=(1+A+b).^(N);% where 1 descibes no occupency, a is the
%statistical weight of activator binding to site, C is cooperativity
%factor and b is statistical weight for repressor binding to site
w=N:-1:1;%weight of activation for binding states corresponding to location.

T=0;
for n=1:N
    if w(n)<4
    T=T + w(n)*((nchoosek(N,(n-1)).*a.^w(n))./z);%calculating each alpha
    elseif w(n)>=4
      T=T + 4*((nchoosek(N,(n-1)).*a.^w(n))./z);%calculating each alpha  
end




dlog_Out_a=diff(log(T),1,2);%differentiating log of output acros columns - rows are different betas
dlog_In_a=diff(log(A),1,2);%differentiating log of input acros columns - rows are different betas

Sa=dlog_Out_a./dlog_In_a;%calc. sensitivity

dlog_Out_b=diff(log(T));%differentiating log of output acros columns - rows are different betas
dlog_In_b=diff(meshgrid(b)');%differentiating log of input acros rows - columns are different betas
Sb=abs(dlog_Out_b./dlog_In_b);%calc. sensitivity for beta

ind_a=Sa<1;

ind_b=Sb<1;

phase=zeros([size(T)]);
ind_a(:,size((T),2))=0;
ind_b(size((T),1),:)=0;
ind_cor=logical(~ind_a.*~ind_b);
phase(~ind_a)=1;
phase(~ind_b)=2;
phase(ind_cor)=3;

end

