function [z,T,Sa,Sb,phase] = nth_csl_mth_sps(N,M)

%this script is meant to solve numerically the partition function and sensitivity
%function of nth csl sites and mth sps sites notch gene.
%alpha and beta ranges logarithmically from 0 to 5
%this code assumes repressor prescence completly inhibits transcription,
%therefore states that occupy b result in 0 transcription. even if there
%are huge amounts of a.


%  a=[0, logspace(-3,0.7,1001)];%alpha vec - spanned logarithmically to improve memory
%  b=[0, (logspace(-3,0.7,1001))]';%beta vec - spanned logarithmically to improve memory
 a=0:0.001:1;
 b=(0:0.001:1)';
A=meshgrid(a);%creating a metrix of alpha in order to test for each beta

C=20;%choosing cooperativity factor of 20

%first, calc. how many terms would be in the polynomial sum (1+2a+Ca^2)
num_of_sps_terms=nchoosek(M+2,2);%there are 3 "veriables" in the sum and M is the power of that sum - spss

arrangements=mnrnd(M,ones(1,3).*(1/3));%randomly creating a term

while size(arrangements,1)<num_of_sps_terms
    current=mnrnd(M,ones(1,3).*(1/3));%randomly craeting a term
    
    if sum(ismember(arrangements,current,'rows'))%cheking whether the term exists in the array
        continue
    else
        arrangements=[arrangements; current];%adding it to the array if it doesn't
    end
end

num_of_csl_terms=nchoosek(N+1,1);%there are 2 "veriables" in the sum and N is the power of that sum - csls


if N>=0 && M>0
    
    %combining into a partition function
    z=(1+A+b).^(N).*(1+2*A+(C*A.^2)+2*b+(b.^2)+2*A.*b).^(M);% where 1 descibes no occupency, a is the
    %statistical weight of activator binding to site, C is cooperativity
    %factor and b is statistical weight for repressor binding to site
    
    w=0:N;%weight for csl sites
    
    
    %in order to calculate all of the relevant probabilities for transcription
    %I will calculate the coefficient for the csl sites, than moltiply by 1,
    %2a, and Ca^2 to get all of the RELEVANT probabilities for active
    %transcription
    
    T=0;
%     arrangements(ismember(arrangements,[M, 0, 0],'rows'),:)=[];%erasing the row which results in 1^M. it doesnt need to be included in the calc.
    %note that the previous line mandate the use of 1 as the first out of
    %the three monomials (1+2a+Ca^2)
    
    for itr1=1:(num_of_sps_terms)%going through all of the terms from the sps part
        
        powers=arrangements(itr1,:);%the current term
        fa_powers=factorial(powers);
        current_term_co=(factorial(M)/(fa_powers(1)*fa_powers(2)*fa_powers(3)));%this is the formula for the term's coefficient
        current_term=current_term_co.*(1.^powers(1)).*((2*a).^powers(2)).*((C*a.^2).^powers(3));%formula for the appropriate term from the sps region of the script
        
        for itr2=1:num_of_csl_terms
            if (w(itr2)+powers(2)+2*powers(3))>=4%summing the power of the current alpha for both csl and sps occupation
                T=T + (4)*(((nchoosek(N,(w(itr2))).*a.^w(itr2)).*current_term)./z);%calculating term when its power is bigger or equal to than 4
            elseif (w(itr2)+powers(2)+2*powers(3))<4
                T=T + (w(itr2)+powers(2)+2*powers(3))*(((nchoosek(N,(w(itr2))).*a.^w(itr2)).*current_term)./z);%calculating term when its power is smaller than 4
            end
        end
    end
    
elseif N>=0 && M==0
    [z,T,~,~,~] = nth_csl_sites(N);
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
