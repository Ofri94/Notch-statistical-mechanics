%this script is meant to solve numerically the partition function and sensitivity
%function of double sites notch gene + an Fr site which can be modulated in strength.



%% both activator and repressor

a=[0, logspace(-3,0.7,1001)];%alpha vec - spanned logarithmically to improve memory
b=[0, (logspace(-3,0.7,1001))]';%beta vec - spanned logarithmically to improve memory
A=meshgrid(a);%creating a metrix of alpha in order to test for each beta

for Fr=[1 0.7 0]%indicating site functionality - 0being non functional and 1 fully functional.
for C=[1]%solving for non cooperative and cooperative
    
    z=1 + 3*A + (C+2)*A.^2 + C*A.^3 + 3*b + 3*b.^2 + b.^3 + 6*A.*b + (C+2)*b.*A.^2 + 3*A.*b.^2;%the partition function is the same as in 3 sites
    
    %"probabilities"
    p0=1./z;%probability of 0 sites occupied
    p1=(3*a)./z;%probability of 1 site occupied by ACTIVATOR
    p1_1_a=4*A.*b./z;%prob. of 1 rep and 1 activator, dist Fr
    p1_1_b=2*A.*b./z;%prob. of 1 rep and 1 activator, dist 2*Fr
    p1_2_a=A.*b.^2./z;% prob. of 2 rep and 1 act. activator in the middle of both rep
    p1_2_b=2*A.*b.^2./z;% prob. of 2 rep and 1 act. adjacent repressors
    p2=((C+2)*a.^2)./z;%probability of 2 sites occupied by ACTIVATOR
    p2_1_a=b.*A.^2./z;%probability of 2 sites occupied by ACTIVATOR and 1 repressor, repressor in the middle
    p2_1_b=2*b.*A.^2./z;%probability of 2 sites occupied by ACTIVATOR and 1 repressor, repressor on the side
    p3=C*a.^3./z;%probability of 3 sites occupied by ACTIVATOR
    %IMPORTANT NOTE - in p3 and p4 each COLUMN is a for a fixed alpha, as
    %opposed to p2 and p3, where I look at the ROWS for fixed betas
    %IMPORTANT NOTE 2 - p1_1 wasn't included in active transcrition since
    %repressor prescence is ASSUMED to inhibit transcrition entirely
    
    %We assume site occuancy to be directly related to transcriptional outcome
    %Hence T=Occ. I'll plot T as a function of alpha
    w=[1 2 3];
    %weight of activation for binding states included in "probabilities"
    %location 1 is the weight for single site bound, location 2 is for 2 sites et. cetra
    T=w(1)*p1+w(2)*p2+w(3)*p3+w(1)*p1_1_a*(1-Fr)+p1_1_b*(1-Fr^2)+w(1)*p1_2_a*(1-Fr)^2+w(1)*p1_2_b*(1-Fr)*(1-Fr^2)+w(2)*p2_1_a*(1-Fr)+p2_1_b*((1-Fr)+(1-Fr^2));%transcription is calculated by acticvator statistical density multiplied by number of activators in the state,
    %with the addition of activator/rep combinations in which Fr symbolizes distance between sites as to impair rep ability to act on activator (David's model)
    
    
    
    dlog_Out_a=diff(log(T),1,2);%differentiating log of output acros columns - rows are different betas
    dlog_In_a=diff(log(A),1,2);%differentiating log of input acros columns - rows are different betas
    Sa=dlog_Out_a./dlog_In_a;%calc. sensitivity for alpha
    
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

     [~, t2,~,~,~]=nth_csl_mth_sps(2,0);
    t=T./t2;
    t(813,813)
    mean(t(:,2:end),'all')
    
      plot_T(A,b,t, Fr)
%     plot_Sa(A,b,Sa, Fr)
%     plot_Sb(A,b,Sb, Fr)
%     plot_phase(A,b,phase, Fr)

end
end
%% plotting transcription

function [] = plot_T(A,b,T, Fr)
%plotting Transcription
figure('WindowState', 'maximized');
s=surf(A, meshgrid(b)', T);
view(2)
colorbar
s.EdgeColor = 'none';
title({'Transcription as a function of activator and repressor conc.', ['3 sites with distance Fr=' num2str(Fr)]});
xlabel('alpha');ylabel('beta');zlabel('Transcription')
caxis([0 2])
pbaspect([1 1 1])
xlim([0 6])
ylim([0 6])
end

%%%%%%%%%%%%%%%%%%%%%%
%% plotting sensitivity to a
function [] = plot_Sa(A,b,Sa, Fr)
%plotting sensitivity to a
figure('WindowState', 'maximized');
a=A(1,:);

s=surf(A(1:(end-1),1:(end-1)), meshgrid(b(1:(end-1)))' , Sa(1:(end-1),:));%S is calculated with diff, which shortens the array, so the (end-1) align the sizes
view(2)
colorbar
s.EdgeColor = 'none';
title({'Sensitivity for activator as a function of activator and repressor conc.', ['3 sites with distance Fr=' num2str(Fr)]});
xlabel('alpha');ylabel('beta');zlabel('Sensitivity')
caxis([0 3])
hold on

ind_a=Sa<1;
Pos=logical(stdfilt(ind_a));
[x, y]=ind2sub([size(Pos)],find(Pos));
plot3( b(y)', a(x), Sa(Pos),'k','LineWidth', 2)%marking sensitivity value of 1

hold off

pbaspect([1 1 1])
xlim([0 6])
ylim([0 6])
end
%%%%%%%%%%%%%%%%%%%%%
%% plotting sensitivity to b
function [] = plot_Sb(A,b,Sb, Fr)
%plotting sensitivity to b
figure('WindowState', 'maximized');
a=A(1,:);
s=surf(A(1:(end-1),1:(end-1)), meshgrid(b(1:(end-1)))' , Sb(:,1:(end-1)));%S is calculated with diff, which shortens the array, so the (end-1) align the sizes
view(2)
colorbar
s.EdgeColor = 'none';
title({'Sensitivity for reperssor as a function of activator and repressor conc.', ['3 sites with distance Fr=' num2str(Fr)]});
xlabel('alpha');ylabel('beta');zlabel('Sensitivity')
caxis([0 3])

hold on

ind_b=Sb<1;
Pos=logical(stdfilt(ind_b));
[x, y]=ind2sub([size(Pos)],find(Pos));
plot3( b(y)', a(x), Sb(Pos),'k','LineWidth', 2)%marking sensitivity value of 1

hold off

pbaspect([1 1 1])
xlim([0 6])
ylim([0 6])
end
%%%%%%%%%%%%%%%%%%%%%%
%% plotting phase plane
function [] = plot_phase(A,b,phase,Fr)
%plotting phase plane of sensitivity
figure('WindowState', 'maximized');
s=surf(A(1:(end-1),1:(end-1)), meshgrid(b(1:(end-1)))' , phase(1:(end-1),1:(end-1)));%S is calculated with diff, which shortens the array, so the (end-1) align the sizes
view(2)
s.EdgeColor = 'none';
title({'Sensitivity areas as a function of activator and repressor conc.', ['3 sites with distance Fr=' num2str(Fr)]});
xlabel('alpha');ylabel('beta');
pbaspect([1 1 1])
xlim([0 6])
ylim([0 6])
end
%%%%%%%%%%%%%%%%

