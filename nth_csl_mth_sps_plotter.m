

T=[];
Sa=[];
Sb=[];
Phase=[];
N=[12 0];%csl
M=[0 6];%sps

for itr=1:length(N)
    
    [z,t,sa,sb,phase] = nth_csl_mth_sps(N(itr),M(itr));
    T=cat(3,T,t);
    Sa=cat(3,Sa,sa);
    Sb=cat(3,Sb,sb);
    Phase=cat(3,Phase,phase);
    % drawnow limitrate


% % implay(S)

a=[0, logspace(-3,0.7,1001)];%alpha vec - spanned logarithmically to improve memory
b=[0, (logspace(-3,0.7,1001))]';%beta vec - spanned logarithmically to improve memory
% a=0:0.01:5;
% b=(0:0.01:5)';
A=meshgrid(a);%creating a metrix of alpha in order to test for each beta

% figure('WindowState', 'maximized');

    
        plot_T(A,b,t,N(itr),M(itr))
%         saveas(gcf,['Transcription vs a and b ' num2str(N(itr)) ' csl, ' num2str(M(itr)) ' sps C=20'],'tif')
%         close;
   
        plot_Sa(A,b,sa,N(itr),M(itr))
%         saveas(gcf,['Sensitivity for activator vs a and b ' num2str(N(itr)) ' csl, ' num2str(M(itr)) ' sps C=20'],'tif')
%         close;
%     
        plot_Sb(A,b,sb,N(itr),M(itr))
%         saveas(gcf,['Sensitivity for repressor vs a and b ' num2str(N(itr)) ' csl, ' num2str(M(itr)) ' sps C=20'],'tif')
%         close;
    
        plot_phase(A,b,phase,N(itr),M(itr))
    
%     hold on
%     
%     ind_T1=double(t>0.1);
%     Pos1=logical(stdfilt(ind_T1));
%     [x1, y1]=ind2sub([size(Pos1)],find(Pos1));
%     plot3( b(y1)', a(x1), phase(Pos1),'k','LineWidth', 2)%marking sensitivity value of 1
%     
%     ind_T2=double(t<2);
%     Pos2=logical(stdfilt(ind_T2));
%     [x2, y2]=ind2sub([size(Pos2)],find(Pos2));
%     plot3( b(y2)', a(x2), phase(Pos2),'k','LineWidth', 2)%marking sensitivity value of 1
%     
% %     ind_T=ind_T1.*ind_T2.*4;%this creates an indication of locations where transcriptional value is within a certain range decided in the rows above. 4 is the max z-axis value
% %     
% %     s=surf(A(1:(end-1),1:(end-1)), meshgrid(b(1:(end-1)))' , ind_T(1:end-1,1:end-1));
% %     view(2)
% %     colorbar
% %     s.EdgeColor = 'none';
% %     s.FaceAlpha=0.5;
% %     s.FaceColor='interp';
% %     caxis([0 4])
% %     pbaspect([1 1 1])
% %     xlim([0 6])
% %     ylim([0 6])
%     
%     hold off
    
%         saveas(gcf,['sensitivity phases vs a and b ' num2str(N(itr)) ' csl, ' num2str(M(itr)) ' sps C=20'],'tif')
%         close;
    %     pause(1)
    end



%% plotting transcription

function [] = plot_T(A,b,T,N,M)
%plotting Transcription
figure('WindowState', 'maximized');
s=surf(A, meshgrid(b)', T);
view(2)
colorbar
s.EdgeColor = 'none';
title({'Transcription as a function of activator and repressor conc.', [num2str(N) ' csl, ' num2str(M) ' sps'], 'C=20'});
xlabel('alpha');ylabel('beta');zlabel('Transcription')
caxis([0 2])
pbaspect([1 1 1])
xlim([0 6])
ylim([0 6])
end

%%%%%%%%%%%%%%%%%%%%%%
%% plotting sensitivity to a
function [] = plot_Sa(A,b,Sa,N,M)
%plotting sensitivity to a
figure('WindowState', 'maximized');
a=A(1,:);

s=surf(A(1:(end-1),1:(end-1)), meshgrid(b(1:(end-1)))' , Sa(1:(end-1),:));%S is calculated with diff, which shortens the array, so the (end-1) align the sizes
view(2)
colorbar
s.EdgeColor = 'none';
title({'Sensitivity for activator as a function of activator and repressor conc.', [num2str(N) ' csl, ' num2str(M) ' sps'], 'C=20'});
xlabel('alpha');ylabel('beta');zlabel('Sensitivity')
caxis([0 4])
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
function [] = plot_Sb(A,b,Sb,N,M)
%plotting sensitivity to b
figure('WindowState', 'maximized');
a=A(1,:);
s=surf(A(1:(end-1),1:(end-1)), meshgrid(b(1:(end-1)))' , Sb(:,1:(end-1)));%S is calculated with diff, which shortens the array, so the (end-1) align the sizes
view(2)
colorbar
s.EdgeColor = 'none';
title({'Sensitivity for reperssor as a function of activator and repressor conc.', [num2str(N) ' csl, ' num2str(M) ' sps'], 'C=20'});
xlabel('alpha');ylabel('beta');zlabel('Sensitivity')
caxis([0 4])

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
function [] = plot_phase(A,b,phase,N,M)
%plotting phase plane of sensitivity
figure('WindowState', 'maximized');
s=surf(A(1:(end-1),1:(end-1)), meshgrid(b(1:(end-1)))' , phase(1:(end-1),1:(end-1)));%S is calculated with diff, which shortens the array, so the (end-1) align the sizes
view(2)
s.EdgeColor = 'none';
title({'Sensitivity areas as a function of activator and repressor conc.', [num2str(N) ' csl, ' num2str(M) ' sps'], 'C=20'});
xlabel('alpha');ylabel('beta');
pbaspect([1 1 1])
xlim([0 6])
ylim([0 6])
end
%%%%%%%%%%%%%%%%
