clear
clc
tic
clf
format shortE %�޸���������
addpath(genpath('D:\������\������\�ڰ���\����\AFEM_2D_rectangular\basic_subroutine'));%���ļ����е���س���
addpath(genpath('D:\������\������\�ڰ���\����\AFEM_2D_rectangular\other_subroutine'));
addpath(genpath('D:\������\������\�ڰ���\����\AFEM_2D_rectangular\treat_boudary_subroutine'));
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��Ҫ����Ĳ���%%%%%%%%%%%%%%%%%%%%%%%%%
%% �������(�������±߽�)
left_i=0;right_i=1;bottom_i=0;top_i=1;%����
left=left_i;right=right_i;bottom=bottom_i;top=top_i;%����
%% ʱ���������
theta=1/2;t_initial=0;t_end=2;TT=2*3650;
%% ���������ͣ����ԺͶ���
basis_type_trial_s=241;basis_type_test_s=241;
%% ��˹���ֲ���
Gauss_nodes_number_2D=9;%���ڼ�����ػ���
Gauss_nodes_number_1D=4;%���ڼ����߻���
%% ��ɢ����
n_element=6;
uc=0.8;
dx=1/2^(n_element);dy=1/2^(n_element);
n_x=(right-left)/dx;%x��������Ԫ��
n_y=(top-bottom)/dy;%x��������Ԫ��
t_deta=dx/2;%ʱ���������,
Mm=round((t_end-t_initial)/t_deta);%ʱ���������
n_T=round(TT/Mm);%ʱ��������
u_b=1;0.308;%�߽�Ũ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
n=1;%�����ֵ�Ԫ��Ĵ���
m=0;%ʱ���������
% %% colormap�Ĵ���
% % nr=Mm;%�޸�������Ͷȼ��
% nrgb=Mm; %������������    
% for i=nrgb:10*nrgb
% rgb=64/i;
% mycolormap_r=interp1([1 8 24 40 56 64],[0 0 0 1 1  1 ],1:rgb:64);
% mycolormap_g=interp1([1 8 24 40 56 64],[0 0 1 1 0  0  ],1:rgb:64);
% mycolormap_b=interp1([1 8 24 40 56 64],[0 1 1 0 0  0  ],1:rgb:64);
% if length(mycolormap_r)>=nrgb
% break
% end
% end
% mycolor=[mycolormap_r',mycolormap_g',mycolormap_b'];  
% colormap(mycolor);

while m<Mm-1% ʱ�����
      
     if m==0
        %% �������꼰����Ԫ�ֲ���ź�������ӳ���ϵ����ѡ���Ԫ���޹�,���ѡ��241��
        [P,T]=Generate_P_T_rectangular(241,left,right,bottom,top,n_x,n_y);  %�����㼰ӳ�����
        %% ���ݻ��������Ͷ�����������Ԫ������꼰��Ԫ�ֲ���ź�������ӳ���ϵ����Ԫ���йأ�
        [Nlb_trial_s,Pb_trial_s,Tb_trial_s,N_x_FE_nodes_s,N_y_FE_nodes_s,]=Generate_Pb_Tb_rectangular(basis_type_trial_s,left,right,bottom,top,n_x,n_y);
        [Nlb_test_s,Pb_test_s,Tb_test_s,N_x_FE_nodes_s,N_y_FE_nodes_s]=Generate_Pb_Tb_rectangular(basis_type_test_s,left,right,bottom,top,n_x,n_y);
        %% ����Ԫ��Ϣ
        dx_elment_nodes_s=(right-left)/(N_x_FE_nodes_s-1);%�������ĵ�һ��������x������������Ԫ�ڵ�֮��ľ��룬N_x_FE_nodes_uΪx��������Ԫ�ڵ���
        dy_elment_nodes_s=(top-bottom)/(N_y_FE_nodes_s-1);%�������ĵڶ���������y������������Ԫ�ڵ�֮��ľ���
        XX_s=left:dx_elment_nodes_s:right;%x��������Ԫ�̶�
        YY_s=bottom:dy_elment_nodes_s:top;
        %% ������
          M  = assemble_matrix_2D_s('function_mass_coef',                  Gauss_nodes_number_2D,P,T,Pb_trial_s,Tb_trial_s,Pb_test_s,Tb_test_s,Nlb_trial_s,Nlb_test_s,basis_type_trial_s,0,0,basis_type_test_s,0,0);        
        A_D1 = assemble_matrix_2D_s_t('function_diffusion_t',              Gauss_nodes_number_2D,P,T,Pb_trial_s,Tb_trial_s,Pb_test_s,Tb_test_s,Nlb_trial_s,Nlb_test_s,basis_type_trial_s,1,0,basis_type_test_s,1,0,0);
        A_D2 = assemble_matrix_2D_s_t('function_diffusion_t',              Gauss_nodes_number_2D,P,T,Pb_trial_s,Tb_trial_s,Pb_test_s,Tb_test_s,Nlb_trial_s,Nlb_test_s,basis_type_trial_s,0,1,basis_type_test_s,0,1,0);
        A_C1 = assemble_matrix_2D_s_t('function_convection_x_t',           Gauss_nodes_number_2D,P,T,Pb_trial_s,Tb_trial_s,Pb_test_s,Tb_test_s,Nlb_trial_s,Nlb_test_s,basis_type_trial_s,1,0,basis_type_test_s,0,0,0);
        A_C2 = assemble_matrix_2D_s_t('function_convection_y_t',           Gauss_nodes_number_2D,P,T,Pb_trial_s,Tb_trial_s,Pb_test_s,Tb_test_s,Nlb_trial_s,Nlb_test_s,basis_type_trial_s,0,1,basis_type_test_s,0,0,0);
        A_R  = assemble_matrix_2D_s_t('function_reaction_t',               Gauss_nodes_number_2D,P,T,Pb_trial_s,Tb_trial_s,Pb_test_s,Tb_test_s,Nlb_trial_s,Nlb_test_s,basis_type_trial_s,0,0,basis_type_test_s,0,0,0);
        A=A_D1+A_D2+(A_C1+A_C2)+A_R;        
        %% �߽������Ĵ���
         [boundary_FE_nodes_s,boundary_edges_s]=generate_boundary_nodes_edges_s_t(N_x_FE_nodes_s,N_y_FE_nodes_s,n_x,n_y);
        %% ��ʼ�������趨
        if n==1
            X_old(1:size(Pb_trial_s,2),1)=0;%��һ�α߽�����
        elseif n>1
            clear X_old;
            X_old=initial_boundary_rectangular_t(P_old,T_old,Tb_trial_s_old,X_old0,Nlb_trial_s,Pb_trial_s,basis_type_trial_s);
        end   
    end
    
    %% ��������� ������⣬���ƽ�
    t_m=m*t_deta;
    t_m1=(m+1)*t_deta;
    %% ǰһ��ʱ��
    A_tm=A;
    b_tm=assemble_vector_2D_s_t('function_f_t',Gauss_nodes_number_2D,P,T,Pb_test_s,Tb_test_s,Nlb_test_s,basis_type_test_s,0,0,t_m);
    %% ��ǰʱ��
    A_tm1=A;
    b_tm1=assemble_vector_2D_s_t('function_f_t',Gauss_nodes_number_2D,P,T,Pb_test_s,Tb_test_s,Nlb_test_s,basis_type_test_s,0,0,t_m1);
    %% �߽�������ʩ��
    A_tilde=M/t_deta+theta*A_tm1;
    b_tilde=theta*b_tm1+(1-theta)*b_tm+(M/t_deta-(1-theta)*A_tm)*X_old;
    %% �߽�������ʩ��
    [A_tilde,b_tilde]= treat_boundary_Dirichlet_s_t('function_Dirichlet_s_t',A_tilde,b_tilde,Pb_test_s,boundary_FE_nodes_s,t_m1);
    %% ɱ����Ԫ���ƶ��߽�
     uc=0.92*u_b;%����һά
%      if n>1
         [A_tilde,b_tilde]=killed_element_method(A_tilde,b_tilde,Tb_test_s,P,T,X_old,Tb_trial_s,basis_type_trial_s,Nlb_trial_s,uc);
%      end     
    %% ���Է���������
    X_old=abs(A_tilde\b_tilde);  
    %% �����������
    Z_s=reshape(X_old(1:size(Pb_test_s,2)),N_x_FE_nodes_s,N_y_FE_nodes_s);
    n_jet=7;myjet = jet(64*n_jet);myjet(64*n_jet,:) = [1 1 1];%���䲢�������ֵΪ��ɫ
    figure(n);colormap(myjet);
    contourf(YY_s,XX_s,Z_s);colorbar;
    axis([left_i,right_i,bottom_i,top_i]);
    xlabel('dm');ylabel('dm');
    set(gcf,'color','w')%�趨����
    set(gcf,'position',[500,200,580,500])
    set(gca,'position',[0.09,0.1,0.77,0.82]);
    
    annotation(figure(n),'textbox',[0.355 0.859 0.2 0.15], 'String','t=','LineWidth',1,'LineStyle','none','FitBoxToText','off','FontSize',20);
    gg=annotation(figure(n),'textbox',[0.355 0.859 0.2 0.15], 'String',strcat('t=',num2str((m+1)*n_T),'d'),'LineWidth',1,'LineStyle','none','FitBoxToText','off','FontSize',20);
    h=colorbar;set(get(h,'title'),'string','C(-)');
 
 
%% ĳһ�����ϵ�Ũ��
% figure(n);
% plot(XX_s,Z_s(20,:),'color',mycolor(m+1,:),'LineWidth',2); hold on
% set(gcf,'position',[360,198,560,420]);set(gca,'position',[0.1,0.1,0.86,0.82]);
% xlabel('dm');ylabel('c/%');


%% ����GIF   
    im = frame2im(getframe(gcf));    
    [I, map] = rgb2ind(im,256);    
    address=strcat('D:\������\������\�ڰ���\����\GIF\',num2str(n),'.gif');
    if (m==0)
        imwrite(I,map,address,'gif','Loopcount',inf,'DelayTime',0.2);
    else
        imwrite(I,map,address,'gif','WriteMode','append','DelayTime',0.2);
    end
    
    %% �ﵽ�趨ֵ�������趨�����ʱ�����
%     %��һ�����񻮷��൱�������θ�ʴ�У�����С�����ȴ��������Ũ�ȴﵽ�ٽ�ֵʱ����һ�����񻮷֣����ǵ�һ�������ж�����
%     if n==1
%             X_d1=0.05;%�ٽ������
%            X_c=left_i+n*X_d1;Y_c=(bottom_i+top_i)/2;%���䴦����
% %             Y_c=bottom_i+n*X_d1;X_c=(left_i+right_i)/2;%���䴦����
%             uu=critical_condition(P,T,X_c,Y_c,X_old,Tb_trial_s,basis_type_trial_s);%�����ٽ���Ũ��ֵ
%             cc=0.85*u_b;%�ٽ����Ũ��
%             [uu,cc]
%             time=[];
%         if  uu>cc%�ٽ����Ũ��
%              time(n)=m;
%             [left,right,bottom,top,m,P_old,T_old,Tb_trial_s_old,X_old0,dx,dy,n_x,n_y,t_deta,Mm,n]=space_time_parameter(left_i,right_i,bottom_i,top_i,n,X_d1,m,P,T,Tb_trial_s,X_old,t_end,t_initial,n_element);
%       
%         end
%     % ֮�����񻮷ֵı�׼ʱ���趨���Ũ�ȴﵽ�趨ֵʱ��������ɢ��һ������ˣ�Ϊ�˴ﵽ���õľ��ȣ���Ҫ���»�������  
%     elseif n>1
%             X_d=0.05;%���»�������ĺ��
%              X_c=left_i+n*X_d;Y_c=(bottom_i+top_i)/2;%�ж��Ƿ��趨��
% %             Y_c=bottom_i+n*X_d;X_c=(left_i+right_i)/2;%���䴦����  
%             uu=critical_condition(P,T,X_c,Y_c,X_old,Tb_trial_s,basis_type_trial_s);%�����ٽ���Ũ��ֵ
%             uc1=0.85*u_b;
%              [uu,uc1];
%         if  uu>uc1%���趨���»��������ȴ���Ũ�ȴ������ֵʱ���������񻮷�
%             time(n)=m;
%             [left,right,bottom,top,m,P_old,T_old,Tb_trial_s_old,X_old0,dx,dy,n_x,n_y,t_deta,Mm,n]=space_time_parameter(left_i,right_i,bottom_i,top_i,n,X_d,m,P,T,Tb_trial_s,X_old,t_end,t_initial,n_element);
%         figure(10)
%         plot(1:length(time),time);
%         end
%      end
%  [m,Mm-1,n,size(Pb_test_s,2)]
     m=m+1;
   
 [m,Mm-1,n] 
 delete(gg);
end
gg=annotation(figure(n),'textbox',[0.355 0.859 0.2 0.15], 'String',strcat('t=',num2str((m+1)*n_T),'d'),'LineWidth',1,'LineStyle','none','FitBoxToText','off','FontSize',20);
rmpath(genpath('D:\������\������\�ڰ���\����\AFEM_2D_rectangular\basic_subroutine'));
rmpath(genpath('D:\������\������\�ڰ���\����\AFEM_2D_rectangular\other_subroutine'));
rmpath(genpath('D:\������\������\�ڰ���\����\AFEM_2D_rectangular\treat_boudary_subroutine'));

   