clear
clc
tic
clf
format shortE %修改数字类型
addpath(genpath('D:\大论文\大论文\第八章\程序\AFEM_2D_rectangular\basic_subroutine'));%打开文件夹中的相关程序
addpath(genpath('D:\大论文\大论文\第八章\程序\AFEM_2D_rectangular\other_subroutine'));
addpath(genpath('D:\大论文\大论文\第八章\程序\AFEM_2D_rectangular\treat_boudary_subroutine'));
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%需要输入的参数%%%%%%%%%%%%%%%%%%%%%%%%%
%% 网格参数(左右上下边界)
left_i=0;right_i=1;bottom_i=0;top_i=1;%区域
left=left_i;right=right_i;bottom=bottom_i;top=top_i;%区域
%% 时间迭代参数
theta=1/2;t_initial=0;t_end=2;TT=2*3650;
%% 基函数类型，线性和二次
basis_type_trial_s=241;basis_type_test_s=241;
%% 高斯积分参数
Gauss_nodes_number_2D=9;%用于计算二重积分
Gauss_nodes_number_1D=4;%用于计算线积分
%% 离散参数
n_element=6;
uc=0.8;
dx=1/2^(n_element);dy=1/2^(n_element);
n_x=(right-left)/dx;%x方向网格单元数
n_y=(top-bottom)/dy;%x方向网格单元数
t_deta=dx/2;%时间迭代步长,
Mm=round((t_end-t_initial)/t_deta);%时间迭代步数
n_T=round(TT/Mm);%时间间隔步数
u_b=1;0.308;%边界浓度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
n=1;%代表划分单元格的次数
m=0;%时间迭代次数
% %% colormap的创建
% % nr=Mm;%修改输出饱和度间隔
% nrgb=Mm; %输出线条间隔数    
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

while m<Mm-1% 时间迭代
      
     if m==0
        %% 网格坐标及网格单元局部编号和整体编号映射关系（与选择的元次无关,因此选择241）
        [P,T]=Generate_P_T_rectangular(241,left,right,bottom,top,n_x,n_y);  %网格结点及映射矩阵
        %% 根据基函数类型定义网格及有限元结点坐标及单元局部编号和整体编号映射关系（与元次有关）
        [Nlb_trial_s,Pb_trial_s,Tb_trial_s,N_x_FE_nodes_s,N_y_FE_nodes_s,]=Generate_Pb_Tb_rectangular(basis_type_trial_s,left,right,bottom,top,n_x,n_y);
        [Nlb_test_s,Pb_test_s,Tb_test_s,N_x_FE_nodes_s,N_y_FE_nodes_s]=Generate_Pb_Tb_rectangular(basis_type_test_s,left,right,bottom,top,n_x,n_y);
        %% 有限元信息
        dx_elment_nodes_s=(right-left)/(N_x_FE_nodes_s-1);%向量场的第一个分量在x方向相邻有限元节点之间的距离，N_x_FE_nodes_u为x方向有限元节点数
        dy_elment_nodes_s=(top-bottom)/(N_y_FE_nodes_s-1);%向量场的第二个分量在y方向相邻有限元节点之间的距离
        XX_s=left:dx_elment_nodes_s:right;%x方向有限元刻度
        YY_s=bottom:dy_elment_nodes_s:top;
        %% 质量阵
          M  = assemble_matrix_2D_s('function_mass_coef',                  Gauss_nodes_number_2D,P,T,Pb_trial_s,Tb_trial_s,Pb_test_s,Tb_test_s,Nlb_trial_s,Nlb_test_s,basis_type_trial_s,0,0,basis_type_test_s,0,0);        
        A_D1 = assemble_matrix_2D_s_t('function_diffusion_t',              Gauss_nodes_number_2D,P,T,Pb_trial_s,Tb_trial_s,Pb_test_s,Tb_test_s,Nlb_trial_s,Nlb_test_s,basis_type_trial_s,1,0,basis_type_test_s,1,0,0);
        A_D2 = assemble_matrix_2D_s_t('function_diffusion_t',              Gauss_nodes_number_2D,P,T,Pb_trial_s,Tb_trial_s,Pb_test_s,Tb_test_s,Nlb_trial_s,Nlb_test_s,basis_type_trial_s,0,1,basis_type_test_s,0,1,0);
        A_C1 = assemble_matrix_2D_s_t('function_convection_x_t',           Gauss_nodes_number_2D,P,T,Pb_trial_s,Tb_trial_s,Pb_test_s,Tb_test_s,Nlb_trial_s,Nlb_test_s,basis_type_trial_s,1,0,basis_type_test_s,0,0,0);
        A_C2 = assemble_matrix_2D_s_t('function_convection_y_t',           Gauss_nodes_number_2D,P,T,Pb_trial_s,Tb_trial_s,Pb_test_s,Tb_test_s,Nlb_trial_s,Nlb_test_s,basis_type_trial_s,0,1,basis_type_test_s,0,0,0);
        A_R  = assemble_matrix_2D_s_t('function_reaction_t',               Gauss_nodes_number_2D,P,T,Pb_trial_s,Tb_trial_s,Pb_test_s,Tb_test_s,Nlb_trial_s,Nlb_test_s,basis_type_trial_s,0,0,basis_type_test_s,0,0,0);
        A=A_D1+A_D2+(A_C1+A_C2)+A_R;        
        %% 边界条件的处理
         [boundary_FE_nodes_s,boundary_edges_s]=generate_boundary_nodes_edges_s_t(N_x_FE_nodes_s,N_y_FE_nodes_s,n_x,n_y);
        %% 初始条件的设定
        if n==1
            X_old(1:size(Pb_trial_s,2),1)=0;%第一次边界条件
        elseif n>1
            clear X_old;
            X_old=initial_boundary_rectangular_t(P_old,T_old,Tb_trial_s_old,X_old0,Nlb_trial_s,Pb_trial_s,basis_type_trial_s);
        end   
    end
    
    %% 调用求解器 返回真解，近似解
    t_m=m*t_deta;
    t_m1=(m+1)*t_deta;
    %% 前一段时间
    A_tm=A;
    b_tm=assemble_vector_2D_s_t('function_f_t',Gauss_nodes_number_2D,P,T,Pb_test_s,Tb_test_s,Nlb_test_s,basis_type_test_s,0,0,t_m);
    %% 当前时间
    A_tm1=A;
    b_tm1=assemble_vector_2D_s_t('function_f_t',Gauss_nodes_number_2D,P,T,Pb_test_s,Tb_test_s,Nlb_test_s,basis_type_test_s,0,0,t_m1);
    %% 边界条件的施加
    A_tilde=M/t_deta+theta*A_tm1;
    b_tilde=theta*b_tm1+(1-theta)*b_tm+(M/t_deta-(1-theta)*A_tm)*X_old;
    %% 边界条件的施加
    [A_tilde,b_tilde]= treat_boundary_Dirichlet_s_t('function_Dirichlet_s_t',A_tilde,b_tilde,Pb_test_s,boundary_FE_nodes_s,t_m1);
    %% 杀死单元法移动边界
     uc=0.92*u_b;%对于一维
%      if n>1
         [A_tilde,b_tilde]=killed_element_method(A_tilde,b_tilde,Tb_test_s,P,T,X_old,Tb_trial_s,basis_type_trial_s,Nlb_trial_s,uc);
%      end     
    %% 线性方程组的求解
    X_old=abs(A_tilde\b_tilde);  
    %% 计算结果的输出
    Z_s=reshape(X_old(1:size(Pb_test_s,2)),N_x_FE_nodes_s,N_y_FE_nodes_s);
    n_jet=7;myjet = jet(64*n_jet);myjet(64*n_jet,:) = [1 1 1];%扩充并令其最大值为白色
    figure(n);colormap(myjet);
    contourf(YY_s,XX_s,Z_s);colorbar;
    axis([left_i,right_i,bottom_i,top_i]);
    xlabel('dm');ylabel('dm');
    set(gcf,'color','w')%设定背景
    set(gcf,'position',[500,200,580,500])
    set(gca,'position',[0.09,0.1,0.77,0.82]);
    
    annotation(figure(n),'textbox',[0.355 0.859 0.2 0.15], 'String','t=','LineWidth',1,'LineStyle','none','FitBoxToText','off','FontSize',20);
    gg=annotation(figure(n),'textbox',[0.355 0.859 0.2 0.15], 'String',strcat('t=',num2str((m+1)*n_T),'d'),'LineWidth',1,'LineStyle','none','FitBoxToText','off','FontSize',20);
    h=colorbar;set(get(h,'title'),'string','C(-)');
 
 
%% 某一条线上的浓度
% figure(n);
% plot(XX_s,Z_s(20,:),'color',mycolor(m+1,:),'LineWidth',2); hold on
% set(gcf,'position',[360,198,560,420]);set(gca,'position',[0.1,0.1,0.86,0.82]);
% xlabel('dm');ylabel('c/%');


%% 生成GIF   
    im = frame2im(getframe(gcf));    
    [I, map] = rgb2ind(im,256);    
    address=strcat('D:\大论文\大论文\第八章\程序\GIF\',num2str(n),'.gif');
    if (m==0)
        imwrite(I,map,address,'gif','Loopcount',inf,'DelayTime',0.2);
    else
        imwrite(I,map,address,'gif','WriteMode','append','DelayTime',0.2);
    end
    
    %% 达到设定值后，重新设定网格和时间参数
%     %第一次网格划分相当于硫酸盐腐蚀中，当最小剥落厚度处的硫酸根浓度达到临界值时进行一次网格划分，这是第一层剥落的判断条件
%     if n==1
%             X_d1=0.05;%临界剥落厚度
%            X_c=left_i+n*X_d1;Y_c=(bottom_i+top_i)/2;%剥落处坐标
% %             Y_c=bottom_i+n*X_d1;X_c=(left_i+right_i)/2;%剥落处坐标
%             uu=critical_condition(P,T,X_c,Y_c,X_old,Tb_trial_s,basis_type_trial_s);%计算临界点的浓度值
%             cc=0.85*u_b;%临界剥落浓度
%             [uu,cc]
%             time=[];
%         if  uu>cc%临界剥落浓度
%              time(n)=m;
%             [left,right,bottom,top,m,P_old,T_old,Tb_trial_s_old,X_old0,dx,dy,n_x,n_y,t_deta,Mm,n]=space_time_parameter(left_i,right_i,bottom_i,top_i,n,X_d1,m,P,T,Tb_trial_s,X_old,t_end,t_initial,n_element);
%       
%         end
%     % 之后网格划分的标准时当设定点的浓度达到设定值时，表明扩散有一定深度了，为了达到更好的精度，需要重新划分网格  
%     elseif n>1
%             X_d=0.05;%重新划分网格的厚度
%              X_c=left_i+n*X_d;Y_c=(bottom_i+top_i)/2;%判断是否设定的
% %             Y_c=bottom_i+n*X_d;X_c=(left_i+right_i)/2;%剥落处坐标  
%             uu=critical_condition(P,T,X_c,Y_c,X_old,Tb_trial_s,basis_type_trial_s);%计算临界点的浓度值
%             uc1=0.85*u_b;
%              [uu,uc1];
%         if  uu>uc1%当设定重新划分网格厚度处的浓度大于最大值时，进行网格划分
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
rmpath(genpath('D:\大论文\大论文\第八章\程序\AFEM_2D_rectangular\basic_subroutine'));
rmpath(genpath('D:\大论文\大论文\第八章\程序\AFEM_2D_rectangular\other_subroutine'));
rmpath(genpath('D:\大论文\大论文\第八章\程序\AFEM_2D_rectangular\treat_boudary_subroutine'));

   