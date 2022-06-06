function [A_tilde,b_tilde]=killed_element_method(A_tilde,b_tilde,Tb_test_s,P,T,X_old,Tb_trial_s,basis_type_trial_s,Nlb_trial_s,cu)
%% 杀死单元法
    for N_element=1:size(Tb_test_s,2)% 
        
        vertices=P(:,T(:,N_element));%提取单元的网格结点坐标
        X_o=1/2*(vertices(1,2)+vertices(1,1));Y_o=1/2*(vertices(2,4)+vertices(2,1));%计算单元中心点坐标
        uu=critical_condition(P,T,X_o,Y_o,X_old,Tb_trial_s,basis_type_trial_s);%计算临界点的浓度值
        
        if uu>cu%判断矩形单元中心点浓度是否达到临界值
            for k=1:Nlb_trial_s%循环单元的所有有限元结点
                n_node=Tb_trial_s(k,N_element);%取N_element号单元第k个有限元结点
                A_tilde(n_node,:)=0;%消去第i个方程
                A_tilde(n_node,n_node)=1;
                b_tilde(n_node,1)=1;%换成已知值
            end
        end
    end