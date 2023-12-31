%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                       %%
%%                    Gen_Randperm                       %%
%%                                                       %%
%%                     Main Program                      %%
%%                Version 1.0 ; Aug 2022                 %%
%%                                                       %%
%%                  Author:  Yudi Wang                   %%
%%                Supervisor: Libing Du                  %%
%%                                                       %%
%%       Realized at Southwest Petroleum University      %%
%%                        China                          %%
%%                      Year 2022                        %%
%%                                                       %%
%%            Please read enclosed .txt file             %%
%%                                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%生成随机数矩阵2^1-2^27之间,开机就进行
%每次开机就进行一次运算，保存为一个元胞数组

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% A. SMALL DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
k=569847420;                                                               %输入数据网格数量
if k<=2^27
rand_mat_1=cell(27,1);
for i=1:27
    c=single(2^i);
    rand_mat_1{i}=single(randperm(c));
end
save E:\code\code\rand_mat_small;

%%
%生成随机数矩阵大于2^27后,进行for循环根据数据大小
%以下为10的案例，每次循环成倍增长
% %example：
% %A=randperm(10);
% % i=2;
% % for j=1:i
% % M=rand(1,2^(j-1)*10)>0.5;
% % B=A(A)+2^(j-1)*10;
% % AA=A+M*2^(j-1)*10;
% % BB=B-M(A)*2^(j-1)*10;
% % A=[AA,BB];
% % end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% B. BIG DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif k>2^27
    c=int32(2^27);                                                         %定义数据类型
    rand_mat_1=int32(randperm(c));                                         %生成2^27的乱序数组
    i=ceil(log2(k/(2^27)));                                                %判断循环次数
    for j=1:i
    M=int32(rand(1,(2^(j-1))*(2^27)));                                     %建立（0，1）判断矩阵
    rand_mat_2=rand_mat_1(rand_mat_1)+(2^(j-1))*(2^27);                    %生成下段乱序矩阵
    rand_mat_up=rand_mat_1+M*((2^(j-1))*(2^27));                           %将上段乱序矩阵根据判断矩阵提取部分放下段
    rand_mat_down=rand_mat_2-M(rand_mat_1)*(2^(j-1))*(2^27);               %将下段乱序矩阵根据判断矩阵提取部分放上段
    rand_mat_1=[rand_mat_up,rand_mat_down];                                %更新后上下段矩阵叠合成一个矩阵
    end
end

 save E:\code\code\rand_mat_lager;