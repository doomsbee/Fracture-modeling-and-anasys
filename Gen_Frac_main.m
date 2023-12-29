%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                       %%
%%                      Gen_Frac                         %%
%%                                                       %%
%%                     Main Program                      %%
%%                Version 1.0 ; Aug 2022                 %%
%%                                                       %%
%%                  Author:  Yudi Wang                   %%
%%                Supervisor: Libing Du                  %%
%%                                                       %%
%%       Realized at Southwest Petroleum University      %%
%%                        China                          %%
%%                     Year 2022                         %%
%%                                                       %%
%%            Please read enclosed .txt file             %%
%%                                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%输入约束地震数据作为权值矩阵
%x为线号，y为道号，z从1开始
%要根据相应起始深度添加(z值为输入数据的行数，列数为线道号数量相乘，可验证)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% A. INPUT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all

% load 'G:\code\code\dist_bound';  
% 
% [x,y,z]=size(dist_bound);
% 
% K=x*y*z;
% 
% value=dist_bound;
seismic=read_segy_file('G:\HangZhou\果勒东数据\int_rel_连井_3-8_归一化.sgy');

header=seismic.headers;
xline_step=header(3,2)-header(3,1);                                        %求道号间的步长

x_values=s_gh(seismic,'ffid');                                             %线号道头读取
y_values=s_gh(seismic,'cdp');                                              %道号道头读取

y=max(y_values)/xline_step-min(y_values)/xline_step+1;                     %计算有多少条道

NUM=size(y_values);
x=NUM(1,2)/y;                                                              %计算有多少条线

value=seismic.traces;                                                      %value为每个点对应的约束值
 
[m,n]=size(value);                                                         %m是垂向数据量，n是线道号的乘积
K=m*n;                                                                     %总约束数据量

z=m;                                                                       %z垂向有多少个数据

load 'E:\code\code\rand_mat_small';                                            %导入开机生成的随机排列数
load 'E:\code\code\rand_mat_lager';                                            %导入开机生成的随机排列数

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% B. RANDPERM DATA %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%做个归一化，不同地震数据格式，数据大小可能不同
% value=(value-min(min(value)))/(max(max(value))-min(min(value)));


Log_Num=ceil(log2(K));
if Log_Num<=27
   index_1=rand_mat_1{Log_Num};                                            %提取打乱后index
else
   index_1=rand_mat_2;
end

index_1=index_1(index_1<=K);                                               %提取等于数据量的乱序数据

value_1=value(index_1);                                                    %乱序后约束值

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% C. FIND DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=1000;                                                                    %定义裂缝条数
frac_corner=zeros(N,1);                                                    %建立空矩阵，提前分配内存
 
size_value_1=size(value_1);                                                %生成同样数量随机数矩阵
judge_rand=rand(size_value_1);

threshold=0.2;                                                             %判断大小并设定阈值
tic

judge_difference=value_1-(judge_rand+threshold);

l=find(judge_difference>=0,N);                                             %寻找合适数

frac_corner(:,1)=index_1(l);
frac_corner(:,2)=value_1(l);

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% D. GENERATED COORDINATES %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [frac_x,frac_y,frac_z]=ind2sub([x,y,z],frac_corner(:,1));

frac_x=ceil((ceil(frac_corner(:,1)/z))/y);                                 %将index转化为线号--x坐标
frac_y=mod((ceil(frac_corner(:,1)/z)),y);                                  %将index转化为道号--y坐标
if frac_y==0
    frac_y=header(3,1)+y-1;
end
frac_z=mod(frac_corner(:,1),z);                                            %将index转化为垂向--z坐标
if frac_z==0
    frac_z=z;
end

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% E. GENERATED FRACTURES %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

point=[frac_x,frac_y,frac_z];                                              %中心点位置

%修改中心点位置，全部算完改位置倾斜角之类的有问题
%  x1=20630725;y1=4193550;z1=2440;
 x1=704268.4;y1=4471381.6;z1=5500; 
 
% frac_x=x1+25*frac_x;
% frac_y=y1+12.5*frac_y;
% frac_z=z1+2*frac_z; 
frac_x=x1+25*frac_x;
frac_y=y1+12.5*frac_y;
frac_z=z1+10*frac_z; 
%修改中心点位置，全部算完改位置倾斜角之类的有问题
%  x1=-2;y1=-2;z1=-15.4322;
% 
%  
% frac_x=x1+0.1*frac_x;
% frac_y=y1+0.1*frac_y;
% frac_z=z1+0.1*frac_z; 

point=[frac_x,frac_y,frac_z];                                              %中心点位置

side_num=4;                                                                %定义生成几边形
axial_ratio=1;                                                             %定义生成长短轴比

Ltype=2;                                                                   %裂缝片半长分布类型：均匀分布（上限，下限）正态分布（均值，方差）
Lpara1=20;                                                                  %裂缝片半长分布参数1
Lpara2=40;                                                                  %裂缝片半长分布参数2

Ttype=1;                                                                   %裂缝片走向分布类型
Tpara1=pi/3;                                                               %裂缝片走向分布参数1
Tpara2=pi/3;                                                              %裂缝片走向分布参数2

Itype=1;                                                                   %裂缝片倾角分布类型: 正态分布（均值，方差）
Ipara1=pi/2;                                                               %裂缝片倾角分布参数1
Ipara2=pi/10;                                                              %裂缝片倾角分布参数2

Radius= Random_Sampling(Ltype,Lpara1,Lpara2,N,1);                          %裂缝片的半长、倾角、走向随机生成
Trend = Random_Sampling(Ttype,Tpara1,Tpara2,N,1);
Inclination=Random_Sampling(Itype,Ipara1,Ipara2,N,1);

[new_plys_plot,X_polt,Y_polt,Z_polt]=Disk_gen_111(Radius,frac_x,frac_y,frac_z,Trend,Inclination,side_num,axial_ratio,N);

toc

% for i=1:N
%     [Disk_data,R,Fit]=Disk_gen(Radius(i),frac_x(i),frac_y(i),frac_z(i),Trend(i),Inclination(i),side_num);
% end

% plys=cell(N,1);                                                          %建立裂缝空元胞矩阵
% new_plys=cell(N,1);
% for i=1:N
%     [Disk_data,R,Fit]=Disk_gen_New_1(Radius(i),frac_x(i),frac_y(i),frac_z(i),Trend(i),Inclination(i),side_num,axial_ratio);
%     plys{i}=Fit;
% %     Fit(side_num+1,:)=[Fit(1,1),Fit(1,2),Fit(1,3)];
% %     Fit(side_num+2,:)=[NaN,NaN,NaN];
%     new_plys{i}=Fit;
% end

% save('new_plys_plot.mat','new_plys_plot');
% save('X_polt.mat','X_polt');save('Y_polt.mat','Y_polt');save('Z_polt.mat','Z_polt');

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% F. PLOT FRACTURES %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot3(new_plys_plot(:,1),new_plys_plot(:,2),new_plys_plot(:,3));%一次性绘图
set(gca,'ZDir','reverse');
hold on
fill3(X_polt,Y_polt,Z_polt,'r');
set(gca,'ZDir','reverse');
hold on

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% G. 文章出图用的代码 %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%矩阵三维化

for k=1:z
    for j=1:x
        for i=1:y
            curv(i,j,k)=value(k,y*(j-1)+i);
       end
    end
end 
%%

[X,Y] = meshgrid(0:.05:10);
Z = Y.*sin(X) - X.*cos(Y);
% Z=X+Y-1;
s = surf(Y,X,Z,'EdgeColor','interp','FaceColor','interp');
% xlabel('X'),ylabel('Y'),zlabel('Z');
set(gca,'ZDir','reverse');
set(gcf,'unit','centimeters','position',[10 5 6 4]);
set(gca,'FontName','Times New Roman','FontSize',10.5,'LineWidth',0.5);
%范围绘图（1.判断范围大小绘图；2.指定范围绘图）
%一次性绘图，判断范围大小绘图
% new_plys_plot=cell2mat(new_plys);
if x>200||y>200&&N>10
%     polyin = alphaShape(new_plys_plot);
%     plot(polyin)
plot3(new_plys_plot(:,1),new_plys_plot(:,2),new_plys_plot(:,3));%一次性绘图
set(gca,'ZDir','reverse');
f1=figure;

else
    DrawPolys3D(plys);
end
%原始绘图，指定绘图范围
local_x_min=50;local_x_max=100;
local_y_min=20;local_y_max=40;
local_z_min=10;loacal_z_max=20;
%这些位置都可以进行修改
K=1;
for i=1:N
    if point(i,1)>=local_x_min&&point(i,1)<=local_x_max&&point(i,2)>=local_y_min&&...
            point(i,2)<=local_y_max&&point(i,3)>=local_z_min&&point(i,3)<=loacal_z_max
        local_XYZ=plys{i};
        fill3(local_XYZ(:,1),local_XYZ(:,2),local_XYZ(:,3),'red');
        hold on
        local_plys{K}=plys{i};
        K=K+1;
    else
        continue
    end
end

%  DrawPolys3D(local_plys);


 

[mm,nn]=find(isnan(new_plys_plot));    % 找出NaN数据位置
new_plys_plot(mm,:)=[]; %删除含有NaN的行

%  x1=20630725;y1=4193550;z1=2440;
%  x2=20646950;y2=4202737;z2=3800;
%  
% new_plys_plot(:,1)=x1+25*new_plys_plot(:,1);
% new_plys_plot(:,2)=y1+12.5*new_plys_plot(:,2);
% new_plys_plot(:,3)=z1+2*new_plys_plot(:,3); 

add_num=length(new_plys_plot)/N;
for j=1:N
    for jj=1:add_num
       new_plys_plot((add_num)*(j-1)+jj,4)=j;
       new_plys_plot((add_num)*(j-1)+jj,5)=j;
    end
end

xlswrite('H:\code\code\new_plys_plot.xls',new_plys_plot);
save H:\code\code\new_plys_plot.txt -ascii new_plys_plot