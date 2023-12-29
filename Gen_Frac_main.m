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
%����Լ������������ΪȨֵ����
%xΪ�ߺţ�yΪ���ţ�z��1��ʼ
%Ҫ������Ӧ��ʼ������(zֵΪ�������ݵ�����������Ϊ�ߵ���������ˣ�����֤)

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
seismic=read_segy_file('G:\HangZhou\���ն�����\int_rel_����_3-8_��һ��.sgy');

header=seismic.headers;
xline_step=header(3,2)-header(3,1);                                        %����ż�Ĳ���

x_values=s_gh(seismic,'ffid');                                             %�ߺŵ�ͷ��ȡ
y_values=s_gh(seismic,'cdp');                                              %���ŵ�ͷ��ȡ

y=max(y_values)/xline_step-min(y_values)/xline_step+1;                     %�����ж�������

NUM=size(y_values);
x=NUM(1,2)/y;                                                              %�����ж�������

value=seismic.traces;                                                      %valueΪÿ�����Ӧ��Լ��ֵ
 
[m,n]=size(value);                                                         %m�Ǵ�����������n���ߵ��ŵĳ˻�
K=m*n;                                                                     %��Լ��������

z=m;                                                                       %z�����ж��ٸ�����

load 'E:\code\code\rand_mat_small';                                            %���뿪�����ɵ����������
load 'E:\code\code\rand_mat_lager';                                            %���뿪�����ɵ����������

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% B. RANDPERM DATA %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%������һ������ͬ�������ݸ�ʽ�����ݴ�С���ܲ�ͬ
% value=(value-min(min(value)))/(max(max(value))-min(min(value)));


Log_Num=ceil(log2(K));
if Log_Num<=27
   index_1=rand_mat_1{Log_Num};                                            %��ȡ���Һ�index
else
   index_1=rand_mat_2;
end

index_1=index_1(index_1<=K);                                               %��ȡ��������������������

value_1=value(index_1);                                                    %�����Լ��ֵ

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% C. FIND DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=1000;                                                                    %�����ѷ�����
frac_corner=zeros(N,1);                                                    %�����վ�����ǰ�����ڴ�
 
size_value_1=size(value_1);                                                %����ͬ���������������
judge_rand=rand(size_value_1);

threshold=0.2;                                                             %�жϴ�С���趨��ֵ
tic

judge_difference=value_1-(judge_rand+threshold);

l=find(judge_difference>=0,N);                                             %Ѱ�Һ�����

frac_corner(:,1)=index_1(l);
frac_corner(:,2)=value_1(l);

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% D. GENERATED COORDINATES %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [frac_x,frac_y,frac_z]=ind2sub([x,y,z],frac_corner(:,1));

frac_x=ceil((ceil(frac_corner(:,1)/z))/y);                                 %��indexת��Ϊ�ߺ�--x����
frac_y=mod((ceil(frac_corner(:,1)/z)),y);                                  %��indexת��Ϊ����--y����
if frac_y==0
    frac_y=header(3,1)+y-1;
end
frac_z=mod(frac_corner(:,1),z);                                            %��indexת��Ϊ����--z����
if frac_z==0
    frac_z=z;
end

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% E. GENERATED FRACTURES %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

point=[frac_x,frac_y,frac_z];                                              %���ĵ�λ��

%�޸����ĵ�λ�ã�ȫ�������λ����б��֮���������
%  x1=20630725;y1=4193550;z1=2440;
 x1=704268.4;y1=4471381.6;z1=5500; 
 
% frac_x=x1+25*frac_x;
% frac_y=y1+12.5*frac_y;
% frac_z=z1+2*frac_z; 
frac_x=x1+25*frac_x;
frac_y=y1+12.5*frac_y;
frac_z=z1+10*frac_z; 
%�޸����ĵ�λ�ã�ȫ�������λ����б��֮���������
%  x1=-2;y1=-2;z1=-15.4322;
% 
%  
% frac_x=x1+0.1*frac_x;
% frac_y=y1+0.1*frac_y;
% frac_z=z1+0.1*frac_z; 

point=[frac_x,frac_y,frac_z];                                              %���ĵ�λ��

side_num=4;                                                                %�������ɼ�����
axial_ratio=1;                                                             %�������ɳ������

Ltype=2;                                                                   %�ѷ�Ƭ�볤�ֲ����ͣ����ȷֲ������ޣ����ޣ���̬�ֲ�����ֵ�����
Lpara1=20;                                                                  %�ѷ�Ƭ�볤�ֲ�����1
Lpara2=40;                                                                  %�ѷ�Ƭ�볤�ֲ�����2

Ttype=1;                                                                   %�ѷ�Ƭ����ֲ�����
Tpara1=pi/3;                                                               %�ѷ�Ƭ����ֲ�����1
Tpara2=pi/3;                                                              %�ѷ�Ƭ����ֲ�����2

Itype=1;                                                                   %�ѷ�Ƭ��Ƿֲ�����: ��̬�ֲ�����ֵ�����
Ipara1=pi/2;                                                               %�ѷ�Ƭ��Ƿֲ�����1
Ipara2=pi/10;                                                              %�ѷ�Ƭ��Ƿֲ�����2

Radius= Random_Sampling(Ltype,Lpara1,Lpara2,N,1);                          %�ѷ�Ƭ�İ볤����ǡ������������
Trend = Random_Sampling(Ttype,Tpara1,Tpara2,N,1);
Inclination=Random_Sampling(Itype,Ipara1,Ipara2,N,1);

[new_plys_plot,X_polt,Y_polt,Z_polt]=Disk_gen_111(Radius,frac_x,frac_y,frac_z,Trend,Inclination,side_num,axial_ratio,N);

toc

% for i=1:N
%     [Disk_data,R,Fit]=Disk_gen(Radius(i),frac_x(i),frac_y(i),frac_z(i),Trend(i),Inclination(i),side_num);
% end

% plys=cell(N,1);                                                          %�����ѷ��Ԫ������
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
plot3(new_plys_plot(:,1),new_plys_plot(:,2),new_plys_plot(:,3));%һ���Ի�ͼ
set(gca,'ZDir','reverse');
hold on
fill3(X_polt,Y_polt,Z_polt,'r');
set(gca,'ZDir','reverse');
hold on

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% G. ���³�ͼ�õĴ��� %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%������ά��

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
%��Χ��ͼ��1.�жϷ�Χ��С��ͼ��2.ָ����Χ��ͼ��
%һ���Ի�ͼ���жϷ�Χ��С��ͼ
% new_plys_plot=cell2mat(new_plys);
if x>200||y>200&&N>10
%     polyin = alphaShape(new_plys_plot);
%     plot(polyin)
plot3(new_plys_plot(:,1),new_plys_plot(:,2),new_plys_plot(:,3));%һ���Ի�ͼ
set(gca,'ZDir','reverse');
f1=figure;

else
    DrawPolys3D(plys);
end
%ԭʼ��ͼ��ָ����ͼ��Χ
local_x_min=50;local_x_max=100;
local_y_min=20;local_y_max=40;
local_z_min=10;loacal_z_max=20;
%��Щλ�ö����Խ����޸�
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


 

[mm,nn]=find(isnan(new_plys_plot));    % �ҳ�NaN����λ��
new_plys_plot(mm,:)=[]; %ɾ������NaN����

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