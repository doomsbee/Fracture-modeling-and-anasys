%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                       %%
%%                   DF Model and Anasys                 %%
%%                                                       %%
%%                     Main Program                      %%
%%                Version 1.0 ; July 2022                %%
%%                                                       %%
%%                  Author:  Junfeng Lu                  %%
%%                Supervisor: Junfeng Lu                 %%
%%                                                       %%
%%            Realized at An Hui University              %%
%%              of Science and Technology                %%
%%                     Year 2022                         %%
%%                                                       %%
%%            Please read Detailed.docx file             %%
%%                                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Disk_data为8个参量的数组，前4个为裂缝片所在面的Ax+By+Cz+D=0参数，后四个为中心点坐标与裂缝片宽度


clc
clear
close all
RGB_color
tic

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% A1.参数输入区 %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Range=32;                              %%模型尺寸,为2^N

View_type=0;                           %%视图确定形式，0为采用倾角走向，1为给定法向量
View_data=[-60,30];                     %%[倾角，走向]
View_center=[Range,Range,Range]/2;     %%视图面的中心点  

RQD_angle=145;                          %%视图平面上用于测试RQD值的测试线角度，0~180
RQD_Center=[1,1,1]/2*Range;            %%默认测线旋转点假定为视图平面的中心点。
RQD_judge=5;                           %%完整岩样的设定标准值，规范假定为10

Num=50;                               %%裂缝片数量
% % Ltype=1;                               %%裂缝片半长分布类型：均匀分布（上限，下限）正态分布（均值，方差）
% % Lpara1=5;                              %%裂缝片半长分布参数1
% % Lpara2=5;                              %%裂缝片半长分布参数2

% % Itype=1;                               %%裂缝片倾角分布类型: 正态分布（均值，方差）
% % Ipara1=pi/3;                           %%裂缝片倾角分布参数1
% % Ipara2=pi/50;                          %%裂缝片倾角分布参数2
% % 
% % Ttype=1;                               %%裂缝片走向分布类型
% % Tpara1=pi/6;                           %%裂缝片走向分布参数1
% % Tpara2=pi/50;                          %%裂缝片走向分布参数2

Radius=repmat(5,Num,1);
Trend=repmat(5*pi/6,Num,1);
Inclination=repmat(pi/2,Num,1);

RQD_Vec=[cos(RQD_angle/180*pi),sin(RQD_angle/180*pi)];   %RQD测线向量
Rect=[[-1 -1 1 1]',[1 -1 -1 1]',[0 0 0 0]']*Range/2;     %%初始视图方框
Rect_test=Rect+View_center; 

%% 根据裂缝片的随机分布形式，分布参数，生成裂缝片
%%%根据视图参数产生视图平面
Center_point=rand(Num,3)*Range;                           %%裂缝片中心点随机生成
% % Radius= Random_Sampling(Ltype,Lpara1,Lpara2,Num,1);        %%裂缝片的半长、倾角、走向随机生成
% % Trend = Random_Sampling(Ttype,Tpara1,Tpara2,Num,1);
% % Inclination=Random_Sampling(Itype,Ipara1,Ipara2,Num,1);

Color=RGB_Color_index(randi(40,1),:)/255;
figure('Name','Random disk','WindowStyle','docked')
Array=cell(1,Num);
Disk_Fit=cell(1,Num);
Fracture=zeros(Num,8);
for i=1:Num
    x=Center_point(i,1);
    y=Center_point(i,2);
    z=Center_point(i,3);
    %%%生成裂缝片
    [Disk_data,R_Matrix,Fit]=Disk_gen(Radius(i),x,y,z,Trend(i),Inclination(i));  
    Array{i}=R_Matrix;
    Fracture(i,:)=Disk_data;
    Disk_Fit{i}=Fit;
end
hold on
plotcube1([Range Range Range],[0 0 0],0.05,Color);
axis square
xlabel('X')
ylabel('Y')
zlabel('Z')
set(gca, 'fontsize',15)

%% 初始视图平面与xy平面平行，视图中心点为View_center
%%%根据视图参数产生视图平面
if View_type==0                                     %%通过输入倾角，走向，确定视图面的法向量
    View_alfa=View_data(1)*pi/180;                  %%视图面的倾角
    View_beta=View_data(2)*pi/180;                  %%视图面的走向
    Aview=-sin(View_alfa)*sin(View_beta);
    Bview=sin(View_alfa)*cos(View_beta);
    Cview=cos(View_alfa);                           %%A,B,C裂缝片的方向，D裂缝片位置
else                                                %%直接输入视图面的法向量
    Aview=View_data(1);
    Bview=View_data(2);
    Cview=View_data(3);
end
N_Vec=[Aview,Bview,Cview];                          %%裂缝面的法向量
Dview=dot(N_Vec,View_center);                       %%裂缝片的位置
N_Vec0=[0,0,1];                                     %%裂缝面初始法向量，默认与x,y平面平行
Ro_axis=cross(N_Vec0,N_Vec);                        %%轴向量
Ro_axis=Ro_axis/norm(Ro_axis);                      %%向量归一化
theta = acos(dot(N_Vec,N_Vec0))/2;
q = [cos(theta) sin(theta)*Ro_axis];
R=[2*q(1).^2-1+2*q(2)^2  2*(q(2)*q(3)+q(1)*q(4)) 2*(q(2)*q(4)-q(1)*q(3));
    2*(q(2)*q(3)-q(1)*q(4)) 2*q(1)^2-1+2*q(3)^2 2*(q(3)*q(4)+q(1)*q(2));
    2*(q(2)*q(4)+q(1)*q(3)) 2*(q(3)*q(4)-q(1)*q(2)) 2*q(1)^2-1+2*q(4)^2];
Rect2=Rect*R+View_center;                           %%旋转后的裂缝方框


%% 任意视图切面上的裂缝显示与统计
%%%计算所有圆心点到裂缝片与视图面交线的距离
figure('Name','Fracture in View plane','WindowStyle','docked')
hold on
A=[Aview,Bview,Cview];                             %%视图切面的法向量
In_num=0;
Index=[];
Line_view={};                                      %%空间中与视图面相切形成的裂缝
Line_test={};                                      %%转换为平面上相切形成的裂缝
for i=1:Num
    B=Fracture(i,1:3);                             %%裂缝片法向量
    C=cross(A,B);                                  %%C为交线的方向
    sin_temp=norm(C);                              %%C的模为裂缝片法向量B与视图面法向量A的夹角正弦
    Dplane=dot(C,Fracture(i,5:7));                 %%Dplane计算垂面的位置
    Dist=dot(A,Fracture(i,5:7))-Dview;             %%裂缝片中心点到视图面的法向距离
    Length=Dist/norm(C);                           %%裂缝片中心点到裂缝片与视图面交线的距离
    if abs(Length) <  Radius(i)                    %%裂缝片与视图面存在交线，否则不相交
        In_num=In_num+1;
        Index=cat(1,Index,i);
        Value=[Dview;Fracture(i,4);Dplane];
        Matrix=[A;B;C];
        Dmaxtrix0=Matrix;
        Dmaxtrix1=Matrix;
        Dmaxtrix2=Matrix;
        Dmaxtrix3=Matrix;
        Dmaxtrix1(:,1)=Value;
        Dmaxtrix2(:,2)=Value;
        Dmaxtrix3(:,3)=Value;
        D0=det(Dmaxtrix0);
        D1=det(Dmaxtrix1);
        D2=det(Dmaxtrix2);
        D3=det(Dmaxtrix3);
        Pcen=[D1/D0;D2/D0;D3/D0];
        Distance=sqrt(Radius(i)^2-Length^2);
        Line1=Pcen'+C*Distance/sin_temp;
        line2=Pcen'-C*Distance/sin_temp;
        Line=cat(1,Line1,line2);
        Line_R=(Line-View_center)*R';
        Line_R=Line_R+View_center;        
        Line_R(:,3)=0;
        Line_view=cat(1,Line_view,Line);
        Line_test=cat(1,Line_test,Line_R);

        plot3(Line(:,1),Line(:,2),Line(:,3),'linewidth',1.5,'color','r')               %%空间裂缝片与视图面的相交裂缝
        fill3(Disk_Fit{i}(:,1),Disk_Fit{i}(:,2),Disk_Fit{i}(:,3),'b','facealpha',0.1)  %%空间裂缝片
        plot3(Line_R(:,1),Line_R(:,2),Line_R(:,3),'linewidth',1.5,'color','b')         %%转换到xy平面下的相交裂缝
    end
end

fill3(Rect2(:,1),Rect2(:,2),Rect2(:,3),'g','facealpha',0.1)
plotcube1([Range Range Range],[0 0 0],0.02,Color);
axis square


%% 根据视图面面边框，获得边框内的裂缝
Rect_test=cat(1,Rect_test,Rect_test(1,:));       %闭合计算边框
Poly2=Rect_test(:,1:2);
Line_D=Line_test;
Line_index=[];
Point_plot=[];
for i=1:In_num
    %线段点是否都在边框内部
    %是则无需修正，否则需进行修正,采用向量计算的方式计算直线交点，并判定交点是否在线段上
    Poly1=Line_test{i}(:,1:2);
    PV1=Poly1;
    PV2=Rect_test(1:end-1,1:2);
    Vec1=Poly1(2:end,:)-Poly1(1:end-1,:);
    Vec2=Poly2(2:end,:)-Poly2(1:end-1,:);
    P_tempx=Poly1(1,1):(Poly1(2,1)-Poly1(1,1))/6:Poly1(2,1);
    P_tempy=Poly1(1,2):(Poly1(2,2)-Poly1(1,2))/6:Poly1(2,2);   %%裂缝简化为多个点组成集合
    In1 = inpolygon(Poly1(:,1),Poly1(:,2),PV2(:,1),PV2(:,2));  %%裂缝端点是否在计算边框内部
    In2 = inpolygon(P_tempx,P_tempy,PV2(:,1),PV2(:,2));        %%裂缝内点是否在计算边框内部
    if sum(In2)==0                 %裂缝段与计算方格无交集
        continue
    end
    if sum(In1)==2                 %裂缝段都在计算方框内
        Line_D{i}=Poly1;
        Line_index=cat(1,Line_index,i);
        Point_plot=cat(1,Point_plot,Poly1);
        Point_plot=cat(1,Point_plot,[nan,nan]);
        continue
    end
    if sum(In2)~=0                 %裂缝段部分在计算方框内
        Vtemx1=bsxfun(@minus,PV1(:,1),PV2(:,1)');
        Vtemy1=bsxfun(@minus,PV1(:,2),PV2(:,2)');
        Vtemx2=bsxfun(@minus,PV2(:,1),PV1(:,1)');
        Vtemy2=bsxfun(@minus,PV2(:,2),PV1(:,2)');
        t1=(diag(Vec1(:,2))*Vtemx1-diag(Vec1(:,1))*Vtemy1)./(Vec1(:,2)*Vec2(:,1)'-Vec1(:,1)*Vec2(:,2)');
        t2=(diag(Vec2(:,2))*Vtemx2-diag(Vec2(:,1))*Vtemy2)./(Vec2(:,2)*Vec1(:,1)'-Vec2(:,1)*Vec1(:,2)');
        Id=intersect(find(t2'<1&t2'>0),find(t1<1&t1>0));
        [~,J]=ind2sub(size(t1),Id);
        Point=PV2(J,:)+Vec2(J,:).*t1(Id)';
        if length(Id)==1
            In_temp = inpolygon(Poly1(:,1),Poly1(:,2),PV2(:,1),PV2(:,2));
            Point_temp=Poly1(In_temp,:);
            Point=cat(1,Point,Point_temp);
        end
        Line_D{i}=Point;
        Line_index=cat(1,Line_index,i);
        Point_plot=cat(1,Point_plot,Point);
        Point_plot=cat(1,Point_plot,[nan,nan]);        
    end
end
Line_D=Line_D(Line_index);
figure('Name','Counted Fracture in view plane','WindowStyle','docked')
hold on
plot(Poly2(:,1),Poly2(:,2),'linewidth',0.5,'color','b')
plot(Point_plot(:,1),Point_plot(:,2),'linewidth',1.5,'color','r')
axis square

%% 计算裂缝分形维数，先统计不同大小方框内，大于方框尺寸的裂缝数量
%%%采用盒计数法
t=log2(Range)+1;                   %叠代次数
N_dimension=zeros(t,1);
figure('Name','Count the number of fracture','WindowStyle','docked')
hold on
for E=1:t
    Ne=0;                          %累积覆盖信号的格子的总数
    cellsize=2^(E-1);              %每次的格子大小
    NumSeg=Range/cellsize;         %横轴划分成的段数
    for i=1:NumSeg
        for j=1:NumSeg             %由横轴第一个段起通过计算纵轴跨越的格子数累积N(e)
            Corner=[(i-1),(j-1)]*cellsize;
            Rect_grid=[[0 0 1 1]',[1 0 0 1]']*cellsize+Corner;
            patch(Rect_grid(:,1),Rect_grid(:,2),'r','facealpha',0.02)
            Num=Count_line(Rect_grid,Line_D,cellsize);
            Ne=Ne+Num;             %累加每一个格子数内所含裂缝数量（满足长度特征）
        end
    end
    N_dimension(E)=Ne;            %记录每尺度上的裂缝个数
end
axis square

%% 对log(N(e))和log(k/e)进行最小二乘的一次曲线拟合,斜率就是D
T=2.^(0:t-1);Te=T;
r=-diff(log2(N_dimension));     %去掉r超过2的野点数据
id=find(abs(r)>2);              
Ne=N_dimension;
Ne(id+1)=[];   Te(id+1)=[];
L0=2.^id;
P=polyfit(log2(Te),log2(Ne),1);       %一次曲线拟合返回斜率和截距
disp(['Dimension = ', num2str(abs(P(1)))])

x1 = linspace(0,5);
y1 = polyval(P,x1);
figure('Name','Caculated fractal dimension ','WindowStyle','docked')
plot(log2(T),log2(N_dimension),'o')
hold on
plot(x1,y1)
axis square
xlabel('log(L)')
ylabel('log(Num)')
set(gca, 'fontsize',15)

%% 视图上任意测线的RQD值计算
%%测线假定绕视图中心点旋转
PV2=Rect_test(1:end-1,1:2);
Poly2=Rect_test(:,1:2);
Vec2=Poly2(2:end,:)-Poly2(1:end-1,:);
Vtemx1=bsxfun(@minus,RQD_Center(:,1),PV2(:,1)');
Vtemy1=bsxfun(@minus,RQD_Center(:,2),PV2(:,2)');
t1=(diag(RQD_Vec(:,2))*Vtemx1-diag(RQD_Vec(:,1))*Vtemy1)./(RQD_Vec(:,2)*Vec2(:,1)'-RQD_Vec(:,1)*Vec2(:,2)');
Id=find(t1'<1&t1'>0);
[~,J]=ind2sub(size(t1),Id);
RQD_line=PV2(J,:)+Vec2(J,:).*t1(Id)';      %%测线与边框之间的交点

figure('Name','the ROD of a desired line','WindowStyle','docked')
hold on
plot(Poly2(:,1),Poly2(:,2),'linewidth',1,'color','k')
plot(Point_plot(:,1),Point_plot(:,2),'linewidth',0.5,'color','r')
plot(RQD_line(:,1),RQD_line(:,2),'linewidth',1.5,'color','b')
Point_RQD=RQD_line;
for i=1:size(Line_D,1)
    %线段点是否都在边框内部
    %是则无需修正，否则需进行修正,采用向量计算的方式计算直线交点，并判定交点是否在线段上
    PV1=Line_D{i}(:,1:2);
    PV2=RQD_line;
    Vec1=PV1(2:end,:)-PV1(1:end-1,:);
    Vec2=PV2(2:end,:)-PV2(1:end-1,:);
    Vtemx1=bsxfun(@minus,PV1(:,1),PV2(:,1)');
    Vtemy1=bsxfun(@minus,PV1(:,2),PV2(:,2)');
    Vtemx2=bsxfun(@minus,PV2(:,1),PV1(:,1)');
    Vtemy2=bsxfun(@minus,PV2(:,2),PV1(:,2)');
    t1=(diag(Vec1(:,2))*Vtemx1-diag(Vec1(:,1))*Vtemy1)./(Vec1(:,2)*Vec2(:,1)'-Vec1(:,1)*Vec2(:,2)');
    t2=(diag(Vec2(:,2))*Vtemx2-diag(Vec2(:,1))*Vtemy2)./(Vec2(:,2)*Vec1(:,1)'-Vec2(:,1)*Vec1(:,2)');
    Id=intersect(find(t2'<1&t2'>0),find(t1<1&t1>0));
    if isempty(Id)==1
        continue
    end
    [~,J]=ind2sub(size(t1),Id);
    P_temp=PV2(J,:)+Vec2.*t1(Id)';
    plot(P_temp(:,1),P_temp(:,2),'ro')
    Point_RQD=cat(1,P_temp,Point_RQD);
end
[~,Index_RQD]=unique(Point_RQD(:,1));
Point_RQD=Point_RQD(Index_RQD,:);
RQD_len=Point_RQD(2:end,:)-Point_RQD(1:end-1,:);
RQD_len=sqrt(sum(RQD_len.^2,2));
RQD=sum(RQD_len(RQD_len>RQD_judge))/sum(RQD_len);
disp(['RQD = ', num2str(RQD)])

plot(Point_RQD(:,1),Point_RQD(:,2),'ro')
for i=1:size(Point_RQD,1)-1
    plot(Point_RQD(i:i+1,1),Point_RQD(i:i+1,2),'linewidth',3)
end
axis square
xlim([0,Range])
ylim([0,Range])

toc

