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
%%%Disk_dataΪ8�����������飬ǰ4��Ϊ�ѷ�Ƭ�������Ax+By+Cz+D=0���������ĸ�Ϊ���ĵ��������ѷ�Ƭ���


clc
clear
close all
RGB_color
tic

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% A1.���������� %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Range=32;                              %%ģ�ͳߴ�,Ϊ2^N

View_type=0;                           %%��ͼȷ����ʽ��0Ϊ�����������1Ϊ����������
View_data=[-60,30];                     %%[��ǣ�����]
View_center=[Range,Range,Range]/2;     %%��ͼ������ĵ�  

RQD_angle=145;                          %%��ͼƽ�������ڲ���RQDֵ�Ĳ����߽Ƕȣ�0~180
RQD_Center=[1,1,1]/2*Range;            %%Ĭ�ϲ�����ת��ٶ�Ϊ��ͼƽ������ĵ㡣
RQD_judge=5;                           %%�����������趨��׼ֵ���淶�ٶ�Ϊ10

Num=50;                               %%�ѷ�Ƭ����
% % Ltype=1;                               %%�ѷ�Ƭ�볤�ֲ����ͣ����ȷֲ������ޣ����ޣ���̬�ֲ�����ֵ�����
% % Lpara1=5;                              %%�ѷ�Ƭ�볤�ֲ�����1
% % Lpara2=5;                              %%�ѷ�Ƭ�볤�ֲ�����2

% % Itype=1;                               %%�ѷ�Ƭ��Ƿֲ�����: ��̬�ֲ�����ֵ�����
% % Ipara1=pi/3;                           %%�ѷ�Ƭ��Ƿֲ�����1
% % Ipara2=pi/50;                          %%�ѷ�Ƭ��Ƿֲ�����2
% % 
% % Ttype=1;                               %%�ѷ�Ƭ����ֲ�����
% % Tpara1=pi/6;                           %%�ѷ�Ƭ����ֲ�����1
% % Tpara2=pi/50;                          %%�ѷ�Ƭ����ֲ�����2

Radius=repmat(5,Num,1);
Trend=repmat(5*pi/6,Num,1);
Inclination=repmat(pi/2,Num,1);

RQD_Vec=[cos(RQD_angle/180*pi),sin(RQD_angle/180*pi)];   %RQD��������
Rect=[[-1 -1 1 1]',[1 -1 -1 1]',[0 0 0 0]']*Range/2;     %%��ʼ��ͼ����
Rect_test=Rect+View_center; 

%% �����ѷ�Ƭ������ֲ���ʽ���ֲ������������ѷ�Ƭ
%%%������ͼ����������ͼƽ��
Center_point=rand(Num,3)*Range;                           %%�ѷ�Ƭ���ĵ��������
% % Radius= Random_Sampling(Ltype,Lpara1,Lpara2,Num,1);        %%�ѷ�Ƭ�İ볤����ǡ������������
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
    %%%�����ѷ�Ƭ
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

%% ��ʼ��ͼƽ����xyƽ��ƽ�У���ͼ���ĵ�ΪView_center
%%%������ͼ����������ͼƽ��
if View_type==0                                     %%ͨ��������ǣ�����ȷ����ͼ��ķ�����
    View_alfa=View_data(1)*pi/180;                  %%��ͼ������
    View_beta=View_data(2)*pi/180;                  %%��ͼ�������
    Aview=-sin(View_alfa)*sin(View_beta);
    Bview=sin(View_alfa)*cos(View_beta);
    Cview=cos(View_alfa);                           %%A,B,C�ѷ�Ƭ�ķ���D�ѷ�Ƭλ��
else                                                %%ֱ��������ͼ��ķ�����
    Aview=View_data(1);
    Bview=View_data(2);
    Cview=View_data(3);
end
N_Vec=[Aview,Bview,Cview];                          %%�ѷ���ķ�����
Dview=dot(N_Vec,View_center);                       %%�ѷ�Ƭ��λ��
N_Vec0=[0,0,1];                                     %%�ѷ����ʼ��������Ĭ����x,yƽ��ƽ��
Ro_axis=cross(N_Vec0,N_Vec);                        %%������
Ro_axis=Ro_axis/norm(Ro_axis);                      %%������һ��
theta = acos(dot(N_Vec,N_Vec0))/2;
q = [cos(theta) sin(theta)*Ro_axis];
R=[2*q(1).^2-1+2*q(2)^2  2*(q(2)*q(3)+q(1)*q(4)) 2*(q(2)*q(4)-q(1)*q(3));
    2*(q(2)*q(3)-q(1)*q(4)) 2*q(1)^2-1+2*q(3)^2 2*(q(3)*q(4)+q(1)*q(2));
    2*(q(2)*q(4)+q(1)*q(3)) 2*(q(3)*q(4)-q(1)*q(2)) 2*q(1)^2-1+2*q(4)^2];
Rect2=Rect*R+View_center;                           %%��ת����ѷ췽��


%% ������ͼ�����ϵ��ѷ���ʾ��ͳ��
%%%��������Բ�ĵ㵽�ѷ�Ƭ����ͼ�潻�ߵľ���
figure('Name','Fracture in View plane','WindowStyle','docked')
hold on
A=[Aview,Bview,Cview];                             %%��ͼ����ķ�����
In_num=0;
Index=[];
Line_view={};                                      %%�ռ�������ͼ�������γɵ��ѷ�
Line_test={};                                      %%ת��Ϊƽ���������γɵ��ѷ�
for i=1:Num
    B=Fracture(i,1:3);                             %%�ѷ�Ƭ������
    C=cross(A,B);                                  %%CΪ���ߵķ���
    sin_temp=norm(C);                              %%C��ģΪ�ѷ�Ƭ������B����ͼ�淨����A�ļн�����
    Dplane=dot(C,Fracture(i,5:7));                 %%Dplane���㴹���λ��
    Dist=dot(A,Fracture(i,5:7))-Dview;             %%�ѷ�Ƭ���ĵ㵽��ͼ��ķ������
    Length=Dist/norm(C);                           %%�ѷ�Ƭ���ĵ㵽�ѷ�Ƭ����ͼ�潻�ߵľ���
    if abs(Length) <  Radius(i)                    %%�ѷ�Ƭ����ͼ����ڽ��ߣ������ཻ
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

        plot3(Line(:,1),Line(:,2),Line(:,3),'linewidth',1.5,'color','r')               %%�ռ��ѷ�Ƭ����ͼ����ཻ�ѷ�
        fill3(Disk_Fit{i}(:,1),Disk_Fit{i}(:,2),Disk_Fit{i}(:,3),'b','facealpha',0.1)  %%�ռ��ѷ�Ƭ
        plot3(Line_R(:,1),Line_R(:,2),Line_R(:,3),'linewidth',1.5,'color','b')         %%ת����xyƽ���µ��ཻ�ѷ�
    end
end

fill3(Rect2(:,1),Rect2(:,2),Rect2(:,3),'g','facealpha',0.1)
plotcube1([Range Range Range],[0 0 0],0.02,Color);
axis square


%% ������ͼ����߿򣬻�ñ߿��ڵ��ѷ�
Rect_test=cat(1,Rect_test,Rect_test(1,:));       %�պϼ���߿�
Poly2=Rect_test(:,1:2);
Line_D=Line_test;
Line_index=[];
Point_plot=[];
for i=1:In_num
    %�߶ε��Ƿ��ڱ߿��ڲ�
    %���������������������������,������������ķ�ʽ����ֱ�߽��㣬���ж������Ƿ����߶���
    Poly1=Line_test{i}(:,1:2);
    PV1=Poly1;
    PV2=Rect_test(1:end-1,1:2);
    Vec1=Poly1(2:end,:)-Poly1(1:end-1,:);
    Vec2=Poly2(2:end,:)-Poly2(1:end-1,:);
    P_tempx=Poly1(1,1):(Poly1(2,1)-Poly1(1,1))/6:Poly1(2,1);
    P_tempy=Poly1(1,2):(Poly1(2,2)-Poly1(1,2))/6:Poly1(2,2);   %%�ѷ��Ϊ�������ɼ���
    In1 = inpolygon(Poly1(:,1),Poly1(:,2),PV2(:,1),PV2(:,2));  %%�ѷ�˵��Ƿ��ڼ���߿��ڲ�
    In2 = inpolygon(P_tempx,P_tempy,PV2(:,1),PV2(:,2));        %%�ѷ��ڵ��Ƿ��ڼ���߿��ڲ�
    if sum(In2)==0                 %�ѷ������㷽���޽���
        continue
    end
    if sum(In1)==2                 %�ѷ�ζ��ڼ��㷽����
        Line_D{i}=Poly1;
        Line_index=cat(1,Line_index,i);
        Point_plot=cat(1,Point_plot,Poly1);
        Point_plot=cat(1,Point_plot,[nan,nan]);
        continue
    end
    if sum(In2)~=0                 %�ѷ�β����ڼ��㷽����
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

%% �����ѷ����ά������ͳ�Ʋ�ͬ��С�����ڣ����ڷ���ߴ���ѷ�����
%%%���úм�����
t=log2(Range)+1;                   %��������
N_dimension=zeros(t,1);
figure('Name','Count the number of fracture','WindowStyle','docked')
hold on
for E=1:t
    Ne=0;                          %�ۻ������źŵĸ��ӵ�����
    cellsize=2^(E-1);              %ÿ�εĸ��Ӵ�С
    NumSeg=Range/cellsize;         %���Ữ�ֳɵĶ���
    for i=1:NumSeg
        for j=1:NumSeg             %�ɺ����һ������ͨ�����������Խ�ĸ������ۻ�N(e)
            Corner=[(i-1),(j-1)]*cellsize;
            Rect_grid=[[0 0 1 1]',[1 0 0 1]']*cellsize+Corner;
            patch(Rect_grid(:,1),Rect_grid(:,2),'r','facealpha',0.02)
            Num=Count_line(Rect_grid,Line_D,cellsize);
            Ne=Ne+Num;             %�ۼ�ÿһ���������������ѷ����������㳤��������
        end
    end
    N_dimension(E)=Ne;            %��¼ÿ�߶��ϵ��ѷ����
end
axis square

%% ��log(N(e))��log(k/e)������С���˵�һ���������,б�ʾ���D
T=2.^(0:t-1);Te=T;
r=-diff(log2(N_dimension));     %ȥ��r����2��Ұ������
id=find(abs(r)>2);              
Ne=N_dimension;
Ne(id+1)=[];   Te(id+1)=[];
L0=2.^id;
P=polyfit(log2(Te),log2(Ne),1);       %һ��������Ϸ���б�ʺͽؾ�
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

%% ��ͼ��������ߵ�RQDֵ����
%%���߼ٶ�����ͼ���ĵ���ת
PV2=Rect_test(1:end-1,1:2);
Poly2=Rect_test(:,1:2);
Vec2=Poly2(2:end,:)-Poly2(1:end-1,:);
Vtemx1=bsxfun(@minus,RQD_Center(:,1),PV2(:,1)');
Vtemy1=bsxfun(@minus,RQD_Center(:,2),PV2(:,2)');
t1=(diag(RQD_Vec(:,2))*Vtemx1-diag(RQD_Vec(:,1))*Vtemy1)./(RQD_Vec(:,2)*Vec2(:,1)'-RQD_Vec(:,1)*Vec2(:,2)');
Id=find(t1'<1&t1'>0);
[~,J]=ind2sub(size(t1),Id);
RQD_line=PV2(J,:)+Vec2(J,:).*t1(Id)';      %%������߿�֮��Ľ���

figure('Name','the ROD of a desired line','WindowStyle','docked')
hold on
plot(Poly2(:,1),Poly2(:,2),'linewidth',1,'color','k')
plot(Point_plot(:,1),Point_plot(:,2),'linewidth',0.5,'color','r')
plot(RQD_line(:,1),RQD_line(:,2),'linewidth',1.5,'color','b')
Point_RQD=RQD_line;
for i=1:size(Line_D,1)
    %�߶ε��Ƿ��ڱ߿��ڲ�
    %���������������������������,������������ķ�ʽ����ֱ�߽��㣬���ж������Ƿ����߶���
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

