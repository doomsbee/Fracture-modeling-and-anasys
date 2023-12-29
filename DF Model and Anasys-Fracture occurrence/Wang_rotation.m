Alfa=30/180*pi;   %倾角
Beta=60/180*pi;   %走向

%采用Ax+By+Cz+D=0的方式描述裂缝片
Aview=-sin(Alfa)*sin(Beta);     
Bview=sin(Alfa)*cos(Beta);
Cview=cos(Alfa);
Dview=3;                              %%A,B,C裂缝片的方向，D裂缝片位置
N_Vec=[Aview,Bview,Cview];            %%裂缝面的法向量
N_Vec0=[0,0,1];                       %%裂缝面初始法向量，默认与x,y平面平行
%%即旋转矩阵R将初始法向量N_Vec0，旋转为目标法向量N_Vec
%%采用轴转角方法，求解旋转矩阵R。
%%首先，求解旋转轴Ro_axis
%%其次，求解旋转角度theta(初始法向量N_Vec0与目标法向量N_Vec的夹角)
Ro_axis=cross(N_Vec0,N_Vec);                   %%轴向量
sin_theta=norm(Ro_axis);                       %%旋转角正弦
cos_theta=dot(N_Vec0,N_Vec);                   %%旋转角余弦
Ro_axis=Ro_axis/norm(Ro_axis);                 %%向量归一化
Rinve=[0,-Ro_axis(3),Ro_axis(2);Ro_axis(3),0,-Ro_axis(1);-Ro_axis(2),Ro_axis(1),0];          %%旋转轴向量反对称矩阵
RR=eye(3) + sin_theta*Rinve + (1- cos_theta)*Rinve*Rinve;                                    %%RR（列向量）与R（行向量）互为转置矩阵
N_Vec
N_Vec0*R
%由四元数构造旋转矩阵
theta = acos(dot(N_Vec,N_Vec0))/2;
q = [cos(theta) sin(theta)*Ro_axis];
R=[2*q(1).^2-1+2*q(2)^2  2*(q(2)*q(3)+q(1)*q(4)) 2*(q(2)*q(4)-q(1)*q(3));
    2*(q(2)*q(3)-q(1)*q(4)) 2*q(1)^2-1+2*q(3)^2 2*(q(3)*q(4)+q(1)*q(2));
    2*(q(2)*q(4)+q(1)*q(3)) 2*(q(3)*q(4)-q(1)*q(2)) 2*q(1)^2-1+2*q(4)^2];


%%案例验证，对中心点位于原点，方向与xy平面一致的矩形裂缝片，进行旋转。
Rect=[[-1 -1 1 1]',[1 -1 -1 1]',[0 0 0 0]']*20;
Rect2=Rect*R;
figure
hold on
fill3(Rect2(:,1),Rect2(:,2),Rect2(:,3),'r','facealpha',0.1)
plot3(Rect(:,1),Rect(:,2),Rect(:,3))

% 
% 
% 
% B=NVec1;                          %%裂缝面法向量,裂缝中心点[0,0,0]
% C=cross(A,NVec1);                 %%交线向量，其模等于两平面法向量夹角的正弦
% 
% Dist=dot(A,[0,0,0])-Dview;        %%裂缝中心点到视图面的垂直距离
% Length=Dist/norm(C);              %%裂缝中心点到交线的距离
