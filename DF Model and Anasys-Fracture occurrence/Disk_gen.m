function [Disk_data,R,Fit]=Disk_gen(Len,x,y,z,Trend,Inclination)
%%根据倾角、倾向计算裂缝片的法向量，通过法向量计算旋转矩阵，将裂缝片进行旋转，最后再平移[xs,ys,zs].是内旋还是外旋
%%%Len裂缝片宽度
%%%x,y,z裂缝片中心点
%%%Trend倾向
%%%Inclination倾角
%%%Disk_data为8个参量的数组，前4个为裂缝片所在面的Ax+By+Cz+D=0参数，后四个为中心点坐标与裂缝片宽度
%%%R裂缝片旋转矩阵
Disk_data=zeros(8,1);
Aview=-sin(Inclination)*sin(Trend);     
Bview=sin(Inclination)*cos(Trend);
Cview=cos(Inclination);               %%A,B,C裂缝片的方向，D裂缝片位置
N_Vec=[Aview,Bview,Cview];            %%裂缝面的法向量
N_Vec0=[0,0,1];                       %%裂缝面初始法向量，默认与x,y平面平行
%%即旋转矩阵R将初始法向量N_Vec0，旋转为目标法向量N_Vec
%%采用轴转角方法，求解旋转矩阵R。
%%首先，求解旋转轴Ro_axis
%%其次，求解旋转角度theta(初始法向量N_Vec0与目标法向量N_Vec的夹角)
Ro_axis=cross(N_Vec0,N_Vec);                   %%轴向量
Ro_axis=Ro_axis/norm(Ro_axis);                 %%向量归一化
%由四元数构造旋转矩阵
theta = acos(dot(N_Vec,N_Vec0))/2;
q = [cos(theta) sin(theta)*Ro_axis];
R=[2*q(1).^2-1+2*q(2)^2  2*(q(2)*q(3)+q(1)*q(4)) 2*(q(2)*q(4)-q(1)*q(3));
    2*(q(2)*q(3)-q(1)*q(4)) 2*q(1)^2-1+2*q(3)^2 2*(q(3)*q(4)+q(1)*q(2));
    2*(q(2)*q(4)+q(1)*q(3)) 2*(q(3)*q(4)-q(1)*q(2)) 2*q(1)^2-1+2*q(4)^2];

Dview=dot([x,y,z],N_Vec);
Disk_data(1:3)=N_Vec;
Disk_data(4)=Dview;
Disk_data(5)=x;
Disk_data(6)=y;
Disk_data(7)=z;
Disk_data(8)=Len;

%%%生成的裂缝；片展示效果
t=-pi:pi/2:pi;                      %%裂缝片边界简化为离散点,进行展示
xc=Len*sin(t);
yc=Len*cos(t);
zc=0*t;
XYZ=[xc', yc', zc']*R;
xc=XYZ(:,1)+x;
yc=XYZ(:,2)+y;
zc=XYZ(:,3)+z;
Fit=[xc,yc,zc];
plot3(xc,yc,zc,'linewidth',1.5,'color','k')
hold on;
fill3(xc,yc,zc,'r','facealpha',0.3);
rotate3d on