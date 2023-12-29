function [Disk_data,R,Fit]=Disk_gen(Len,x,y,z,Trend,Inclination)
%%������ǡ���������ѷ�Ƭ�ķ�������ͨ��������������ת���󣬽��ѷ�Ƭ������ת�������ƽ��[xs,ys,zs].��������������
%%%Len�ѷ�Ƭ���
%%%x,y,z�ѷ�Ƭ���ĵ�
%%%Trend����
%%%Inclination���
%%%Disk_dataΪ8�����������飬ǰ4��Ϊ�ѷ�Ƭ�������Ax+By+Cz+D=0���������ĸ�Ϊ���ĵ��������ѷ�Ƭ���
%%%R�ѷ�Ƭ��ת����
Disk_data=zeros(8,1);
Aview=-sin(Inclination)*sin(Trend);     
Bview=sin(Inclination)*cos(Trend);
Cview=cos(Inclination);               %%A,B,C�ѷ�Ƭ�ķ���D�ѷ�Ƭλ��
N_Vec=[Aview,Bview,Cview];            %%�ѷ���ķ�����
N_Vec0=[0,0,1];                       %%�ѷ����ʼ��������Ĭ����x,yƽ��ƽ��
%%����ת����R����ʼ������N_Vec0����תΪĿ�귨����N_Vec
%%������ת�Ƿ����������ת����R��
%%���ȣ������ת��Ro_axis
%%��Σ������ת�Ƕ�theta(��ʼ������N_Vec0��Ŀ�귨����N_Vec�ļн�)
Ro_axis=cross(N_Vec0,N_Vec);                   %%������
Ro_axis=Ro_axis/norm(Ro_axis);                 %%������һ��
%����Ԫ��������ת����
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

%%%���ɵ��ѷ죻ƬչʾЧ��
t=-pi:pi/2:pi;                      %%�ѷ�Ƭ�߽��Ϊ��ɢ��,����չʾ
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