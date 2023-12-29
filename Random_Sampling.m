function Sam= Random_Sampling(Rtype,Param1,Param2,Number1,Number2)

%Random_Sampling�����������ͣ�����������
%�˴���ʾ��ϸ˵��
switch nargin
    case 2
        %%%ָ���ֲ�
        Sam=random('exp',Param1,1,1);
    case 4
        %%%ָ���ֲ�
        Sam=random('exp',Param1,Number1,Number2);
    case 3
        switch Rtype
            case 0
                Sam=Param1(discreteinvrnd(Param2,1,1));                
            case 1          %%%��̫�ֲ�
                Sam=random('Normal',Param1,Param2,1,1);
            case 2          %%%���ȷֲ�
                Sam=random( 'Uniform',Param1,Param2,1,1);
            case 3          %%%������̫�ֲ�
                Sig=sqrt(log(1+Param2^2/Param1^2));      %ת��Ϊ��̫�ֲ�
                Mu=log(Param1)-1/2*Sig^2;              %���ݶ�����̫�ֲ�ͳ�Ƴ������ݣ�������Ӧ��̫�ֲ��Ĳ���
                Sam=random('logn',Mu,Sig,1,1);
            case 4          %%%Weibull�ֲ�
                Sam=random('wbl',Param1,Param2,1,1);
            otherwise        %%%ָ���ֲ�
                disp('Sorry the value is invalid');
        end
    case 5        
        switch Rtype
            case 0
                Sam=Param1(discreteinvrnd(Param2,Number1,Number2));                
            case 1          %%%��̫�ֲ�
                Sam=random('Normal',Param1,Param2,Number1,Number2);
            case 2          %%%���ȷֲ�
                Sam=random( 'Uniform',Param1,Param2,Number1,Number2);
            case 3          %%%������̫�ֲ�
                Sig=sqrt(log(1+Param2^2/Param1^2));      %ת��Ϊ��̫�ֲ�
                Mu=log(Param1)-1/2*Sig^2;              %���ݶ�����̫�ֲ�ͳ�Ƴ������ݣ�������Ӧ��̫�ֲ��Ĳ���
                Sam=random('logn',Mu,Sig,Number1,Number2);
            case 4          %%%Weibull�ֲ�
                Sam=random('wbl',Param1,Param2,Number1,Number2);
            otherwise        %%%ָ���ֲ�
                disp('Sorry the value is invalid');
        end
end
end


