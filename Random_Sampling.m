function Sam= Random_Sampling(Rtype,Param1,Param2,Number1,Number2)

%Random_Sampling根据输入类型，输出随机变量
%此处显示详细说明
switch nargin
    case 2
        %%%指数分布
        Sam=random('exp',Param1,1,1);
    case 4
        %%%指数分布
        Sam=random('exp',Param1,Number1,Number2);
    case 3
        switch Rtype
            case 0
                Sam=Param1(discreteinvrnd(Param2,1,1));                
            case 1          %%%正太分布
                Sam=random('Normal',Param1,Param2,1,1);
            case 2          %%%均匀分布
                Sam=random( 'Uniform',Param1,Param2,1,1);
            case 3          %%%对数正太分布
                Sig=sqrt(log(1+Param2^2/Param1^2));      %转换为正太分布
                Mu=log(Param1)-1/2*Sig^2;              %根据对数正太分布统计出的数据，计算相应正太分布的参数
                Sam=random('logn',Mu,Sig,1,1);
            case 4          %%%Weibull分布
                Sam=random('wbl',Param1,Param2,1,1);
            otherwise        %%%指数分布
                disp('Sorry the value is invalid');
        end
    case 5        
        switch Rtype
            case 0
                Sam=Param1(discreteinvrnd(Param2,Number1,Number2));                
            case 1          %%%正太分布
                Sam=random('Normal',Param1,Param2,Number1,Number2);
            case 2          %%%均匀分布
                Sam=random( 'Uniform',Param1,Param2,Number1,Number2);
            case 3          %%%对数正太分布
                Sig=sqrt(log(1+Param2^2/Param1^2));      %转换为正太分布
                Mu=log(Param1)-1/2*Sig^2;              %根据对数正太分布统计出的数据，计算相应正太分布的参数
                Sam=random('logn',Mu,Sig,Number1,Number2);
            case 4          %%%Weibull分布
                Sam=random('wbl',Param1,Param2,Number1,Number2);
            otherwise        %%%指数分布
                disp('Sorry the value is invalid');
        end
end
end


