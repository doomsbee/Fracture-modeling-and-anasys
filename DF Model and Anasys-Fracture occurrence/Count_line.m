function Num=Count_line(Rect_test,Line_test,Lsize)
Rect_test=cat(1,Rect_test,Rect_test(1,:));       %闭合计算边框
Poly2=Rect_test(:,1:2);
Point_plot=[];
Num=0;
for i=1:size(Line_test,1)
    %需进行进一步的筛选，以提高效率
    %线段点是否都在边框内部
    %是则无需修正，否则需进行修正,采用向量计算的方式计算直线交点，并判定交点是否在线段上
    Poly1=Line_test{i}(:,1:2);
    PV1=Poly1(1:end-1,1:2);
    PV2=Rect_test(1:end-1,1:2);
    Vec1=Poly1(2:end,:)-Poly1(1:end-1,:);
    Vec2=Poly2(2:end,:)-Poly2(1:end-1,:);
    D_num=ceil(norm(Vec1)/Lsize);
    P_tempx=Poly1(1,1):(Poly1(2,1)-Poly1(1,1))/D_num:Poly1(2,1);
    P_tempy=Poly1(1,2):(Poly1(2,2)-Poly1(1,2))/D_num:Poly1(2,2);   %%裂缝简化为多个点组成集合
    In1 = inpolygon(Poly1(:,1),Poly1(:,2),PV2(:,1),PV2(:,2));  %%裂缝端点是否在计算边框内部
    In2 = inpolygon(P_tempx,P_tempy,PV2(:,1),PV2(:,2));        %%裂缝内点是否在计算边框内部
    if sum(In2)==0                 %裂缝段与计算方格无交集
        continue
    end
    if sum(In1)==2                 %裂缝段所有点都在计算方框内
        if norm(Vec1)>Lsize
            Point_plot=cat(1,Point_plot,Poly1);
            Point_plot=cat(1,Point_plot,[nan,nan]);
            Num=Num+1;
        end
        continue
    end
    if sum(In2)~=0                 %裂缝段部分在计算方框内，需计算线段与方框的交点
        Vtemx1=bsxfun(@minus,PV1(:,1),PV2(:,1)');
        Vtemy1=bsxfun(@minus,PV1(:,2),PV2(:,2)');
        Vtemx2=bsxfun(@minus,PV2(:,1),PV1(:,1)');
        Vtemy2=bsxfun(@minus,PV2(:,2),PV1(:,2)');
        t1=(diag(Vec1(:,2))*Vtemx1-diag(Vec1(:,1))*Vtemy1)./(Vec1(:,2)*Vec2(:,1)'-Vec1(:,1)*Vec2(:,2)');
        t2=(diag(Vec2(:,2))*Vtemx2-diag(Vec2(:,1))*Vtemy2)./(Vec2(:,2)*Vec1(:,1)'-Vec2(:,1)*Vec1(:,2)');
        Id=intersect(find(t2'<1&t2'>0),find(t1<1&t1>0));
        [~,J]=ind2sub(size(t1),Id);
        Point=PV2(J,:)+Vec2(J,:).*t1(Id)';

        if length(Id)>=1                            %裂缝与计算方格相较
            In_temp =  inpolygon(Poly1(:,1),Poly1(:,2),PV2(:,1),PV2(:,2));
            Point_temp=Poly1(In_temp,:);
            Point=cat(1,Point,Point_temp);
%             Point_plot=cat(1,Point_plot,Point);
%             Point_plot=cat(1,Point_plot,[nan,nan]);
        end
        if norm(Point(2,:)-Point(1,:))>Lsize        %相较裂缝长度大于方框长度
            Point_plot=cat(1,Point_plot,Point);
            Point_plot=cat(1,Point_plot,[nan,nan]);
            plot(Poly1(:,1),Poly1(:,2))
            Num=Num+1;
        end
    end
end
if Num~=0
    plot(Point_plot(:,1),Point_plot(:,2),'linewidth',3)
end