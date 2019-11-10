for i=1:nel
    elab=['(',num2str(i),')'];
    if ElemTypes==1
        Exy=[gcoord(nodes(i,1),:);gcoord(nodes(i,2),:);gcoord(nodes(i,3),:);...
                gcoord(nodes(i,4),:)];
        text((Exy(1,1)+Exy(3,1))/2,(Exy(1,2)+Exy(3,2))/2,elab,...
             'color',[1 0 0],'fontsize',12);
    else
        Exy=[gcoord(nodes(i,1),:);gcoord(nodes(i,2),:);gcoord(nodes(i,3),:)];
        text((Exy(1,1)+Exy(2,1)+Exy(3,1))/3,(Exy(1,2)+Exy(2,2)+Exy(3,2))/3,...
            elab,'color',[1 0 0],'fontsize',12);
    end
end