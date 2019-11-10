for i=1:nnode
    plot(gcoord(i,1),gcoord(i,2),'.',gcoord(i,1),gcoord(i,2),'markersize',28);
    text(gcoord(i,1)+0.05,gcoord(i,2)+0.05,num2str(i),'color',[1 0 0],'fontsize',12);
end