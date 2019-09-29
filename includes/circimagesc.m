function circimagesc(v,markerSize)
    
    bg = repmat(min(min(v)),size(v,1),size(v,2));
    imagesc(bg);
    hold on

    [I,J] = ind2sub(size(v),1:length(v(:)));

    units = get(gca,'units');
   
    set(gca,'units','pixels')
    pos = get(gca,'position');
   
    markerSize = (pos(4)*pos(3))/(2.5*(size(v,1)*size(v,2)));
    
    scatter(J,I,markerSize,v(:),'filled');

    set(gca,'units',units);
    set(gca,'ydir','reverse');

    axis off;
end