% Ryan Calmus, 2019
% Laboratory of Comparative Neuropsychology
% Newcastle University

function VSBIND_Figure3(seed)

    disp('Figure 3 diagrams (3 windows) rendering. Wait for the figure to render, and then resize the figure window slowly until the aspect ratio is correct.');

    vectorSize = 64;
    usecircles = true;

    randState = rng;
    if nargin < 1  % choose the prettiest vectors of certain sizes to plot by setting the random number generator seed
        if vectorSize == 64
            seed = 57;
        elseif vectorSize == 1024
            seed = 15;
        end
    end

    rng(seed);

    % NB Key element names are all followed by "_ " in the title because this forces
    % all titles to be rendered using the subscript-titled sizes, in all
    % figures, and thus aspect ratios match when all figures are rendered.
    
    A1 = SymbolicVector(vectorSize);
    A1.name = '\bf{ A_ }';

    X1 = SymbolicVector(vectorSize);
    X1.name = '\bf{ X_ }';

    X2 = SymbolicVector(vectorSize);
    X2.name = '\bf{ X_ }';

    A2 = SymbolicVector(vectorSize);
    A2.name = '\bf{ A2_ }';

    B2 = SymbolicVector(vectorSize);
    B2.name = '\bf{ B2_ }';

    B1 = SymbolicVector(vectorSize);
    B1.name = '\bf{ B_ }';
    B1.visualisationLineWidth = 10;

    X2.visualisationLineWidth = B1.visualisationLineWidth*0.75;
    A1.visualisationLineWidth = X2.visualisationLineWidth;
    X1.visualisationLineWidth = X2.visualisationLineWidth*0.75;
    A2.visualisationLineWidth = X1.visualisationLineWidth;
    B2.visualisationLineWidth = X2.visualisationLineWidth;

    Primary = SymbolicVector(vectorSize);
    Primary.name = '\color{black}\bf{1^{\circ}_I}';
    Primary.visualisationColor = [173 185 202]/255;

    Secondary = SymbolicVector(vectorSize);
    Secondary.name = '\color{black}\bf{2^{\circ}_I}';
    Secondary.visualisationColor = [173 185 202]/255;

    Context1 = SymbolicVector(vectorSize);
    Context1.name = '\color{black}\bf{Context 1}';

    Context2 = SymbolicVector(vectorSize);
    Context2.name = '\color{black}\bf{Context 2}';

    Primary1 = Primary;
    Primary1.name = '\color{black}\bf{1^{\circ}_I}';
    Primary1.visualisationColor = [173 185 202]/255;

    Secondary1 = Secondary;
    Secondary1.name = '\color{black}\bf{2^{\circ}_I}';
    Secondary1.visualisationColor = [173 185 202]/255;

    Primary2 = Primary*Context2;
    Primary2.name = '\color{black}\bf{1^{\circ}_C}';
    Primary2.visualisationColor = [173 185 202]/255;

    Secondary2 = Secondary*Context2;
    Secondary2.name = '\color{black}\bf{2^{\circ}_C}';
    Secondary2.visualisationColor = [173 185 202]/255;

    Item1 = Primary*A1;
    Item1.name = '1^ary*A1';
    Item1.visualisationColor = 'k';

    Item2 = Secondary*B1;
    Item2.name = '2^ary*A1';
    Item2.visualisationColor = 'k';

    Item1b = Primary*A2;
    Item1b.name = '1^ary*A2';
    Item1b.visualisationColor = 'k';

    Item2b = Secondary*B2;
    Item2b.name = '2^ary*A2';
    Item2b.visualisationColor = 'k';

    Dependencies = Item1+Item2;
    Dependencies.name = '(1^{\circ}*A1) + (2^{\circ}*B1)';
    Dependencies.visualisationColor = 'k';

    Dependencies2 = Item1b+Item2b;
    Dependencies2.name = '(1^{\circ}*A2) + (2^{\circ}*B2)';
    Dependencies2.visualisationColor = 'k';

    Rule1 = Primary2*Dependencies;
    Rule1.name = '1^{\circ}*Dependency 1';
    Rule1.visualisationColor = 'k';

    Rule2 = Secondary2*Dependencies2;
    Rule2.name = '1^{\circ}*Dependency 2';
    Rule2.visualisationColor = 'k';

    Nested = Rule1+Rule2;
    Nested.name = '(1^{\circ}*Dependency 1) + (2^{\circ}*Dependency 2)';
    Nested.visualisationColor = 'k';

    %% Adjacent dependency diagram

    h = scalingFigure([8,7]);

    A1Plot = multispy(A1,true,usecircles,8,7,52);
    B1Plot = multispy(B1,true,usecircles,8,7,54);

    FirstPlot = multispy(Primary,true,usecircles,8,7,37);
    LastPlot  = multispy(Secondary,true,usecircles,8,7,41);

    Item1Plot = multispy(Item1,true,usecircles,8,7,38,'*');
    Item2Plot = multispy(Item2,true,usecircles,8,7,40,'*');

    DepPlot = multispy(Dependencies,true,usecircles,8,7,32,'+');

    arrowBetweenSubplots(A1Plot,Item1Plot,'yx');
    arrowBetweenSubplots(FirstPlot,Item1Plot,'yx');
    arrowBetweenSubplots(B1Plot,Item2Plot,'yx');
    arrowBetweenSubplots(LastPlot,Item2Plot,'yx');

    arrowBetweenSubplots(Item1Plot,DepPlot,'yx');
    arrowBetweenSubplots(Item2Plot,DepPlot,'yx');

    %% Nonadjacency dependency diagram
    A1.visualisationLineWidth = X1.visualisationLineWidth*0.75;

    h = scalingFigure([8,7]);

    A1Plot = multispy(A1,true,usecircles,8,7,51);
    B1Plot = multispy(B1,true,usecircles,8,7,55);
    X1Plot = multispy(X1,true,usecircles,8,7,52);
    X2Plot = multispy(X2,true,usecircles,8,7,54);

    FirstPlot = multispy(Primary,true,usecircles,8,7,36);
    LastPlot  = multispy(Secondary,true,usecircles,8,7,42);

    Item1Plot = multispy(Item1,true,usecircles,8,7,37,'*');
    Item2Plot = multispy(Item2,true,usecircles,8,7,41,'*');

    DepPlot = multispy(Dependencies,true,usecircles,8,7,32,'+');

    arrowBetweenSubplots(A1Plot,Item1Plot,'yx');
    arrowBetweenSubplots(FirstPlot,Item1Plot,'yx');
    arrowBetweenSubplots(B1Plot,Item2Plot,'yx');
    arrowBetweenSubplots(LastPlot,Item2Plot,'yx');

    arrowBetweenSubplots(Item1Plot,DepPlot,'yx');
    arrowBetweenSubplots(Item2Plot,DepPlot,'yx');

    %% Nested dependency diagram

    Primary.name = '\color{black}\bf{1^{\circ}_I}';
    Secondary.name = '\color{black}\bf{2^{\circ}_I}';

    Dependencies.visualisationColor = [255 192 0]/255;
    Item1.visualisationColor = Dependencies.visualisationColor;
    Item2.visualisationColor = Dependencies.visualisationColor;
    Dependencies2.visualisationColor = [0 176 240]/255;
    Item1b.visualisationColor = Dependencies2.visualisationColor;
    Item2b.visualisationColor = Dependencies2.visualisationColor;

    A1.visualisationLineWidth = X1.visualisationLineWidth*0.75;

    h = scalingFigure([8,7]);

    A1.name = '\bf{A1_ }';
    B1.name = '\bf{B1_ }';

    A1Plot = multispy(A1,true,usecircles,8,7,51);
    B1Plot = multispy(B1,true,usecircles,8,7,55);
    A2Plot = multispy(A2,true,usecircles,8,7,52);
    B2Plot = multispy(B2,true,usecircles,8,7,54);

    FirstPlot = multispy(Primary,true,usecircles,8,7,36);
    LastPlot  = multispy(Secondary,true,usecircles,8,7,42);

    Item1Plot = multispy(Item1,true,usecircles,8,7,37,'*');
    Item2Plot = multispy(Item2,true,usecircles,8,7,41,'*');

    Item1bPlot = multispy(Item1b,true,usecircles,8,7,31,'*');
    Item2bPlot = multispy(Item2b,true,usecircles,8,7,33,'*');

    DepPlot = multispy(Dependencies,true,usecircles,8,7,17,'+');
    Dep2Plot = multispy(Dependencies2,true,usecircles,8,7,26,'+');

    Rule1Plot = multispy(Rule1,true,usecircles,8,7,10,'*');
    Rule2Plot = multispy(Rule2,true,usecircles,8,7,12,'*');

    First2Plot = multispy(Primary2,true,usecircles,8,7,8);
    Last2Plot  = multispy(Secondary2,true,usecircles,8,7,14);

    NestedPlot = multispy(Nested,true,usecircles,8,7,4,'+');

    arrowBetweenSubplots(A1Plot,Item1Plot,'yx');
    arrowBetweenSubplots(FirstPlot,Item1Plot,'yx');
    arrowBetweenSubplots(B1Plot,Item2Plot,'yx');
    arrowBetweenSubplots(LastPlot,Item2Plot,'yx');

    arrowBetweenSubplots(A2Plot,Item1bPlot,'yx');
    arrowBetweenSubplots(FirstPlot,Item1bPlot,'yx');
    arrowBetweenSubplots(B2Plot,Item2bPlot,'yx');
    arrowBetweenSubplots(LastPlot,Item2bPlot,'yx');

    arrowBetweenSubplots(Item1Plot,DepPlot,'yx','post');
    arrowBetweenSubplots(Item2Plot,DepPlot,'yx','post');

    arrowBetweenSubplots(Item1bPlot,Dep2Plot,'yx','post');
    arrowBetweenSubplots(Item2bPlot,Dep2Plot,'yx','post');

    arrowBetweenSubplots(DepPlot,Rule1Plot,'yx');
    arrowBetweenSubplots(First2Plot,Rule1Plot,'yx');

    arrowBetweenSubplots(Dep2Plot,Rule2Plot,'yx');
    arrowBetweenSubplots(Last2Plot,Rule2Plot,'yx');

    arrowBetweenSubplots(Rule1Plot,NestedPlot,'yx');
    arrowBetweenSubplots(Rule2Plot,NestedPlot,'yx');

    rng(randState);
end

% Ryan Calmus, Petkov Lab, 2019 - Draw diagonal or elbow arrows to connect subplots
% sensibly.
function arrowBetweenSubplots(startPlotHandle,finishPlotHandle,elbowType,prePostColorMatch,varargin)
    startX = startPlotHandle.Position(1)+(startPlotHandle.Position(3))/2;
    finishX = finishPlotHandle.Position(1)+(finishPlotHandle.Position(3))/2;
    startY = startPlotHandle.Position(2)+(startPlotHandle.Position(4))/2;
    finishY = finishPlotHandle.Position(2)+(finishPlotHandle.Position(4))/2;

    yOffset = -0.01;
    finishY = finishY + yOffset;
    startY = startY + yOffset;

    if ~isempty(get(get(startPlotHandle,'Title'),'String'))
        startTitleHeight = 0.02;
    else
        startTitleHeight = -0.02;
    end

    if ~isempty(get(get(finishPlotHandle,'Title'),'String'))
        finishTitleHeight = 0.02;
    else
        finishTitleHeight = -0.02;
    end

    if nargin < 3 || isempty(elbowType)
        elbowType = 'xy';
    end

    if nargin < 4 || isempty(prePostColorMatch)
        prePostColorMatch = 'match pre';
    end

    if strcmpi(elbowType,'xy')
        if finishX>startX
            startX = startX+(startPlotHandle.Position(3)/2);
        elseif finishX<startX
            startX = startX-(startPlotHandle.Position(3)/2);
        end

        if startY>finishY
            finishY = finishY+(finishPlotHandle.Position(4)/2)+finishTitleHeight;
        elseif finishY>startY
            finishY = finishY-(finishPlotHandle.Position(4)/2);
        end
    else
        if finishX>startX
            finishX = finishX-(finishPlotHandle.Position(3)/2);
        elseif finishX<startX
            finishX = finishX+(finishPlotHandle.Position(3)/2);
        end

        if startY>finishY
            startY = startY-(startPlotHandle.Position(4)/2);
        elseif finishY>startY
            startY = startY+(startPlotHandle.Position(4)/2)+startTitleHeight;
        end
    end

    if strcmpi(prePostColorMatch,'match pre') || strcmpi(prePostColorMatch,'pre')
        colorToMatch = startPlotHandle.XColor;
    elseif strcmpi(prePostColorMatch,'match post') || strcmpi(prePostColorMatch,'post')
        colorToMatch = finishPlotHandle.XColor;
    end

    headsize = (10/3)*startPlotHandle.LineWidth;
    extraArgsLine = [{'LineWidth',startPlotHandle.LineWidth,'Color',colorToMatch},varargin];
    extraArgsArrow = [{'LineWidth',startPlotHandle.LineWidth,'HeadWidth',headsize,'HeadLength',headsize,'Color',colorToMatch},varargin];
    capExtraArgs = [{'LineWidth',1,'LineStyle','none','FaceColor',colorToMatch},varargin];

    if strcmpi(elbowType,'diagonal')
        annotation2(startPlotHandle,'arrow', [startX finishX],[startY finishY],extraArgsArrow{:});
    else
        capRadius = 0.0007*startPlotHandle.LineWidth; % prettify the elbow of the arrow
        if finishX-startX ~= 0
            if finishY-startY ~= 0
                if strcmpi(elbowType,'xy')
                    annotation2(startPlotHandle,'line', [startX finishX],[startY startY],extraArgsLine{:});
                    circle2(finishX,startY,capRadius,capExtraArgs{:});  % to add corner
                    annotation2(startPlotHandle,'arrow', [finishX finishX],[startY finishY],extraArgsArrow{:});
                else
                    annotation2(startPlotHandle,'line', [startX startX],[startY finishY],extraArgsLine{:});
                    circle2(startX,finishY,capRadius,capExtraArgs{:});  % to add corner
                    annotation2(startPlotHandle,'arrow', [startX finishX],[finishY finishY],extraArgsArrow{:});
                end
            else
                if startX>finishX
                    startX = startX-(startPlotHandle.Position(3)/2);
                    annotation2(startPlotHandle,'arrow', [startX finishX],[startY finishY],extraArgsArrow{:});
                elseif startX<finishX
                    startX = startX+(startPlotHandle.Position(3)/2);
                    annotation2(startPlotHandle,'arrow', [startX finishX],[startY finishY],extraArgsArrow{:});
                else
                end
            end
        else
            if startY>finishY
                startY = startY-(startPlotHandle.Position(4)/2);
                annotation2(startPlotHandle,'arrow', [startX finishX],[startY finishY],extraArgsArrow{:});
            elseif startY<finishY
                finishY = finishY-(finishPlotHandle.Position(4)/2);
                annotation2(startPlotHandle,'arrow', [startX finishX],[startY finishY],extraArgsArrow{:});
            else
            end
        end

    end

    uistack(startPlotHandle,'top');
    uistack(finishPlotHandle,'top');
end

function h = circle2(x,y,r,varargin)
    d = r*2;
    px = x-r;
    py = y-r;

    scaleFactors = pbaspect;

    h = annotation2(gca,'Ellipse',[px py d d*(scaleFactors(2)/scaleFactors(1))],varargin{:});
end

function h = annotation2(startPlotHandle,anType,Xs,Ys,varargin)
    h = annotation(anType, Xs, Ys, varargin{:});
    return
    if strcmpi(anType,'arrow')
        h = annotation(anType, Xs, Ys, varargin{:});
        propsToRemove = {'headwidth','headlength'};
        for iRemoveProps = 1:length(propsToRemove)
            propInd = find(cellfun(@(x) strcmpi(x,propsToRemove{iRemoveProps}),varargin));
            try
                varargin([propInd,propInd+1]) = [];
            catch
            end
        end

    elseif strcmpi(anType,'line')
        datX = [];
        datY = [];

        [datX(1),datY(1)] = figureSpaceToDataSpaceCoords(startPlotHandle,Xs(1),Ys(1));
        [datX(2),datY(2)] = figureSpaceToDataSpaceCoords(startPlotHandle,Xs(2),Ys(2));

        h = line(startPlotHandle,datX,datY,varargin{:});
    else
        h = annotation(anType, Xs, Ys, varargin{:});
    end
end

% Ryan Calmus, Petkov Lab, 2019 - Create dynamically rendering figure
% capable of scaling diagrammatic components with window size.

function h = scalingFigure(subplotDims)
    if nargin < 1
        subplotDims = [1,1];
    end

    function Resize_clbk(hObject,Event)
        try
            if hObject.updatingSizes
                return
            else
                hObject.updatingSizes = true;
            end
            allHandles = findall(hObject,'-property', 'Position');
            for i = 1:length(allHandles)
                try
                    set(allHandles(i),'Units','Normalized');            
                catch
                end
            end
    
            if hObject.firstResize
                hObject.firstResize = false;

                graphicHandles = findall(hObject,'-property', 'LineWidth');
                arrayfun(@(x) set(x,'LineWidth',get(x,'LineWidth')*hObject.annotationLineScaleFactor),graphicHandles);
            end

            graphicHandles = findall(hObject,'-property', 'LineWidth');
            scale = hObject.Position(3)/hObject.origPos(3);

            arrayfun(@(x) set(x,'LineWidth',get(x,'LineWidth')*scale),graphicHandles);

            graphicHandles = findall(hObject,'-property', 'HeadWidth');
            arrayfun(@(x) set(x,'HeadWidth',get(x,'HeadWidth')*scale),graphicHandles);        
            arrayfun(@(x) set(x,'HeadLength',get(x,'HeadLength')*scale),graphicHandles);

            graphicHandles = findall(hObject,'-property', 'FontSize');
            arrayfun(@(x) set(x,'FontSize',get(x,'FontSize')*scale),graphicHandles);

            ellipseHandles = findall(hObject,'Type', 'Ellipse','-or','Type','TextBox');
            for iEll = 1:length(ellipseHandles)
                
                curPos = ellipseHandles(iEll).Position;
                newElPos = curPos;

                newElPos(3) = curPos(3) * scale;
                newElPos(1) = (curPos(1)+curPos(3)/2) - (ellipseHandles(iEll).Position(3)/2);
                newElPos(4) = curPos(4) * scale;
                newElPos(2) = (curPos(2)-curPos(4)/2) + (ellipseHandles(iEll).Position(4)/2);

                ellipseHandles(iEll).Position = newElPos;
            end

            graphicHandles = findall(hObject,'Type', 'scatter');
            arrayfun(@(x) set(x,'SizeData',get(x,'SizeData')*scale^2),graphicHandles);

            hObject.origPos = hObject.Position;
            hObject.updatingSizes = false;
        catch me
            me
            printstruct(me.stack)
            hObject.origPos = hObject.Position;
            hObject.updatingSizes = false;
        end
    end

    h = figure;
    set(h,'Units','Normalized');
    addprop(h,'origPos');
    
     try
         w = warning();
         warning('off','all');
         jFrame = get(h,'JavaFrame');
         pause(0.3);
         set(jFrame,'Maximized',true);
         warning(w);
     catch
        h.OuterPosition = [0 0 1 1];
     end
    drawnow();

    set(h,'Units','Pixels');
    newPos = h.OuterPosition;
    newPos(3) = newPos(4);
    h.OuterPosition = newPos;
    
    drawnow();

    anAx = subplot(subplotDims(1),subplotDims(2),1);
    drawnow();
    as = pbaspect(anAx);

    addprop(h,'annotationLineScaleFactor');
    anAn = rectangle(anAx,'Position',[0 0 0.1 0.1],'LineWidth',7);

    set(anAx,'Units','pixels');
    set(h,'Units','pixels');
    anPos = anAn.Position;
    axlim = axis(anAx);
    axwidth = diff(axlim(1:2));
    axpos = get(anAx,'Position');
    set(h,'Units','Normalized');

    h.annotationLineScaleFactor = 0.0008*axpos(3)/axwidth;

    delete(anAn);

    delete(anAx);

    drawnow();

    set(h,'Units','pixels');
    newPos = h.Position;
    newPos(3) = newPos(4)*as(2)/as(1);
    h.Position = newPos;

    set(h,'Units','Normalized');

    drawnow();

    h.origPos = h.Position;

    addprop(h,'updatingSizes');
    h.updatingSizes = false;

    addprop(h,'firstResize');
    h.firstResize = true;

    set(h,'SizeChangedFcn',@Resize_clbk);
  
    drawnow();
end