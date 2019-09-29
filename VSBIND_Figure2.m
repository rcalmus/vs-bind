% Ryan Calmus, 2019
% Laboratory of Comparative Neuropsychology
% Newcastle University

function VSBIND_Figure2(seed)
    
    disp('Figure 2 diagrams (2 windows) rendering. Wait for the figure to render, and then resize the figure window slowly until the aspect ratio is correct.');

    vectorSize = 256;
    usecircles = true;
    lineWidth = 7;
    
    randState = rng;
    sparsity = 0.05;
    if nargin < 1  % choose the prettiest vectors of certain sizes to plot by setting the random number generator seed
        if vectorSize == 64
            seed = 57;
        elseif vectorSize == 1024
            seed = 15;
        elseif vectorSize == 256
            seed = 5;
            sparsity = 0.1;
        end
    end

    rng(seed);

    A = SymbolicVector(vectorSize,sparsity);
    A.name = '\bf{A}';
    A.visualisationLineWidth = lineWidth;

    B = SymbolicVector(vectorSize,sparsity);
    B.name = '\bf{B}';
    B.visualisationLineWidth = lineWidth;
    
    aSum = A+B;
    aSum.name = '\color{black}A+B';
    aSum.visualisationColor = 'k';
    aSum.visualisationLineWidth = lineWidth;

    aBinding = A*B;
    aBinding.name = '\color{black}A*B';
    aBinding.name = strrep(aBinding.name,'*',sprintf(' \x2297')); % to show actual circular convolution symbol
    aBinding.visualisationColor = 'k';
    aBinding.visualisationLineWidth = lineWidth;

    plot3d([A,B,aSum,aBinding])
    
    anInverse = ~B;
    anInverse.name = '\color{black}¬B';
    anInverse.visualisationColor = 'k';
    anInverse.visualisationLineWidth = lineWidth;

    reconstructA = aBinding*anInverse;
    reconstructA.name = '\color{black}(A*B)*¬B \approx A';
    reconstructA.name = strrep(reconstructA.name,'*',sprintf(' \x2297')); % to show actual circular convolution symbol
    reconstructA.visualisationColor = 'k';
    reconstructA.visualisationLineWidth = lineWidth;

    h = scalingFigure();

    APlot = subplot(4,5,2);
    multispy(A,true,usecircles);

    BPlot = subplot(4,5,12);
    multispy(B,true,usecircles);
    
    SumPlot = subplot(4,5,6);
    multispy(aSum,true,usecircles);

    BindingPlot = subplot(4,5,8);
    multispy(aBinding,true,usecircles);

    InversePlot = subplot(4,5,18);
    multispy(anInverse,true,usecircles);

    ReconstructPlot = subplot(4,5,14);
    multispy(reconstructA,true,usecircles);
    
    CleanupPlot = subplot(4,5,20);
    multispy(A,true,usecircles);

    arrowBetweenSubplots(APlot,BindingPlot);
    arrowBetweenSubplots(BPlot,BindingPlot);
    arrowBetweenSubplots(APlot,SumPlot);
    arrowBetweenSubplots(BPlot,SumPlot);
    arrowBetweenSubplots(BindingPlot,ReconstructPlot);
    arrowBetweenSubplots(InversePlot,ReconstructPlot);
    arrowBetweenSubplots(BPlot,InversePlot,'yx');
    arrowBetweenSubplots(ReconstructPlot,CleanupPlot,[],'LineStyle',':');
    rng(randState);
end

% Ryan Calmus, Petkov Lab, 2019 - Draw diagonal or elbow arrows to connect subplots
% sensibly.
function arrowBetweenSubplots(startPlotHandle,finishPlotHandle,elbowType,varargin)
    startX = startPlotHandle.Position(1)+(startPlotHandle.Position(3))/2;
    finishX = finishPlotHandle.Position(1)+(finishPlotHandle.Position(3))/2;
    startY = startPlotHandle.Position(2)+(startPlotHandle.Position(4))/2;
    finishY = finishPlotHandle.Position(2)+(finishPlotHandle.Position(4))/2;

    if ~isempty(get(get(startPlotHandle,'Title'),'String'))
        startTitleHeight = 0.02;
    else
        startTitleHeight = -0.03;
    end

    if ~isempty(get(get(finishPlotHandle,'Title'),'String'))
        finishTitleHeight = 0.02;
    else
        finishTitleHeight = -0.03;
    end

    if nargin < 3 || isempty(elbowType)
        elbowType = 'xy';
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

    headsize = (10/3)*startPlotHandle.LineWidth;
    extraArgsLine = [{'LineWidth',startPlotHandle.LineWidth,'Color',startPlotHandle.XColor},varargin];
    extraArgsArrow = [{'LineWidth',startPlotHandle.LineWidth,'HeadWidth',headsize,'HeadLength',headsize,'Color',startPlotHandle.XColor},varargin];
    capExtraArgs = [{'LineWidth',0.1,'LineStyle','-','FaceColor',startPlotHandle.XColor,'EdgeColor',startPlotHandle.XColor},varargin];

    if strcmpi(elbowType,'diagonal')
        annotation('arrow', [startX finishX],[startY finishY],extraArgsArrow{:});
    else
        capRadius = 0.0006*startPlotHandle.LineWidth;
        if finishX-startX ~= 0
            if finishY-startY ~= 0
                if strcmpi(elbowType,'xy')
                    annotation('line', [startX finishX],[startY startY],extraArgsLine{:});
                    circle2(finishX,startY,capRadius,capExtraArgs{:});  % to add corner
                    annotation('arrow', [finishX finishX],[startY finishY],extraArgsArrow{:});
                else
                    annotation('line', [startX startX],[startY finishY],extraArgsLine{:});
                    circle2(startX,finishY,capRadius,capExtraArgs{:});  % to add corner
                    annotation('arrow', [startX finishX],[finishY finishY],extraArgsArrow{:});
                end
            else
                if startX>finishX
                    startX = startX-(startPlotHandle.Position(3)/2);
                    annotation('arrow', [startX finishX],[startY finishY],extraArgsArrow{:});
                elseif startX<finishX
                    startX = startX+(startPlotHandle.Position(3)/2);
                    annotation('arrow', [startX finishX],[startY finishY],extraArgsArrow{:});
                else
                end
            end
        else
            if startY>finishY
                startY = startY-(startPlotHandle.Position(4)/2);
                annotation('arrow', [startX finishX],[startY finishY],extraArgsArrow{:});
            elseif startY<finishY
                finishY = finishY-(finishPlotHandle.Position(4)/2);
                annotation('arrow', [startX finishX],[startY finishY],extraArgsArrow{:});
            else
            end
        end

    end
end

function h = circle2(x,y,r,varargin)
    d = r*2;
    px = x-r;
    py = y-r;

    scaleFactors = pbaspect;

    h = annotation('Ellipse',[px py-0.001 d d*(scaleFactors(2)/scaleFactors(1))],varargin{:});
end

% Ryan Calmus, Petkov Lab, 2019 - Create dynamically rendering figure
% capable of scaling diagrammatic components with window size.

function h = scalingFigure()
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
            graphicHandles = findall(hObject,'-property', 'LineWidth');
            scale = hObject.Position(3)/hObject.origPos(3);
            scaleY = hObject.Position(4)/hObject.origPos(4);

            arrayfun(@(x) set(x,'LineWidth',get(x,'LineWidth')*scale),graphicHandles);

            graphicHandles = findall(hObject,'-property', 'HeadWidth');
            arrayfun(@(x) set(x,'HeadWidth',get(x,'HeadWidth')*scale),graphicHandles);        
            arrayfun(@(x) set(x,'HeadLength',get(x,'HeadLength')*scale),graphicHandles);

            graphicHandles = findall(hObject,'-property', 'FontSize');
            arrayfun(@(x) set(x,'FontSize',get(x,'FontSize')*scale),graphicHandles);

            hObject.origPos = hObject.Position;
            hObject.updatingSizes = false;
        catch
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
    newPos = h.InnerPosition;
    newPos(3) = newPos(4);
    h.InnerPosition = newPos;
    set(h,'Units','Normalized');

    drawnow();

    h.origPos = h.Position;
    addprop(h,'updatingSizes');
    h.updatingSizes = false;

    set(h,'SizeChangedFcn',@Resize_clbk);
end