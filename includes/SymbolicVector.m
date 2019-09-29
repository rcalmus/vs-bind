% Ryan Calmus, 2019
% Laboratory of Comparative Neuropsychology
% Newcastle University

classdef SymbolicVector
    properties
        vector
        dimensions
        sparsity = 0.05
        name
        visualisationColor = [65 113 156]/255
        visualisationLineWidth = 7
    end
    properties (Constant, Access = private)
        default_dimensions = 4096;
    end
    properties (Access = private)
        NormaliseVectorsToUnitLength = true
        preferredDirectionVectors = [];
        alphas = [];
        betas = [];
    end
    methods
        
        function out = cleanup(obj,varargin)
            out = obj.associate([varargin{:}],[varargin{:}]);
        end
        
        function out = autoassociate(obj,varargin)
            out = obj.associate([varargin{:}],[varargin{:}]);
        end
        
        % See Crawford (2014)
        function out = associate(obj,keys,values)
            thres = 0.05;
            out = SymbolicVector.zero();
            for k = 1:length(keys)
                simil = dot(keys(k),obj);
                scale = double(simil > thres);
                out = out + scale*values(k);
            end
        end
        
        % Concatenate the content of 2 representations (i.e. concatenate
        % their dimensions), and then (optionally) normalise the result
        % This builds up more complex representations from individual parts
        % where individual parts of the vector are functionally independent
        % from each other.
        function newVSAPopulation = join(obj,obj2)
            newVSAPopulation = SymbolicVector(obj.dimensions+obj2.dimensions,obj.sparsity,[]);
            newVSAPopulation.vector = [obj.vector,obj2.vector];
            newVSAPopulation.normaliseVector();
        end
        
        % Create a new VSA representation from selected dimensions of the
        % current representation and optionally normalise.
        % Part and Join are complementary functions - split up a vector
        % using part, and then concatenate these parts using join.
        % They are carefully named to avoid confusion with the
        % superposition and binding operators.
        function newVSAPopulation = part(obj,indices)
            newVSAPopulation = SymbolicVector(length(indices),obj.sparsity,[],length(indices));
            newVSAPopulation.vector = obj.vector(indices);
            newVSAPopulation.normaliseVector();
        end
        
        % Add two vectors together but only after offsetting one against
        % the other by n dimensions. This allows the encoding space for one
        % vector to be overlapped with a hypothetically different encoding
        % space for another vector. Contrast this with either completely
        % segregated (never combined) vectors, or completely overlapping
        % (normal VSA output) vectors, that can be visualised as being part
        % of a non-overlapping tiling of the cortex.
        function obj3 = offsetSuperposition(obj,obj2,n)
            obj3 = obj;
            if nargin < 3 || isempty(n)
                n = 0;
            end
            padding = repmat(0,1,n);
            obj3.vector = [obj3.vector(:);padding(:)] + [padding(:);obj2.vector(:)];
            length(obj3.vector)
            obj3.dimensions = length(obj3.vector);
        end
        
        function obj = withNormalisationOff(obj)
            for i = 1:length(obj)
                obj(i).NormaliseVectorsToUnitLength = false;
            end
        end
        
        function obj = withNormalisationOn(obj)
            for i = 1:length(obj)
                obj(i).NormaliseVectorsToUnitLength = true;
            end
        end

        % Plot one or more vectors mapped to a 3d coordinate system,
        % projected onto a sphere, and with the option to plot several in
        % one go, with angles between them
        function plot3d(obj,labels)
            try
                if nargin < 2
                    labels = [];
                end
                figure;
                set(gcf,'color','w');
                [x,y,z] = sphere(40);
                hold on
                axis equal

                pxy = patch([1 -1 -1 1], [1 1 -1 -1], [0 0 0 0], [1 1 -1 -1]);
                set(pxy, 'FaceAlpha', 0.1);
                set(pxy, 'LineStyle', 'None');

                pxz = patch([1 -1 -1 1], [0 0 0 0], [1 1 -1 -1], [1 1 -1 -1]);
                set(pxz, 'FaceAlpha', 0.1);
                set(pxz, 'LineStyle', 'None');

                pyz = patch([0 0 0 0],  [-1 1 1 -1], [1 1 -1 -1], [1 1 -1 -1]);
                set(pyz, 'FaceAlpha', 0.1);
                set(pyz, 'LineStyle', 'None');

                if isscalar(obj)
                    error('plot3d makes no sense in an isolated context; provide a vector of SymbolicVectors to begin plotting');
                    return
                end

                for i = 1:length(obj)
                    X(i,:) = obj(i).vector;
                end
                
                %curDir = pwd;
                %statsDir = toolboxdir('stats');
                %cd(fullfile(statsDir,'stats'));
                
                [P,mapp] = pca(X,3); %stats toolbox required %%drtoolbox required
                
                %cd(curDir);
                
                for i = 1:length(obj)
                    p = [P(i,1),P(i,2),P(i,3)];
                    Pt = P*1.2;
                    if ~isempty(labels)
                        t = text(Pt(i,1),Pt(i,2),Pt(i,3),labels{i});
                    else
                        t = text(Pt(i,1),Pt(i,2),Pt(i,3),obj(i).name);
                    end
                    t.FontSize = 30;
                end

                global ColorOrder;
                ColorOrder = zeros(length(obj),3);

                charValues = 'rgbcmywk'.';
                rgbValues = [eye(3); 1-eye(3); 1 1 1; 0 0 0];

                for i = 1:length(obj)
                    if ischar(obj(i).visualisationColor)
                        [isColor,colorIndex] = ismember(obj(i).visualisationColor(:),charValues);
                        assert(all(isColor),'convert_color:badInputContents',...
                               'String input can only contain the characters ''rgbcmywk''.');
                        outColor = rgbValues(colorIndex,:);
                    else
                        outColor = obj(i).visualisationColor;
                    end
                    ColorOrder(i,:) = outColor;
                end

                set(gca,'ColorOrder',ColorOrder);
                w = warning;
                warning('off'); % hide ColorOrder warnings from arrow3
                a = arrow3(repmat([0,0,0],length(obj),1),P,'o5',5,7,1);
                warning(w);
                view(50,20);
                ax = gca;
                ax.Visible = false;
                title(sprintf('3-Dimensional Embeddings of High-Dimensional\nSparse Distributed Representations'));
            catch me
                me
               disp('Something went wrong displaying 3D vector embeddings. Please ensure you have the data reduction toolbox (drtoolbox) installed to undertake PCA.'); 
            end
        end
        
        function obj = set.dimensions(obj,v)
            if v ~= length(obj.vector)
                obj.vector = zeros(1,v);
            end
            obj.dimensions = v;
        end
        
        function result = mul(objs)
            result = objs(1);
            for i = 2:length(objs)
                result = result*objs(i);
            end
        end
        
        function obj2 = cos(obj)
            obj2 = obj;
            obj2.vector = cos(obj.vector);
        end

        function obj = circshift(obj,shiftAmount)
            if nargin < 2 || isempty(shiftAmount)
                shiftAmount = 1;
            end
            obj.vector = circshift(obj.vector,shiftAmount);
        end

        function obj = SymbolicVector(dimensions,sparsity,name,normalise)
            global allSymbolicVectors;
            
            if isempty(allSymbolicVectors)
                allSymbolicVectors = {};
            end
            
            if nargin < 3|| isempty(name)
                obj.name = 'Unnamed VSA population';
            else
                obj.name = name;
            end
            
            if nargin < 1 || isempty(dimensions)
                obj.dimensions = SymbolicVector.default_dimensions; 
            else
                obj.dimensions = dimensions;
            end
            
            if nargin >= 2 && ~isempty(sparsity)
                if isnumeric(sparsity) % if numeric, it's a proportional figure
                    obj.sparsity = sparsity;
                elseif isstring(sparsity) % if a string, let us say it's an absolute number of neurons being specified
                    obj.sparsity = str2double(sparsity)/obj.dimensions;
                end
            end
            
            if nargin >= 4 && ~isempty(normalise)
                obj.NormaliseVectorsToUnitLength = normalise;
            end
            
            obj.vector =  full(sprandn(1,obj.dimensions,obj.sparsity));
            
            obj = obj.normaliseVector();
        end

        % Like spy, but with multiple marker sizes; an intermediate for
        % graphics less complex than a full imagesc, but more detailed than
        % "spy"
        function h = multispy(obj,plotImagesc,littleCircles,varargin)

            if nargin < 2
                plotImagesc = true;
            end
            if nargin < 3
                littleCircles = true;
            end
            
            symbolArg = 0;
            symbolDisplayText = [];

            if nargin > 3
                try
                    h = subplot(varargin{1:3});
                    symbolArg = 4;
                catch
                    h = gca;
                    symbolArg = 1;
                end

                if length(varargin) >= symbolArg
                    symbolDisplayText = varargin{symbolArg};
                end
            else
                h = gca;
            end

            squareSize = floor(sqrt(obj.dimensions));
            displayVector = obj.vector(1:squareSize^2);
            displayMat = reshape(displayVector,squareSize,squareSize); % we now have a standard matrix
            displayMat = abs(round(displayMat,4));

            if plotImagesc
                if littleCircles
                    circimagesc(displayMat);
                else
                    imagesc(displayMat);
                end
            else
                levels = 10;

                hold off

                units = get(gca,'units');
                set(gca,'units','pixels');
                pos = get(gca,'position');

                markersize = max(4,min(14,round(6*min(pos(3:4))/max(squareSize+1,squareSize+1))));
                set(gca,'units','Normalized');
                markersize = 0.5/max(squareSize+1,squareSize+1);

                co = get(gca,'colororder');
                color = co(1,:);

                for i = 1:levels
                    binMat = displayMat>prctile(displayMat(displayMat~=0),(100/levels)*i);
                    [row, col] = find(binMat);
                    hold on;
                    scatter(col,row,markersize*2*(i/levels),'b','MarkerEdgeColor',color,'MarkerFaceColor',color);
                    set(gca,'ydir','Reverse')
                    ylim([0 squareSize+1]);
                    xlim([0 squareSize+1]);
                end
            end

            scaleFactors = pbaspect;

            if isempty(obj.visualisationColor)
                obj.visualisationColor = [65 113 156]/255;
            end

            xMax = xlim;
            yMax = ylim;
            inset = obj.visualisationLineWidth*0.024;

            rectangle('Position',[0+inset+0.05 0+inset+0.05 xMax(2)+inset-0.1 yMax(2)+inset-0.1],'LineWidth',obj.visualisationLineWidth,'FaceColor','none','EdgeColor',obj.visualisationColor);
            

            if isempty(symbolDisplayText)
                textInterpreter = 'tex';
                try
                    if strcmpi(obj.name,'Unnamed VSA population') || isempty(obj.name)
                        titleText = inputname(1);
                    else
                        titleText = obj.name;
                    end
                catch
                    titleText = obj.name;
                end

                titleText = strrep(titleText,'\*','\{asterisk}');
                titleText = strrep(titleText,'*',sprintf(' \x2297'));
                titleText = strrep(titleText,'\{asterisk}','*');

                tH = title(titleText,'interpreter',textInterpreter,'Color',obj.visualisationColor);
            end

            if ~isempty(symbolDisplayText)
                symbolDisplayText = strrep(symbolDisplayText,'\*','\{asterisk}');
                symbolDisplayText = strrep(symbolDisplayText,'*',sprintf('\x2297'));
                symbolDisplayText = strrep(symbolDisplayText,'\{asterisk}','*');

                symbolDisplayText = strcat('\bf{',[' ',symbolDisplayText],'}'); % make bold and nudge right with a space as otherwise, when copied from Matlab, it looks too far left

                outPos = get(gca,'OuterPosition');

                
                boxPos = get(gca,'InnerPosition');
                fontSize = 18;
                
                titleHeight = abs(outPos(4)-boxPos(4));

                scaleRelativeToFont = 0.003;
                boxWidth = 0.06*(scaleFactors(1)/scaleFactors(2));

                boxPos(2) = boxPos(2) + boxPos(4) - boxWidth*(scaleFactors(2)/scaleFactors(1));
                boxPos(1) = boxPos(1) + boxPos(3) - boxWidth/2;
                boxPos(3) = boxWidth;
                
                boxPos(4) = boxPos(3)*(scaleFactors(2)/scaleFactors(1));
           
                tY = 1;
                if strcmpi(get(gca,'ydir'),'Reverse')
                    tY= 0;
                end

                boxToFontRatio = fontSize/10;
                set(gca,'Units','Normalized');
                r = rectangle(gca,'Position',[0-boxToFontRatio tY-boxToFontRatio, 2*boxToFontRatio, 2*boxToFontRatio]);
                r.LineWidth = obj.visualisationLineWidth;
                r.FaceColor = 'white';
                r.EdgeColor = obj.visualisationColor;

                tH = text(0,tY,symbolDisplayText);
                tH.HorizontalAlignment = 'center';
                tH.VerticalAlignment = 'middle';
                tH.FontSize = fontSize;
                try
                    tH.FontName = 'Cambria Math';
                catch
                end

            end

            set(gca,'FontSize',20);
            set(gca,'xtick',[]);
            set(gca,'xticklabel',[]);
            set(gca,'ytick',[]);
            set(gca,'yticklabel',[]);
            set(gca,'xlabel',[]);
            set(gca,'ylabel',[]);
            set(gca,'LineWidth',obj.visualisationLineWidth);
            set(gca,'Box','off');
            set(gca,'Clipping','off');

            set(gca,'XColor',obj.visualisationColor);
            set(gca,'YColor',obj.visualisationColor);
            set(gca,'ZColor',obj.visualisationColor);

            boundPos = get(gca,'InnerPosition');

            capExtraArgs = {'LineWidth',1,'LineStyle','none','FaceColor',obj.visualisationColor};
            % Correct a Matlab "bug" with visualisation of big boxes around the
            % plot, where the top left corner is not cleanly rendered:
            capRadius = 0.0007*obj.visualisationLineWidth; % prettify the elbow of the arrow
            axesHandles = findobj(get(gcf,'Children'), 'flat','Type','axes');
            axis('off');
            colormap(flipud(colormap('gray')));
        end

        function spy(obj,varargin)
             squareSize = floor(sqrt(obj.dimensions));
            displayVector = obj.vector(1:squareSize^2);
            displayMat = reshape(displayVector,squareSize,squareSize); % we now have a standard matrix

            binarisedMat = abs(round(displayMat,5))>0;
            spy(binarisedMat,varargin{:});

            try
                if strcmpi(obj.name,'Unnamed VSA population') || isempty(obj.name)
                    title(inputname(1));
                else
                    title(obj.name)
                end
            catch
                title(obj.name);
            end

            set(gca,'xtick',[]);
            set(gca,'xticklabel',[]);
            set(gca,'ytick',[]);
            set(gca,'yticklabel',[]);
            set(gca,'xlabel',[]);
            set(gca,'ylabel',[]);
        end
        
        function result = isSquare(obj)
            result = (round(sqrt(obj.dimensions))==sqrt(obj.dimensions));
        end
        
        function result = add(obj,anotherPopulation)
            if ~isscalar(obj)
                if size(obj) == size(anotherPopulation)
                    for i = 1:length(obj)
                        result(i) = add(obj(i),anotherPopulation(i));
                    end    
                else
                    for i = 1:length(obj)
                        result(i) = add(obj(i),anotherPopulation);
                    end
                end
                return
            end
            if ~isscalar(anotherPopulation)
                for i = 1:length(anotherPopulation)
                    result(i) = add(obj,anotherPopulation(i));
                end
                return
            end
            otherVector = anotherPopulation.vector;
            if obj.dimensions ~= anotherPopulation.dimensions
                otherVector = resample(anotherPopulation.vector,obj.dimensions,anotherPopulation.dimensions);
            end
            result = SymbolicVector(obj.dimensions);
            
            result.NormaliseVectorsToUnitLength = obj.NormaliseVectorsToUnitLength;
            
            result.vector = obj.vector + otherVector;
            
            result = result.normaliseVector();
        end
        
        function result = subtract(obj,anotherPopulation)
            
            if ~isscalar(obj)
                if size(obj) == size(anotherPopulation)
                    for i = 1:length(obj)
                        result(i) = subtract(obj(i),anotherPopulation(i));
                    end    
                else
                    for i = 1:length(obj)
                        result(i) = subtract(obj(i),anotherPopulation);
                    end
                end
                return
            end
            if ~isscalar(anotherPopulation)
                for i = 1:length(anotherPopulation)
                    result(i) = subtract(obj,anotherPopulation(i));
                end
                return
            end
            
            otherVector = anotherPopulation.vector;
            if obj.dimensions ~= anotherPopulation.dimensions
                otherVector = resample(anotherPopulation.vector,obj.dimensions,anotherPopulation.dimensions);
            end
            result = SymbolicVector(obj.dimensions);
            
            result.NormaliseVectorsToUnitLength = obj.NormaliseVectorsToUnitLength;
            
            result.vector = obj.vector - otherVector;
            
            result = result.normaliseVector();
        end
        
        % Threshold the representation so elements above the threshold
        % remain unchanged, and the others are zeroed.
        function obj2 = gt(obj,thres)
            obj2 = obj;
            obj2.vector = (obj.vector>thres) .* obj.vector;
        end
        
        function obj2 = lt(obj,thres)
            obj2 = obj;
            obj2.vector = (obj.vector<thres) .* obj.vector;
        end
        
        function obj2 = lte(obj,thres)
            obj2 = obj;
            obj2.vector = (obj.vector<=thres) .* obj.vector;
        end
        
        function obj2 = gte(obj,thres)
            obj2 = obj;
            obj2.vector = (obj.vector>=thres) .* obj.vector;
        end
        
        % More like "trinarise" the representation into trits (-1,0,1).
        function obj2 = binarise(obj)
            obj2 = obj;
            obj2.vector = (-1 .* (obj.vector<0)) + (1 .* (obj.vector>0));
            obj2.normaliseVector();
        end
        
        % Take absolute of VSA vector
        function obj2 = abs(obj)
            obj2 = obj;
            obj2.vector = abs(obj.vector);
        end
        
        % Calculate the difference between each vector symbol of a list
        % of separate vectors
        function result = diff(inputs)
            result = repmat(SymbolicVector.zero(inputs(1).dimensions),1,length(inputs)-1);
            for i = 1:(length(inputs)-1)
                result(i) = inputs(i+1)-inputs(i);
            end
        end
        
        function result = dot(obj,obj2)
            result = dot(obj.vector,obj2.vector);
        end
        
        function result = norm(obj)
            result = norm(obj.vector);
        end
        
        % sort the features of obj.vector in ascending order
        function [obj2,i] = sort(obj)
            obj2 = obj;
            [obj2.vector,i] = sort(obj2.vector);
        end
        
        % Tag each item of a list of separate vectors, each one tagged
        % by binding with a power of positionTag; return as a mean result
        function result = tag(inputs,positionTag)
            result = inputs;
            for i = 1:(length(inputs))
                result(i) = inputs(i)*(positionTag^i);
            end
            result = mean(result);
        end
        
        function obj3 = prod(obj,obj2)
            obj3 = and(obj,obj2);
        end
        
        function obj3 = div(obj,obj2)
            obj3 = obj;
            obj3.vector = obj.vector ./ obj2.vector;
            obj3.normaliseVector();
        end
        
        function obj3 = and(obj,obj2)
            obj3 = obj;
            obj3.vector = obj.vector .* obj2.vector;
            obj3.normaliseVector();
        end
        
        % Reduce dimensionality of representation by adding random
        % permutations of the representation together and only taking first
        % N elements of the result, where N = round(obj.dimension *
        % proportion)
        function obj2 = reduce(obj,proportion)
            newSize = round(proportion*obj.dimensions);
            obj2 = SymbolicVector.zero(newSize);
            for i = 0:round(1/proportion)-1
                obj2 = obj2 + part(permute(obj,i),1:newSize);
            end
        end
        
        % Permute the vector in a consistent way according to a given seed,
        % permutationSeed. All positive values are unique permutations. All
        % negative values are the inverse permutation of their absolute
        % positive. Vectors of seeds can be given, in which case each
        % permutation is performed serially. As a consequence, providing a
        % seed of [s,-s] will always result in the original vector.
        function obj2 = permute(obj,permutationSeed)
            obj2 = obj;
            if nargin < 2 || isempty(permutationSeed)
                permutationSeed = 1;
            end
            if ~isscalar(permutationSeed)
                for i = 1:length(permutationSeed)
                    obj2 = permute(obj2,permutationSeed(i));
                end
            else
                if permutationSeed >= 0
                    seed = permutationSeed;
                    invert = false;
                else
                    seed = -permutationSeed;
                    invert = true;
                end
                if seed ~= 0
                    s = rng;
                    rng(seed);

                    if ~invert
                        obj2.vector = obj2.vector(randperm(obj2.dimensions));
                    else
                        obj2.vector(randperm(obj2.dimensions)) = obj2.vector;
                    end

                    rng(s);
                end
            end
        end
        
        function result = superpose(obj,anotherPopulation)
            result = obj.add(anotherPopulation);
        end
        
        function result = bundle(obj,anotherPopulation)
            result = obj.add(anotherPopulation);
        end
        
        function result = bind(obj,anotherPopulation)
            otherVector = anotherPopulation.vector;
            if obj.dimensions ~= anotherPopulation.dimensions
                otherVector = resample(anotherPopulation.vector,obj.dimensions,anotherPopulation.dimensions);
            end
            result = SymbolicVector(obj.dimensions);
            result.NormaliseVectorsToUnitLength = obj.NormaliseVectorsToUnitLength;
            result.vector = cconv(obj.vector,otherVector,obj.dimensions);

            result = result.normaliseVector();
        end
        
        function result = convolve(obj,anotherPopulation)
            result = obj.bind(obj,anotherPopulation);
        end
        
        % Proof this works:
        % selfBindResult = @(x) (x*~x);
        % (mean(arrayfun(@(~) selfBindResult(SymbolicVector),1:2000)) == SymbolicVector.identity)
        % = 0.9998
        function result = involute(obj) 
            
            if ~isscalar(obj)
                for i = 1:length(obj)
                    result(i) = involute(obj(i),anotherPopulation);
                end
                return
            end
            
             result = SymbolicVector(obj.dimensions);
             result.NormaliseVectorsToUnitLength = obj.NormaliseVectorsToUnitLength;

             result.vector = zeros(1,obj.dimensions);
             for i = 0:length(result.vector)-1
                result.vector(i+1) = obj.vector(mod(-i,obj.dimensions)+1);
             end
             
            result = result.normaliseVector();
        end
        
        function result = invert(obj)
            result = obj.involute();
        end
        
        % Misleading name. Correlate or use dot product to find similarity.
        function r = correlate(obj,anotherPopulation)
            if ~isscalar(obj)
                for i = 1:length(obj)
                    r(i) = correlate(obj(i),anotherPopulation);
                end
                return
            end
            if ~isscalar(anotherPopulation)
                for i = 1:length(anotherPopulation)
                    r(i) = correlate(obj,anotherPopulation(i));
                end
                return
            end
            otherVector = anotherPopulation.vector;
            if obj.dimensions ~= anotherPopulation.dimensions
                otherVector = resample(anotherPopulation.vector,obj.dimensions,anotherPopulation.dimensions);
            end
            
            useCorrCoef = true;
            if useCorrCoef
                Rmatrix = corrcoef(obj.vector,otherVector);
                r = Rmatrix(2,1);
            else
                r = dot(obj.vector,otherVector);
            end
        end
        
        function r = similarity(obj,anotherPopulation)
            r = obj.correlate(anotherPopulation);
        end
        
        function r = compare(obj,anotherPopulation)
            r = obj.correlate(anotherPopulation);
        end
        
%         function [cleaned maxRval] = cleanup(obj,varargin)
%             comparisonResults = [];
%             for i = 1:length(varargin)
%                 comparisonResults(end+1) = obj.correlate(varargin{i});
%             end
%             [maxRval,i] = max(comparisonResults);
%             cleaned = SymbolicVector(obj.dimensions,obj.sparsity,obj.name);
%             cleaned.NormaliseVectorsToUnitLength = obj.NormaliseVectorsToUnitLength;
%             
%             otherVector = varargin{i}.vector;
%             
%             if obj.dimensions ~= varargin{i}.dimensions
%                 otherVector = resample(varargin{i}.vector,obj.dimensions,varargin{i}.dimensions);
%             end
%             cleaned.vector = otherVector;
%             
%             cleaned = cleaned.normaliseVector();
%         end
        
        function r = compareAll(obj,varargin)
            if length(obj) > 1
                Z = zeros(length(obj),length(varargin));
                for iObj = 1:length(obj)
                    Z(iObj,:) = obj(iObj).compareAll(varargin{:});
                end

                if nargout == 0
                figure;
                bar3(Z);
                objName = string(unique(categorical({obj.name})))
                title(sprintf('Similarity of ''%s'' to specified representations over time',objName(1)));
                zlabel('Similarity');
                ylabel('t (steps)');

                populationLabels = {};
                for i = 1:length(varargin)
                    try
                        if ~isempty(inputname(1+i))
                            populationLabels{end+1} = inputname(1+i);
                        else
                            populationLabels{end+1} = varargin{i}.name;
                        end
                    catch
                        populationLabels{end+1} = varargin{i}.name;
                    end
                end

                xticks(1:length(populationLabels));
                xticklabels(populationLabels);
                xtickangle(45);
                end
                
                r = Z;
                return;
            end
            comparisonResults = [];
            populationLabels = {};
            for i = 1:length(varargin)
                comparisonResults(end+1) = obj.correlate(varargin{i});
                try
                    if ~isempty(inputname(1+i))
                        populationLabels{end+1} = inputname(1+i);
                    else
                        populationLabels{end+1} = varargin{i}.name;
                    end
                catch
                    populationLabels{end+1} = varargin{i}.name;
                end
            end
            if nargout == 0
                figure;
                bar(comparisonResults);
                ylim([-1 1]);
                xticks(1:length(comparisonResults));
                xticklabels(populationLabels);
                xtickangle(45);
                try
                    objName = inputname(1);
                catch
                    objName = obj.name;
                end
                title(sprintf('Similarity of ''%s'' to specified representations',objName));
            end
            r = comparisonResults;
        end
        
        function h = plot(obj)
            if ~isscalar(obj)
                for i = 1:length(obj)
                    obj(i).plot;
                end
                return;
            end
            if nargout == 0
                figure;
            end
            h = gca;
            squareSize = floor(sqrt(obj.dimensions));
            
            if ~obj.isSquare()
                fprintf('SymbolicVector - Warning: choosing a non-square population size (%d) necessitates approximate visualisation (%d x %d).\n',obj.dimensions,squareSize,squareSize);
            end
            
            displayVector = obj.vector(1:squareSize^2);
            imagesc(reshape(displayVector,squareSize,squareSize));
            colorbar;
            title(obj.name);
        end
        
        function obj = normaliseVector(obj)
            if obj.NormaliseVectorsToUnitLength
                
                mag = norm(obj.vector);
                obj.vector = bsxfun(@rdivide, obj.vector, mag);
            end
            obj.vector = real(obj.vector);
            obj.vector(isnan(obj.vector)) = 0;
        end
        
        function printSummary(obj,title)
            fprintf('%s: Max = %f; Min = %f, Mean = %f\n',title,max(obj.vector),min(obj.vector),mean(obj.vector));
        end
        
        %% operator overloading:
        
        function r = plus(obj,obj2)
           r = obj.add(obj2);
        end
        function r = minus(obj,obj2)
           r = obj.subtract(obj2);
        end
        function r = uminus(obj)
           z = SymbolicVector.zero();
           r = z.subtract(obj);
        end

        function r = power(obj,exponent)
            r = obj.mpower(exponent);
        end
        
        function r = mpower(obj,exponent)
            % see Plate (1995) for an explanation:
            r = obj;
            
            r.vector = ifft(fft(obj.vector).^exponent);
            
            r = r.normaliseVector();
        end

        function r = mtimes(firstTerm,secondTerm)
            if isnumeric(secondTerm)
                if isscalar(secondTerm)
                    r = SymbolicVector(firstTerm.dimensions);
                    r.NormaliseVectorsToUnitLength = firstTerm.NormaliseVectorsToUnitLength;
                    
                    r.vector = firstTerm.vector .* secondTerm;
                else
                    if size(firstTerm) == size(secondTerm)
                        for i = 1:length(secondTerm)
                            r(i) = SymbolicVector(firstTerm(i).dimensions);
                            
                            r(i).NormaliseVectorsToUnitLength = firstTerm(i).NormaliseVectorsToUnitLength;
                            
                            r(i).vector = firstTerm(i).vector .* secondTerm(i);
                        end    
                    else
                        for i = 1:length(secondTerm)
                            r(i) = SymbolicVector(firstTerm.dimensions);
                            
                            r(i).NormaliseVectorsToUnitLength = firstTerm.NormaliseVectorsToUnitLength;
                            
                            r(i).vector = firstTerm.vector .* secondTerm(i);
                        end
                    end
                end
            elseif isnumeric(firstTerm)
                if isscalar(firstTerm)
                    r = SymbolicVector(secondTerm.dimensions);
                    r.NormaliseVectorsToUnitLength = secondTerm.NormaliseVectorsToUnitLength;
                    
                    r.vector = secondTerm.vector .* firstTerm;
                else
                    if size(firstTerm) == size(secondTerm)
                        for i = 1:length(firstTerm)
                            r(i) = SymbolicVector(secondTerm(i).dimensions);
                            r(i).NormaliseVectorsToUnitLength = secondTerm(i).NormaliseVectorsToUnitLength;
                            
                            r(i).vector = secondTerm(i).vector .* firstTerm(i);
                        end    
                    else
                        for i = 1:length(firstTerm)
                            r(i) = SymbolicVector(secondTerm.dimensions);
                            
                            r(i).NormaliseVectorsToUnitLength = secondTerm.NormaliseVectorsToUnitLength;
                            
                            r(i).vector = secondTerm.vector .* firstTerm(i);
                        end
                    end
                end
            else
                if ~isscalar(firstTerm)
                    if size(firstTerm) == size(secondTerm)
                        for i = 1:length(firstTerm)
                            r(i) = firstTerm(i).bind(secondTerm(i));
                        end
                    else
                        for i = 1:length(firstTerm)
                            r(i) = firstTerm(i).bind(secondTerm);
                        end
                    end
                    return
                end
                if ~isscalar(secondTerm)
                    for i = 1:length(secondTerm)
                        r(i) = firstTerm.bind(secondTerm(i));
                    end
                    return
                end
                r = firstTerm.bind(secondTerm);
            end
        end
        
        function r = not(obj)
           r = obj.involute();
        end
        function r = times(obj,secondTerm)
            r = mtimes(obj,secondTerm);
        end
        function r = eq(obj,obj2)
            r = obj.correlate(obj2);
        end
        
        function mu = mean(obj,varargin)
            allVectors = zeros(length(obj),obj(1).dimensions);

            for iObj = 1:length(obj)
                if obj(1).dimensions ~= obj(iObj).dimensions
                    allVectors(iObj,:) = resample(obj(iObj).vector,obj(1).dimensions,obj(iObj).dimensions);
                else
                    allVectors(iObj,:) = obj(iObj).vector;
                end
            end
            
            muVector = mean(allVectors,1,varargin{:});
            
            mu = SymbolicVector(obj(1).dimensions,obj(1).sparsity,[],obj(1).dimensions);
            mu.NormaliseVectorsToUnitLength = obj.NormaliseVectorsToUnitLength;
            
            mu.vector = muVector;
            mu = mu.normaliseVector();
        end
    end
    methods (Static)
        
        % Convolutional identity vector (i.e. the approximate result of A * ~A)
        function aPopulation = identity(dimensions,name)
            if nargin < 2
                name = [];
            end
            if nargin < 1
                dimensions = [];
            end
            aPopulation = SymbolicVector(dimensions,[],name);
            aPopulation.vector = zeros(1,aPopulation.dimensions);
            aPopulation.vector(1) = 1;
        end
        
        function aPopulation = zero(dimensions,name)
            if nargin < 2
                name = [];
            end
            if nargin < 1
                dimensions = [];
            end
            aPopulation = SymbolicVector(dimensions,[],name);
            aPopulation.vector = zeros(1,aPopulation.dimensions);
        end
    end
    methods (Static, Access = private)
        function result = endsWithExt(aString,extension)
            result = any((strfind(aString,extension)+(length(extension)-1))==length(aString));
        end

        function h = circle2(x,y,r,varargin)
            d = r*2;
            px = x-r;
            py = y-r;
            h = annotation('Ellipse',[px py d d],varargin{:});
        end
        function h = centredRect(x,y,r,varargin)
            d = r*2;
            px = x-r;
            py = y-r;
            h = annotation('Rectangle',[px py d d],varargin{:});
        end

        function [datX,datY] = figureSpaceToDataSpaceCoords(hAx,x,y)
            axun = get(hAx,'Units');
            set(hAx,'Units','normalized');
            axpos = get(hAx,'Position');
            axlim = axis(hAx);
            axwidth = diff(axlim(1:2));
            axheight = diff(axlim(3:4));

            %% Transform data

            datX = ((x - axpos(1))*axwidth/axpos(3))+axlim(1);
            datY = ((y - axpos(2))*axheight/axpos(4))+axlim(3);

            if strcmpi(get(hAx,'ydir'),'Reverse')
                datY = (axheight-datY);
            end

            %% Restore axes units
            set(hAx,'Units',axun);
        end
    end
end