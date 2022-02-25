classdef HW4 < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure    matlab.ui.Figure
        PCAButton   matlab.ui.control.StateButton
        LoadButton  matlab.ui.control.Button
        UIAxes_6    matlab.ui.control.UIAxes
        UIAxes_5    matlab.ui.control.UIAxes
        UIAxes_4    matlab.ui.control.UIAxes
        UIAxes_3    matlab.ui.control.UIAxes
        UIAxes_2    matlab.ui.control.UIAxes
        UIAxes      matlab.ui.control.UIAxes
    end

    
    properties (Access = private)
        imgs = {};
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: LoadButton
        function LoadButtonPushed(app, event)
            % Reset PCA button
            app.PCAButton.Value = false;
            
            % Load and show each band image on axes
            app.imgs = {};
            fnames = {'Band1.tiff', 'Band2.tiff', 'Band3.tiff', 'Band4.tiff', 'Band5.tiff', 'Band6.tiff'};
            uiaxes = {app.UIAxes, app.UIAxes_2, app.UIAxes_3, app.UIAxes_4, app.UIAxes_5, app.UIAxes_6};
            for i = 1:length(uiaxes)
                img = imread(fnames{i});
                app.imgs{end+1} = img;
                imshow(img, 'parent', uiaxes{i});
            end
        end

        % Value changed function: PCAButton
        function PCAButtonValueChanged(app, event)
            value = app.PCAButton.Value;
            if value
                % Compute covariance matrix and eigenvectors
                X = [app.imgs{1}(:), app.imgs{2}(:), app.imgs{3}(:), ...
                    app.imgs{4}(:), app.imgs{5}(:), app.imgs{6}(:)]';
                X = double(X);
                m = mean(X,2);
                Cov = cov(X');
                [V, D] = eig(Cov);
                
                % Rearrange eigenvectors with descending order
                [d, index] = sort(diag(D), 'descend');
                V = V(:, index);
                
                % Transform
                y = V'*(X-m);

                % Show PCA images on axes
                img_size = size(app.imgs{1});
                uiaxes = {app.UIAxes, app.UIAxes_2, app.UIAxes_3, app.UIAxes_4, app.UIAxes_5, app.UIAxes_6};
                for i = 1:length(uiaxes)
                    img = reshape(y(i, :), img_size);
                    imshow(img, 'parent', uiaxes{i}, 'DisplayRange', []);
                end
            else
                uiaxes = {app.UIAxes, app.UIAxes_2, app.UIAxes_3, app.UIAxes_4, app.UIAxes_5, app.UIAxes_6};
                for i = 1:length(uiaxes)
                    imshow(app.imgs{i}, 'parent', uiaxes{i});
                end
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1000 600];
            app.UIFigure.Name = 'UI Figure';

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            app.UIAxes.Position = [167 311 254 235];

            % Create UIAxes_2
            app.UIAxes_2 = uiaxes(app.UIFigure);
            app.UIAxes_2.Position = [420 311 254 235];

            % Create UIAxes_3
            app.UIAxes_3 = uiaxes(app.UIFigure);
            app.UIAxes_3.FontWeight = 'bold';
            app.UIAxes_3.Position = [673 311 254 235];

            % Create UIAxes_4
            app.UIAxes_4 = uiaxes(app.UIFigure);
            app.UIAxes_4.Position = [167 46 254 235];

            % Create UIAxes_5
            app.UIAxes_5 = uiaxes(app.UIFigure);
            app.UIAxes_5.Position = [420 46 254 235];

            % Create UIAxes_6
            app.UIAxes_6 = uiaxes(app.UIFigure);
            app.UIAxes_6.FontWeight = 'bold';
            app.UIAxes_6.Position = [673 46 254 235];

            % Create LoadButton
            app.LoadButton = uibutton(app.UIFigure, 'push');
            app.LoadButton.ButtonPushedFcn = createCallbackFcn(app, @LoadButtonPushed, true);
            app.LoadButton.Position = [51 311 79 44];
            app.LoadButton.Text = 'Load';

            % Create PCAButton
            app.PCAButton = uibutton(app.UIFigure, 'state');
            app.PCAButton.ValueChangedFcn = createCallbackFcn(app, @PCAButtonValueChanged, true);
            app.PCAButton.Text = 'PCA';
            app.PCAButton.Position = [51 235 79 46];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = HW4

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end