classdef HW1 < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                 matlab.ui.Figure
        ImageThresholdingButton  matlab.ui.control.Button
        ThresholdEditField       matlab.ui.control.NumericEditField
        ThresholdEditFieldLabel  matlab.ui.control.Label
        OpenButton               matlab.ui.control.Button
        UIAxes                   matlab.ui.control.UIAxes
    end

    
    properties (Access = private)
        img % Image
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: OpenButton
        function OpenButtonPushed(app, event)
            app.img = imread('Test1.jpg');              % ??????? img
            imshow(app.img, 'parent', app.UIAxes);      % ??????? UIAxes ?
        end

        % Button pushed function: ImageThresholdingButton
        function ImageThresholdingButtonPushed(app, event)
            threshold = app.ThresholdEditField.Value;   % ???????
            img_new = app.img;                          % ??????????????
            for i = 1:size(app.img, 1)                  % size()????????size(img, 1) ?????
                   for j = 1:size(app.img, 2)           % size(img, 2) ????
                        if(app.img(i,j) > threshold)    % ?????????????????
                            img_new(i,j) = 255;
                        else                            % ??????????
                            img_new(i,j) = 0;
                        end
                   end
            end
            imshow(img_new, 'parent', app.UIAxes);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 640 480];
            app.UIFigure.Name = 'MATLAB App';

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            title(app.UIAxes, 'Title')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.Position = [82 218 300 185];

            % Create OpenButton
            app.OpenButton = uibutton(app.UIFigure, 'push');
            app.OpenButton.ButtonPushedFcn = createCallbackFcn(app, @OpenButtonPushed, true);
            app.OpenButton.Position = [82 157 100 22];
            app.OpenButton.Text = 'Open';

            % Create ThresholdEditFieldLabel
            app.ThresholdEditFieldLabel = uilabel(app.UIFigure);
            app.ThresholdEditFieldLabel.HorizontalAlignment = 'right';
            app.ThresholdEditFieldLabel.Position = [232 157 59 22];
            app.ThresholdEditFieldLabel.Text = 'Threshold';

            % Create ThresholdEditField
            app.ThresholdEditField = uieditfield(app.UIFigure, 'numeric');
            app.ThresholdEditField.Position = [306 157 100 22];

            % Create ImageThresholdingButton
            app.ImageThresholdingButton = uibutton(app.UIFigure, 'push');
            app.ImageThresholdingButton.ButtonPushedFcn = createCallbackFcn(app, @ImageThresholdingButtonPushed, true);
            app.ImageThresholdingButton.Position = [459 157 121 22];
            app.ImageThresholdingButton.Text = 'Image Thresholding';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = HW1

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
