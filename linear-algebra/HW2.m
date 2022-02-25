classdef HW2 < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure       matlab.ui.Figure
        edt_angle      matlab.ui.control.NumericEditField
        angleLabel     matlab.ui.control.Label
        edt_x          matlab.ui.control.NumericEditField
        xLabel         matlab.ui.control.Label
        btn_transform  matlab.ui.control.Button
        edt_y          matlab.ui.control.NumericEditField
        yLabel         matlab.ui.control.Label
        btn_open       matlab.ui.control.Button
        UIAxes         matlab.ui.control.UIAxes
    end

    
    properties (Access = private)
        img % image
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: btn_open
        function btn_openButtonPushed(app, event)
            app.img = imread('Test2.jpg'); % ??????? img
            imshow(app.img, 'parent', app.UIAxes); % ????? Axes ???
        end

        % Button pushed function: btn_transform
        function btn_transformButtonPushed(app, event)
            % ???????
            x = app.edt_x.Value;
            y = app.edt_y.Value;
            angle = app.edt_angle.Value;
            % ???????????
            R = [cosd(angle) -sind(angle) ; sind(angle) cosd(angle)];
            T = [x ; y];
            % ?????????????????
            img_new = zeros(size(app.img));
            img_new = uint8(img_new);
            % ????? pixel ?????
            for r = 1:size(img_new, 1)
                for c = 1:size(img_new, 2)
                    p_new = [c;r];
                    p = round(inv(R) * (p_new - T)); % ????? ? ??????round()?????
                    if(1 <= p(2) && p(2) <= size(app.img, 1) && ... 
                        1 <= p(1) && p(1) <= size(app.img, 2)) % ??????????????
                        img_new(r,c) = app.img(p(2),p(1)); % ???????????
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
            app.UIAxes.Position = [111 192 378 239];

            % Create btn_open
            app.btn_open = uibutton(app.UIFigure, 'push');
            app.btn_open.ButtonPushedFcn = createCallbackFcn(app, @btn_openButtonPushed, true);
            app.btn_open.Position = [59 56 100 22];
            app.btn_open.Text = {'open'; ''};

            % Create yLabel
            app.yLabel = uilabel(app.UIFigure);
            app.yLabel.HorizontalAlignment = 'right';
            app.yLabel.Position = [292 95 25 22];
            app.yLabel.Text = {'y'; ''};

            % Create edt_y
            app.edt_y = uieditfield(app.UIFigure, 'numeric');
            app.edt_y.Position = [332 95 100 22];

            % Create btn_transform
            app.btn_transform = uibutton(app.UIFigure, 'push');
            app.btn_transform.ButtonPushedFcn = createCallbackFcn(app, @btn_transformButtonPushed, true);
            app.btn_transform.Position = [509 56 100 22];
            app.btn_transform.Text = 'transform';

            % Create xLabel
            app.xLabel = uilabel(app.UIFigure);
            app.xLabel.HorizontalAlignment = 'right';
            app.xLabel.Position = [292 137 25 22];
            app.xLabel.Text = 'x';

            % Create edt_x
            app.edt_x = uieditfield(app.UIFigure, 'numeric');
            app.edt_x.Position = [332 137 100 22];

            % Create angleLabel
            app.angleLabel = uilabel(app.UIFigure);
            app.angleLabel.HorizontalAlignment = 'right';
            app.angleLabel.Position = [282 56 35 22];
            app.angleLabel.Text = {'angle'; ''};

            % Create edt_angle
            app.edt_angle = uieditfield(app.UIFigure, 'numeric');
            app.edt_angle.Position = [332 56 100 22];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = HW2

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