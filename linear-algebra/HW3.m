classdef HW3 < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure               matlab.ui.Figure
        btn_regression         matlab.ui.control.Button
        btn_open               matlab.ui.control.Button
        UIAxes_day2_corrected  matlab.ui.control.UIAxes
        UIAxes_day2            matlab.ui.control.UIAxes
        UIAxes_combained       matlab.ui.control.UIAxes
        UIAxes_scatter_plot    matlab.ui.control.UIAxes
        UIAxes_day1            matlab.ui.control.UIAxes
    end

    
    properties (Access = private)
        img1 % day1
        img2 % day2
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: btn_open
        function btn_openButtonPushed(app, event)
            app.img1 = imread('Satellite Image Day1 (B1).tif');
            imshow(app.img1, 'parent', app.UIAxes_day1);
            app.img2 = imread('Satellite Image Day2 (B1).tif');
            imshow(app.img2, 'parent', app.UIAxes_day2);
            
            % scatter(ax, x, y, sz, 'filled')??? x, y ????????ax ???? Axes?sz ??????'filled' ???????
            scatter(app.UIAxes_scatter_plot, app.img2(:), app.img1(:), 10, 'filled');
            % axis(ax, [xmin xmax ymin ymax])????????
            axis(app.UIAxes_scatter_plot, [0 255 0 255]);
            % xticks(), yticks()??? x, y ????0 ? 250 ? 50 ?????
            xticks(app.UIAxes_scatter_plot, 0:50:250);
            yticks(app.UIAxes_scatter_plot, 0:50:250);
            % xlabel(), ylabel()??? x, y ???
            xlabel(app.UIAxes_scatter_plot, 'Satellite Image Day2 (B1)');
            ylabel(app.UIAxes_scatter_plot, 'Satellite Image Day1 (B1)');
        end

        % Button pushed function: btn_regression
        function btn_regressionButtonPushed(app, event)
            % ?? y ??? C ??
            y = double(app.img1(:));
            % ?outlier??????
            thres_map = y > 100;
            y(thres_map) = 0;
            v1 = ones(size(y, 1), 1);
            img2_thrs = app.img2;
            img2_thrs(thres_map) = 0;
            v2 = double(img2_thrs(:));
            C = [v1 v2];
            % ??
            a = inv(C' * C) * C' * y;
            display(a)
            % y = Ca?????DAY2??????thresholding??
            img2_new = [v1 double(app.img2(:))] * a;
            % ?? Day2 ??
            img2_new = uint8(img2_new);
            img2_new = reshape(img2_new, size(app.img2));
            % ??Day1?Day2
            img_comb = app.img1;
            img_comb(floor(size(img2_new, 1)/2):end, :) = img2_new(floor(size(img2_new, 1)/2):end, :);
            imshow(img2_new, 'parent', app.UIAxes_day2_corrected);
            imshow(img_comb, 'parent', app.UIAxes_combained);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 945 502];
            app.UIFigure.Name = 'MATLAB App';

            % Create UIAxes_day1
            app.UIAxes_day1 = uiaxes(app.UIFigure);
            title(app.UIAxes_day1, 'Day1')
            app.UIAxes_day1.Position = [1 291 300 185];

            % Create UIAxes_scatter_plot
            app.UIAxes_scatter_plot = uiaxes(app.UIFigure);
            title(app.UIAxes_scatter_plot, 'Scatter Plot')
            app.UIAxes_scatter_plot.Position = [308 291 300 185];

            % Create UIAxes_combained
            app.UIAxes_combained = uiaxes(app.UIFigure);
            title(app.UIAxes_combained, {'Combained'; ''})
            app.UIAxes_combained.Position = [615 291 300 185];

            % Create UIAxes_day2
            app.UIAxes_day2 = uiaxes(app.UIFigure);
            title(app.UIAxes_day2, 'Day2')
            app.UIAxes_day2.Position = [308 62 300 185];

            % Create UIAxes_day2_corrected
            app.UIAxes_day2_corrected = uiaxes(app.UIFigure);
            title(app.UIAxes_day2_corrected, 'Day2 Corrected')
            app.UIAxes_day2_corrected.Position = [615 62 300 185];

            % Create btn_open
            app.btn_open = uibutton(app.UIFigure, 'push');
            app.btn_open.ButtonPushedFcn = createCallbackFcn(app, @btn_openButtonPushed, true);
            app.btn_open.Position = [45 143 100 22];
            app.btn_open.Text = {'Open'; ''};

            % Create btn_regression
            app.btn_regression = uibutton(app.UIFigure, 'push');
            app.btn_regression.ButtonPushedFcn = createCallbackFcn(app, @btn_regressionButtonPushed, true);
            app.btn_regression.Position = [184 143 100 22];
            app.btn_regression.Text = 'Regression';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = HW3

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