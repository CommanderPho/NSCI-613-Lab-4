

%% Testing:

% M-type K+ Current:
FigH_K_Current = figure(6);
clf(FigH_K_Current)
K_Current_SeriesConfigs.xlabel = 'time (msec)';
K_Current_SeriesConfigs.ylabel = 'M-type K+ current (microA/cm2)';
K_Current_SeriesConfigs.title = legend_strings;

[FigH_K_Current] = fnPlotInteractiveSlider(time_t, IzTraces, FigH_K_Current, K_Current_SeriesConfigs);

% Membrane Voltage:
FigH_MembraneVoltage = figure(7);
clf(FigH_MembraneVoltage)
K_Current_SeriesConfigs.xlabel = 'time (msec)';
K_Current_SeriesConfigs.ylabel = 'membrane voltage (mV)';
K_Current_SeriesConfigs.title = legend_strings;

[FigH_MembraneVoltage] = fnPlotInteractiveSlider(time_t, voltageTraces, FigH_MembraneVoltage, K_Current_SeriesConfigs);










% 
% function sliderSin
% FigH = figure('position',[360 500 400 400]);
% axes('XLim', [0 4*pi], 'units','pixels', ...
%      'position',[100 50 200 200], 'NextPlot', 'add');
% x     = linspace(0, 4*pi, 400);
% y     = sin(x);
% LineH = plot(x,y);
% TextH = uicontrol('style','text',...
%     'position',[170 340 40 15]);
% SliderH = uicontrol('style','slider','position',[100 280 200 20],...
%     'min', 0, 'max', 4*pi);
% addlistener(SliderH, 'Value', 'PostSet', @callbackfn);
% movegui(FigH, 'center')
%     function callbackfn(source, eventdata)
%     num          = get(eventdata.AffectedObject, 'Value');
%     LineH.YData  = sin(num * x);
%     TextH.String = num2str(num);
%     end
%   end