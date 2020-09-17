

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

% Other Ionic Currents:
FigH_OtherIonicCurrents = figure(8);
clf(FigH_OtherIonicCurrents)
[OtherCurrentTracesInfo.lim] = fnFindSeriesBounds(time_t, [INaTraces, IKdrTraces, INaPTraces, IATraces], false);
OtherCurrentTracesInfo.xlabel = 'time (msec)';
OtherCurrentTracesInfo.ylabel = 'Other Ionic Currents [\mu A/cm^{2}]';
OtherCurrentTracesInfo.title = legend_strings;
OtherCurrentTracesInfo.legend = {'I_{Na+}','I_{Kdr}','I_{NaP}','I_{A}'};

numIterations = length(INaTraces);
for i = 1:numIterations
	curr_OtherIonicCurrents_data{i} = {INaTraces{i}, IKdrTraces{i}, INaPTraces{i}, IATraces{i}};
end

% curr_OtherIonicCurrents_data = {INaTraces, IKdrTraces, INaPTraces, IATraces};
[FigH_OtherIonicCurrents] = fnPlotInteractiveSlider(time_t, curr_OtherIonicCurrents_data, FigH_OtherIonicCurrents, OtherCurrentTracesInfo);
