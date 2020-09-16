

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


% 	figure(2)
% 	plot(t,INa)
% 	hold on
% 	plot(t,IKdr)
% 	plot(t,INaP)
% 	plot(t,Iz)
% 	plot(t,IA)
% 	plot(tspan,[0 0],'k-')
% 	legend('INa','IKdr','INaP','IKM', 'IKA', 'zero')
% 	ylim([-5 5])
% 	xlabel('time (msec)')
% 	ylabel('ionic currents (microA/cm2)')
% 	hold off