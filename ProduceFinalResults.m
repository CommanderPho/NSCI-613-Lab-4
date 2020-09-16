% Produce Final Results: Lab 4
% Pho Hale
% After calling run_GolombNeuron_Ca0.m:

ProblemRunIndex = 1;
should_plot_other_ionic_currents = true;

% M-type K+ Current:
[IzTracesInfo.lim] = fnFindSeriesBounds(time_t, IzTraces, false);
IzTracesInfo.xlabel = 'time (msec)';
% IzTracesInfo.ylabel = 'M-type K+ current [\mu A/cm^{2}]';
IzTracesInfo.ylabel = 'I_{MK+} [\mu A/cm^{2}]';
IzTracesInfo.title = legend_strings;

% [FigH_K_Current] = fnPlotInteractiveSlider(time_t, IzTraces, FigH_K_Current, K_Current_SeriesConfigs);

% Membrane Voltage:
[voltageTracesInfo.lim] = fnFindSeriesBounds(time_t, voltageTraces, false);
voltageTracesInfo.xlabel = 'time [msec]';
% voltageTracesInfo.ylabel = 'membrane voltage [mV]';
voltageTracesInfo.ylabel = 'V_{m} [mV]';
voltageTracesInfo.title = legend_strings;
	

% Other Ionic Currents:
[OtherCurrentTracesInfo.lim] = fnFindSeriesBounds(time_t, [INaTraces, IKdrTraces, INaPTraces, IATraces], false);
OtherCurrentTracesInfo.xlabel = 'time (msec)';
OtherCurrentTracesInfo.ylabel = 'Other Ionic Currents [\mu A/cm^{2}]';
OtherCurrentTracesInfo.title = legend_strings;
OtherCurrentTracesInfo.legend = {'I_{Na+}','I_{Kdr}','I_{NaP}','I_{A}'};

% % INa Current:
% [INaTracesInfo.lim] = fnFindSeriesBounds(time_t, INaTraces, false);
% INaTracesInfo.xlabel = 'time (msec)';
% INaTracesInfo.ylabel = 'I_{Na+} [\mu A/cm^{2}]';
% INaTracesInfo.title = legend_strings;
% 
% % IKdr Current:
% [IKdrTracesInfo.lim] = fnFindSeriesBounds(time_t, IKdrTraces, false);
% IKdrTracesInfo.xlabel = 'time (msec)';
% IKdrTracesInfo.ylabel = 'I_{Kdr} [\mu A/cm^{2}]';
% IKdrTracesInfo.title = legend_strings;
% 
% % INaP Current:
% [INaPTracesInfo.lim] = fnFindSeriesBounds(time_t, INaPTraces, false);
% INaPTracesInfo.xlabel = 'time (msec)';
% INaPTracesInfo.ylabel = 'I_{NaP} [\mu A/cm^{2}]';
% INaPTracesInfo.title = legend_strings;
% 
% % IA Current:
% [IATracesInfo.lim] = fnFindSeriesBounds(time_t, IATraces, false);
% IATracesInfo.xlabel = 'time (msec)';
% IATracesInfo.ylabel = 'I_{A} [\mu A/cm^{2}]';
% IATracesInfo.title = legend_strings;


% [FigH_MembraneVoltage] = fnPlotInteractiveSlider(time_t, voltageTraces, FigH_MembraneVoltage, K_Current_SeriesConfigs);

%% Multi-voltage curve plots
% Multi-subplot version:

if ProblemRunIndex == 1
	searchPlotValues = [0.0, 0.4, 0.5, 1.2, 1.3, 1.5];
else
	searchPlotValues = [];
end

% Find the indicies of the relevant values
interestingPlotIndicies = zeros(length(searchPlotValues),1);
for i=1:length(searchPlotValues)
	found_indicies = find(abs(resultsTable.gM-searchPlotValues(i)) < 0.001);
	interestingPlotIndicies(i) = found_indicies;
end

% Do each plot:
num_active_indices = length(interestingPlotIndicies);
curr_voltage_strings = {};
curr_sub_strings = {};

% plot_mode = 'subplot';
plot_mode = 'tiled';
if should_plot_other_ionic_currents
	t = tiledlayout(num_active_indices,3,'TileSpacing','compact');
else
	t = tiledlayout(num_active_indices,2,'TileSpacing','compact');
end

% Loop through the results
for i=1:num_active_indices
	is_last_iteration = (num_active_indices == i);
	active_i = interestingPlotIndicies(i);
	curr_time_t_data = time_t{active_i};
	curr_Vm_data = voltageTraces{active_i};
	curr_Iz_data = IzTraces{active_i};
	
	% {INaTraces, IKdrTraces, INaPTraces, IATraces}
	if should_plot_other_ionic_currents
		curr_OtherIonicCurrents_data = {INaTraces{active_i}, IKdrTraces{active_i}, INaPTraces{active_i}, IATraces{active_i}};
	end

	%% Plot the current data
	
	voltageTracesInfo.title{active_i} = sprintf('%s for (%s)', 'Membrane Voltage', legend_strings{active_i});
	IzTracesInfo.title{active_i} = sprintf('%s for (%s)', 'M-type K+ Current', legend_strings{active_i});
	
	[voltageTracesInfo] = fnPlotData(num_active_indices, i, active_i, voltageTracesInfo, curr_time_t_data, curr_Vm_data, plot_mode);
	[IzTracesInfo] = fnPlotData(num_active_indices, i, active_i, IzTracesInfo, curr_time_t_data, curr_Iz_data, plot_mode);

	if should_plot_other_ionic_currents
		[OtherCurrentTracesInfo] = fnPlotData(num_active_indices, i, active_i, OtherCurrentTracesInfo, curr_time_t_data, curr_OtherIonicCurrents_data, plot_mode, OtherCurrentTracesInfo.legend);
	end
	
% 	OtherCurrentTracesInfos = {INaTracesInfo, IKdrTracesInfo, INaPTracesInfo, IATracesInfo};
% 	[OtherCurrentTracesInfos] = fnPlotData(num_active_indices, i, active_i, OtherCurrentTracesInfos, curr_time_t_data, curr_Iz_data, plot_mode);
	
	
	
end

linkaxes(voltageTracesInfo.ax_handle,'x');
linkaxes(IzTracesInfo.ax_handle,'x');
if should_plot_other_ionic_currents
	linkaxes(OtherCurrentTracesInfo.ax_handle,'x');
end

t.Padding = 'none';
t.TileSpacing = 'none';

base_export_path = '/Users/pho/Dropbox/Classes/Fall 2020/NSCI 613 - Neurophysiology and Computational Neuroscience/Lab 4/Results';

should_export_fig = true;
should_export_pdf = false;
should_export_png = true;


if ProblemRunIndex == 1
	fnSaveFigureForExport(t, fullfile(base_export_path,'2-1'), should_export_fig, should_export_pdf, should_export_png);
elseif ProblemRunIndex == 2
	fnSaveFigureForExport(t, fullfile(base_export_path,'2-2'), should_export_fig, should_export_pdf, should_export_png);
else
	error('Undefined case');
end


function [tracesInfo] = fnPlotData(num_active_indices, i, active_i, tracesInfo, x_data, y_data, plot_mode, legendNames)
	% non-reusable helper function to make plotting easier
	is_last_iteration = (num_active_indices == i);
	
	if strcmp(plot_mode, 'subplot')
		tracesInfo.ax_handle(i) = subplot(num_active_indices,1,i);
	elseif strcmp(plot_mode, 'tiled')
		tracesInfo.ax_handle(i) = nexttile;
	else
		error('Unhandled case!')
	end
	
	% Make a single-element cell array if the y_data is only a single
	% series
	if ~iscell(y_data)
		y_data = {y_data}; % Make a single-element cell array
	end
	
	% Loop through the subseries (if there are any) and plot them with hold
	% on.
	numSubSeries = length(y_data);
	hold off
	for j = 1:numSubSeries
		plot(x_data, y_data{j});
		hold on
	end
	% Add the things common to the plot:
	xlim([tracesInfo.lim(1),tracesInfo.lim(2)]);
	ylim([tracesInfo.lim(3),tracesInfo.lim(4)]);
	ylabel(tracesInfo.ylabel);
	title(tracesInfo.title{active_i});
	if exist('legendNames','var')
		legend(legendNames);
	end
	
	if is_last_iteration
		xlabel(tracesInfo.xlabel); % Only add the xlabel to the last (bottommost) plot
	else
		fnPostSubplotCleanup();
	end

end

function [export_result] = fnSaveFigureForExport(fig_h, figPath, should_export_fig, should_export_pdf, should_export_png)
	% fnSaveFigureForExport: performs export to disk of a provided figure.
	% Default values for optional parameters
	if ~exist('should_export_fig','var')
		should_export_fig = true;
	end
	if ~exist('should_export_png','var')
		should_export_png = false;
	end
	if ~exist('should_export_pdf','var')
		should_export_pdf = true;
	end
	

	% Perform requested exports
	if should_export_fig
		export_result.fig = [figPath '.fig'];
		savefig(export_result.fig)
	end
	if should_export_png
		export_result.png = [figPath '.png'];
		% Requires R2020a or later
		exportgraphics(fig_h, export_result.png,'Resolution',300)
	end
	if should_export_pdf
		export_result.pdf = [figPath '.pdf'];
		% Requires R2020a or later
		exportgraphics(fig_h, export_result.pdf,'ContentType','vector');
	end

end

function fnPostSubplotCleanup()
% 	set(gca,'XTick',[], 'YTick', []);
	set(gca,'XTick',[]);
	xlabel('');
end

