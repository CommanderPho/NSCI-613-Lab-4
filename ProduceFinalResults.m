% Produce Final Results: Lab 4
% Pho Hale
% After calling run_GolombNeuron_Ca0.m:


ProblemRunIndex = 1;

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
t = tiledlayout(num_active_indices,2,'TileSpacing','compact');

% Loop through the results
for i=1:num_active_indices
	is_last_iteration = (num_active_indices == i);
	active_i = interestingPlotIndicies(i);
	curr_time_t_data = time_t{active_i};
	curr_Vm_data = voltageTraces{active_i};
	curr_Iz_data = IzTraces{active_i};
	
	% Plot the current data
	if strcmp(plot_mode, 'subplot')
		voltageTracesInfo.ax_handle(i) = subplot(num_active_indices,1,i);
	elseif strcmp(plot_mode, 'tiled')
% 		nexttile(t,i)
		voltageTracesInfo.ax_handle(i) = nexttile;
	else
		error('Unhandled case!')
	end
	curr_plot_handle = plot(curr_time_t_data, curr_Vm_data);
	xlim([voltageTracesInfo.lim(1),voltageTracesInfo.lim(2)]);
	ylim([voltageTracesInfo.lim(3),voltageTracesInfo.lim(4)]);
	
	ylabel(voltageTracesInfo.ylabel);
	voltageTracesInfo.title{active_i} = sprintf('%s for (%s)', 'Membrane Voltage', legend_strings{active_i});
	title(voltageTracesInfo.title{active_i});
	if is_last_iteration
		xlabel(voltageTracesInfo.xlabel);
	else
		fnPostSubplotCleanup();
	end
	
	if strcmp(plot_mode, 'subplot')
		IzTracesInfo.ax_handle(i) = subplot(num_active_indices,1,i);
	elseif strcmp(plot_mode, 'tiled')
		IzTracesInfo.ax_handle(i) = nexttile;
	else
		error('Unhandled case!')
	end
	
	curr_plot_handle = plot(curr_time_t_data, curr_Iz_data);
	xlim([IzTracesInfo.lim(1),IzTracesInfo.lim(2)]);
	ylim([IzTracesInfo.lim(3),IzTracesInfo.lim(4)]);
	xlabel(IzTracesInfo.xlabel);
	ylabel(IzTracesInfo.ylabel);
	IzTracesInfo.title{active_i} = sprintf('%s for (%s)', 'M-type K+ Current', legend_strings{active_i});
	title(IzTracesInfo.title{active_i});
	if is_last_iteration
		xlabel(voltageTracesInfo.xlabel);
	else
		fnPostSubplotCleanup();
	end
	
	
end

linkaxes(voltageTracesInfo.ax_handle,'x');
linkaxes(IzTracesInfo.ax_handle,'x');

t.Padding = 'none';
t.TileSpacing = 'none';

base_export_path = '/Users/pho/Dropbox/Classes/Fall 2020/NSCI 613 - Neurophysiology and Computational Neuroscience/Lab 4/Results';

should_export_fig = true;
should_export_pdf = false;
should_export_png = true;

fnSaveFigureForExport(t, fullfile(base_export_path,'2-1'), should_export_fig, should_export_pdf, should_export_png);


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

