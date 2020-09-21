% Produce Final Results: Lab 4
% Pho Hale
% After calling run_GolombNeuron_Ca0.m:

% ProblemRunIndex = 1;

if ProblemRunIndex == 1
	should_plot_other_ionic_currents = false;
elseif ProblemRunIndex == 2
	should_plot_other_ionic_currents = false;
else
	error('Undefined case');
end

% M-type K+ Current:
[IzTracesInfo.lim] = fnFindSeriesBounds(time_t, IzTraces, false);
IzTracesInfo.xlabel = 'time (msec)';
% IzTracesInfo.ylabel = 'M-type K+ current [\mu A/cm^{2}]';
IzTracesInfo.ylabel = 'I_{MK+} [\mu A/cm^{2}]';
IzTracesInfo.title = legend_strings;

% Membrane Voltage:
[voltageTracesInfo.lim] = fnFindSeriesBounds(time_t, voltageTraces, false);
voltageTracesInfo.xlabel = 'time [msec]';
% voltageTracesInfo.ylabel = 'membrane voltage [mV]';
voltageTracesInfo.ylabel = 'V_{m} [mV]';
voltageTracesInfo.title = legend_strings;

% Combined M-type K+ Current and Persistent Na+ Currents:
[CombinedCurrentTracesInfo.lim] = fnFindSeriesBounds(time_t, [IzTraces, INaPTraces], false);
CombinedCurrentTracesInfo.xlabel = 'time (msec)';
CombinedCurrentTracesInfo.ylabel = 'Combined Ionic Currents [\mu A/cm^{2}]';
CombinedCurrentTracesInfo.title = legend_strings;
CombinedCurrentTracesInfo.legend = {'I_{MK+}','I_{NaP}'};



% Other Ionic Currents:
[OtherCurrentTracesInfo.lim] = fnFindSeriesBounds(time_t, [INaTraces, IKdrTraces, INaPTraces, IATraces], false);
OtherCurrentTracesInfo.xlabel = 'time (msec)';
OtherCurrentTracesInfo.ylabel = 'Other Ionic Currents [\mu A/cm^{2}]';
OtherCurrentTracesInfo.title = legend_strings;
OtherCurrentTracesInfo.legend = {'I_{Na+}','I_{Kdr}','I_{NaP}','I_{A}'};


%% Build the Frequency Plot:
fig_freqPlot = figure(1);
clf(fig_freqPlot);
hold off
if ProblemRunIndex == 1
	plot(resultsTable.gM, resultsTable.spikeFrequency_last); % Plot of gM vs. Frequency
	xlabel('M-type K+ conductance gM')
	title('Spike Frequency vs gM');
	xlim([resultsTable.gM(1) resultsTable.gM(end)]);

elseif ProblemRunIndex == 2
	plot(resultsTable.gNaP, resultsTable.spikeFrequency_last); % Plot of gNaP vs. Frequency
	xlabel('Persistent Na+ Conductance gNaP')
	title('Spike Frequency vs gNaP');
	xlim([resultsTable.gNaP(1) resultsTable.gNaP(end)]);
else
	error('Undefined case');
end
	
ylabel('Spike Frequency [Hz]');


%% Multi-subplot version:

if ProblemRunIndex == 1
% 	searchPlotValues = [0.0, 0.4, 0.5, 1.2, 1.3];
	searchPlotValues = [0.0, 0.4, 0.5, 0.6, 0.7, 1.2, 1.3];
elseif ProblemRunIndex == 2
	searchPlotValues = [0.06, 0.07, 0.14, 0.2];
else
	searchPlotValues = [];
end

% Find the indicies of the relevant values
interestingPlotIndicies = zeros(length(searchPlotValues),1);
for i=1:length(searchPlotValues)
	if ProblemRunIndex == 1
		found_indicies = find(abs(resultsTable.gM-searchPlotValues(i)) < 0.001);
	elseif ProblemRunIndex == 2
		found_indicies = find(abs(resultsTable.gNaP-searchPlotValues(i)) < 0.001);
	else
		error('Undefined case');
	end

	interestingPlotIndicies(i) = found_indicies;
end

% Generate filtered results:
filteredResultsTable = resultsTable(interestingPlotIndicies,:);
filteredISIs = spikeintervals(interestingPlotIndicies);


% Do each plot:
num_active_indices = length(interestingPlotIndicies);
curr_voltage_strings = {};
curr_sub_strings = {};

fig_mainPlot = figure(2);
% plot_mode = 'subplot';
plot_mode = 'tiled';
if should_plot_other_ionic_currents
	t = tiledlayout(fig_mainPlot, num_active_indices,3,'TileSpacing','compact');
else
	t = tiledlayout(fig_mainPlot, num_active_indices,2,'TileSpacing','compact');
end


% Loop through the results
for i=1:num_active_indices
	is_last_iteration = (num_active_indices == i);
	active_i = interestingPlotIndicies(i);
	curr_time_t_data = time_t{active_i};
	curr_Vm_data = voltageTraces{active_i};
	curr_Iz_data = IzTraces{active_i};

	if should_plot_other_ionic_currents
		curr_OtherIonicCurrents_data = {INaTraces{active_i}, IKdrTraces{active_i}, INaPTraces{active_i}, IATraces{active_i}};
	end

	%% Plot the current data
	voltageTracesInfo.title{active_i} = sprintf('%s for (%s)', 'Membrane Voltage', legend_strings{active_i});
	IzTracesInfo.title{active_i} = sprintf('%s for (%s)', 'M-type K+ Current', legend_strings{active_i});
	
	[voltageTracesInfo] = fnPlotData(num_active_indices, i, active_i, voltageTracesInfo, curr_time_t_data, curr_Vm_data, plot_mode);
	
	
	if ProblemRunIndex == 1
		% Problem 1 Mode Only:
		[IzTracesInfo] = fnPlotData(num_active_indices, i, active_i, IzTracesInfo, curr_time_t_data, curr_Iz_data, plot_mode);

	elseif ProblemRunIndex == 2
		% Problem 2 Mode Only:
		CombinedCurrentTracesInfo.title{active_i} = sprintf('%s for (%s)', 'Combined Currents', legend_strings{active_i});
		curr_CombinedCurrents_data = {curr_Iz_data, INaPTraces{active_i}};
		[CombinedCurrentTracesInfo] = fnPlotData(num_active_indices, i, active_i, CombinedCurrentTracesInfo, curr_time_t_data, curr_CombinedCurrents_data, plot_mode, CombinedCurrentTracesInfo.legend);
	else
		error('Invalid!')
	end

	% Plots the extra ionic currents as an additional column.
	if should_plot_other_ionic_currents
		[OtherCurrentTracesInfo] = fnPlotData(num_active_indices, i, active_i, OtherCurrentTracesInfo, curr_time_t_data, curr_OtherIonicCurrents_data, plot_mode, OtherCurrentTracesInfo.legend);
	end
	
	
% 	curr_ISI_data = spikeintervals{active_i};
% 	figure
% 	plot(curr_ISI_data)
% 	
% 	resultFrequencies(i) = resultsTable(active_i);
	
% 	spikeFrequency_last
% 	OtherCurrentTracesInfos = {INaTracesInfo, IKdrTracesInfo, INaPTracesInfo, IATracesInfo};
% 	[OtherCurrentTracesInfos] = fnPlotData(num_active_indices, i, active_i, OtherCurrentTracesInfos, curr_time_t_data, curr_Iz_data, plot_mode);
	
	
	
end

linkaxes(voltageTracesInfo.ax_handle,'x');

if ProblemRunIndex == 1
	linkaxes(IzTracesInfo.ax_handle,'x');
elseif ProblemRunIndex == 2
	linkaxes([voltageTracesInfo.ax_handle, CombinedCurrentTracesInfo.ax_handle],'x');
	linkaxes(CombinedCurrentTracesInfo.ax_handle,'xy');
else
	error('Invalid!')
end
	
if should_plot_other_ionic_currents
	linkaxes(OtherCurrentTracesInfo.ax_handle,'x');
end

t.Padding = 'none';
t.TileSpacing = 'none';


base_export_path = '/Users/pho/Dropbox/Classes/Fall 2020/NSCI 613 - Neurophysiology and Computational Neuroscience/Lab 4/Results';

should_export_fig = true;
should_export_eps = false;
should_export_pdf = false;
should_export_png = true;

if ProblemRunIndex == 1
	fnSaveFigureForExport(fig_mainPlot, fullfile(base_export_path,'2-1'), should_export_fig, should_export_eps, should_export_pdf, should_export_png);
	fnSaveFigureForExport(fig_freqPlot, fullfile(base_export_path,'2-2'), should_export_fig, should_export_eps, should_export_pdf, should_export_png);
elseif ProblemRunIndex == 2
	fnSaveFigureForExport(fig_mainPlot, fullfile(base_export_path,'3-1'), should_export_fig, should_export_eps, should_export_pdf, should_export_png);
	fnSaveFigureForExport(fig_freqPlot, fullfile(base_export_path,'3-2'), should_export_fig, should_export_eps, should_export_pdf, should_export_png);
	
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

function [export_result] = fnSaveFigureForExport(fig_h, figPath, should_export_fig, should_export_eps, should_export_pdf, should_export_png)
	% fnSaveFigureForExport: performs export to disk of a provided figure.
	% Position the figure:
% 	fig_h.Parent.OuterPosition = [0 0 4 6];
	
	% Default values for optional parameters
	if ~exist('should_export_fig','var')
		should_export_fig = true;
	end
	if ~exist('should_export_png','var')
		should_export_png = false;
	end
	if ~exist('should_export_eps','var')
		should_export_eps = true;
	end
	if ~exist('should_export_pdf','var')
		should_export_pdf = true;
		enable_vector_pdf_output = false; % Explicitly enable vector PDF output if that's desired. It's very slow
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
	
	if should_export_eps
		export_result.eps = [figPath '.eps'];
		exportgraphics(fig_h, export_result.eps)
	end
	
	if should_export_pdf
		export_result.pdf = [figPath '.pdf'];
		% Requires R2020a or later
		if enable_vector_pdf_output
			exportgraphics(fig_h, export_result.pdf,'ContentType','vector','BackgroundColor','none');
		else
			exportgraphics(fig_h, export_result.pdf,'BackgroundColor','none');
		end
	end

end

function fnPostSubplotCleanup()
% 	set(gca,'XTick',[], 'YTick', []);
	set(gca,'XTick',[]);
	xlabel('');
end

