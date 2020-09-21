% code to run simulations of 
% Excitatory CA1 pyramidal model neuron, in 0 [Ca]
% from Golomb et al Contribution of Persistent Na  Current and M-Type K  Current 
%to Somatic Bursting in CA1 Pyramidal Cells: Combined Experimental
% and Modeling Study, J Neurophysiology, 2006

% The portion of the problem I'm working on.
% ProblemRunIndex = 1; % varying gM

clear all;
% ProblemRunIndex = 1; % varying gM
ProblemRunIndex = 2; % varying gNaP

% time span of simulation
tspan = [0 2000];

% initial conditions
VVs0=-71.81327;
hhs0=0.98786;
nns0=0.02457;
bbs0=0.203517;
zzs0=0.00141;
ICs = [VVs0; hhs0; nns0; bbs0; zzs0];

% set applied current pulse
basei = 0;   % baseline applied current
pulsei = 1;  % magnitude of current pulse above baseline
t_on = 100;    % time current pulse turns on
t_off = 2000;  % time current pulse turns off

% set conductances
if ProblemRunIndex == 1
	gNaP=0.0; % persistent Na+ conductance
	gM_vec = 0:0.1:1.5;  % M-type K+ conductance
	num_iterations = length(gM_vec);
elseif ProblemRunIndex == 2
	gM = 1; % M-type K+ conductance
	gNaP_vec = 0:0.01:0.7; % persistent Na+ conductance
	num_iterations = length(gNaP_vec);
else
	error('Unhandled')
end

% Pre-allocate:
spikeFrequency_last = zeros(num_iterations, 1);
spikeCounts = zeros(num_iterations, 1);



for i = 1:num_iterations
	
	if ProblemRunIndex == 1
		gM=gM_vec(i);   % M-type K+ conductance
% 		curr_description_string = sprintf('M-type K+ conductance gM = %.2g', gM);
		legend_strings{i} = sprintf('gM: %.2g', gM_vec(i));
	elseif ProblemRunIndex == 2
		gNaP=gNaP_vec(i);   % Persistent Na+ conductance
% 		curr_description_string = sprintf('Persistent Na+ conductance gNaP = %.2g', gNaP);
		legend_strings{i} = sprintf('gNaP: %.2g', gNaP_vec(i));
	else
		error('Unhandled')
	end
	
	% set up for simulation
	options = odeset('MaxStep',1);
	GNEquations_ftn = @(t,vars)GolombNeuron_Ca0(t, vars, gNaP, gM, basei, pulsei, t_on, t_off);

	% Call the ODE solver ode15s to numerically simulate
	[t,vars] = ode15s(GNEquations_ftn, tspan, ICs, options);

	% save variables
	VVs = vars(:,1);
	hhs = vars(:,2);
	nns = vars(:,3);
	bbs = vars(:,4);
	zzs = vars(:,5);

	% determine spike times and interspike intervals
	[peaks, indxs]=findpeaks(VVs,'MINPEAKHEIGHT',-10);
	spikeTimes{i}=t(indxs);
	spikeintervals{i}=diff(spikeTimes{i});
	spikeCounts(i) = length(spikeTimes{i});
	
	% Compute Frequency:
	% IPI: Inter-peak interval: the duration (in [ms]) between the peaks times.
	% ISI: Inter-spike interval: the duration (in [ms]) between the spike times.
	if ~isempty(spikeintervals{i})
		% Using Two Computation Styles:
		% 1) Dr. Booth uses the last IPI to define the frequency, so I've added this as an alternative frequency metric.
		last_IPI_seconds = spikeintervals{i}(end) / 1000; % Divide by 1000 to convert from [ms] to [sec]
		spikeFrequency_last(i) = 1 ./ last_IPI_seconds;
	else
		spikeFrequency_last(i) = NaN;
	end
	
	% ionic currents
	% inactivating Na
	gNa=35.0; VNa=55.0; thetam=-30.0; sigmam=9.5;
	Minfs=1.0./(1.0+exp(-(VVs-thetam)/sigmam));
	INa=gNa*(Minfs.^3).*hhs.*(VVs-VNa);

	% K+ delayed rectifier
	gKdr=6.0; VK=-90.0;
	IKdr=gKdr*(nns.^4).*(VVs-VK);

	% persistent Na 
	thetap=-47.0; sigmap=3.0;
	Pinfs=1.0./(1.0+exp(-(VVs-thetap)/sigmap));
	INaP=gNaP*Pinfs.*(VVs-VNa);

	% K+ M current
	Iz=gM*zzs.*(VVs-VK);

	% K+ A current
	gA=1.4; thetaa=-50.0; sigmaa=20.0;
	Ainfs=1.0./(1.0+exp(-(VVs-thetaa)/sigmaa));
	IA=gA*Ainfs.^3.*bbs.*(VVs-VK);

	% Save the time and voltage curves for each stimulation current in case we want to plot them.
	time_t{i} = t;
	voltageTraces{i} = VVs;
	IzTraces{i} = Iz;
	
	INaTraces{i} = INa;
	IKdrTraces{i} = IKdr;
	INaPTraces{i} = INaP;
	IATraces{i} = IA;
	
end

% Build a table from the results for easy browsing of the different variables as a function of iteration and current.

if ProblemRunIndex == 1
	resultsTable = table(gM_vec', spikeCounts, spikeFrequency_last,'VariableNames',{'gM','spikeCounts','spikeFrequency_last'});
elseif ProblemRunIndex == 2
	resultsTable = table(gNaP_vec', spikeCounts, spikeFrequency_last,'VariableNames',{'gNaP','spikeCounts','spikeFrequency_last'});
else
	error('Unhandled')
end

