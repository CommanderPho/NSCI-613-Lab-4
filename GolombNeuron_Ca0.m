function dvdt = GolombNeuron_Ca0(t,vars,gNaP, gM, basei,pulsei,t_on,t_off)
% Excitatory CA1 pyramidal model neuron, in 0 [Ca]
% from Golomb et al Contribution of Persistent Na  Current and M-Type K  Current 
%to Somatic Bursting in CA1 Pyramidal Cells: Combined Experimental
% and Modeling Study, J Neurophysiology, 2006



% variables:
    VVs = vars(1);
    hhs = vars(2);
    nns = vars(3);
    bbs = vars(4);
    zzs = vars(5);

% parameter values
    Cm=1.0;
    VNa=55.0; t_tauh=-40.5; t_taun=-27.0;
    thetaa=-50.0; sigmaa=20.0; thetab=-80.0; sigmab=-6.0; tauBs=15.0;
    sigmam=9.5; sigmah=-7.0; sigman=10.0; sigmaz=5.0; 

    % gNaP=0.3; gZ=1.0;
    gNa=35.0; gKdr=6.0; gL=0.05; 
    gA=1.4; 
    thetaz=-39.0; tauZs=75.0;
    phi=10.0;  thetah=-45.0;
    thetam=-30.0; thetan=-35.0; thetap=-47.0; sigmap=3.0;
    VK=-90.0; VL=-70.0;

% equation terms
    Iappx = basei + pulsei*heavyside(t-t_on)*heavyside(t_off-t);

    Minfs=GAMMAF(VVs,thetam,sigmam);
    Pinfs=GAMMAF(VVs,thetap,sigmap);
    Ainfs=GAMMAF(VVs,thetaa,sigmaa);

    INa=gNa*(Minfs^3)*hhs*(VVs-VNa);
    INaP=gNaP*Pinfs*(VVs-VNa);
    IKdr=gKdr*(nns^4)*(VVs-VK);
    IA=gA*Ainfs^3*bbs*(VVs-VK);
    Iz=gM*zzs*(VVs-VK);

% model equations
    dVVsdt = (-gL*(VVs-VL)-INa-INaP-IKdr-IA-Iz+Iappx)/Cm;
    dhhsdt = phi*(GAMMAF(VVs,thetah,sigmah)-hhs)/(1.0+7.5*GAMMAF(VVs,t_tauh,-6.0));
    dnnsdt = phi*(GAMMAF(VVs,thetan,sigman)-nns)/(1.0+5.0*GAMMAF(VVs,t_taun,-15.0));
    dbbsdt = (GAMMAF(VVs,thetab,sigmab)-bbs)/tauBs;
    dzzsdt = (GAMMAF(VVs,thetaz,sigmaz)-zzs)/tauZs;

dvdt = [dVVsdt; dhhsdt; dnnsdt; dbbsdt; dzzsdt];

end

function gf = GAMMAF(VV,theta,sigma)
    gf = 1.0/(1.0+exp(-(VV-theta)/sigma));
end


function hside=heavyside(x)
    if x >= 0
        hside = 1;
    else
        hside = 0;
    end
end
