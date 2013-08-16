
dt=0.01;
time=2; %(s)
numsteps = ceil(time/dt);

% input
I_baseline = 1;
I=I_baseline*ones(1,numsteps);
I(ceil(numsteps/2):end) = 0; % cancel form the middle

alpha_1=1;
tau_1=0.1;

tau_rise = 0.1;
tau_decay = 1.0;
sigma_rise = 0.1*I_baseline;

alpha_2=0.0000001;
tau_2=(tau_1/alpha_1)*alpha_2; % typically alpha_1=1, so we dont really need that in there

% data history
eq1_history = zeros(1,numsteps);
eq2_history = zeros(1,numsteps);

% euler stepping
for t=2:numsteps,
    
    %% time constant swtich
    
    %if(I(t) > 0),
    %    tau = tau_rise;
    %else
    %    tau = tau_decay;
    %end
    
    tau = tau_rise + exp(-(I(t)^2)/(2*sigma_rise^2))*(tau_decay-tau_rise)
    
    eq1_history(t) = eq1_history(t-1) + (dt/tau)*(-eq1_history(t-1) + I(t));
    
    %% Coefficient control
    %eq1_history(t) = eq1_history(t-1) + (dt/tau_1)*(-alpha_1*eq1_history(t-1) + I(t));
    %eq2_history(t) = eq2_history(t-1) + (dt/tau_2)*(-alpha_2*eq2_history(t-1) + I(t));
end

%plot
figure;
hold on;
plot(eq1_history, '-r');
%plot((I_baseline*ones(1,numsteps)/alpha_1)*0.5, 'r');

%plot(eq2_history, '-k');
%plot((I_baseline*ones(1,numsteps)/alpha_2)*0.5, 'k');

%set(gca, 'YScale','log');

ylim([0 1]);
