
dt=0.01;
time=6; %(s)
numsteps = ceil(time/dt);

% input
I_baseline = 1;
I=I_baseline*ones(1,numsteps);
I(ceil(numsteps/2):end) = 0; % cancel form the middle

tau_rise=0.01; % rise time
tau_decay=1.0;
alpha=3;

% data history
eq1_history = zeros(1,numsteps);
eq2_history = zeros(1,numsteps);

% euler stepping
for t=2:numsteps,
    eq1_history(t) = eq1_history(t-1) + (dt/tau_decay)*(-eq1_history(t-1) + alpha*eq2_history(t-1));
    
    eq2_history(t) = eq2_history(t-1) + (dt/tau_rise)*(-eq2_history(t-1) + I(t));
end

%plot
figure;
hold on;
plot(eq1_history, '-r');
%plot((I_baseline*ones(1,numsteps)/alpha_1)*0.5, 'r');

plot(eq2_history, '-k');
%plot((I_baseline*ones(1,numsteps)/alpha_2)*0.5, 'k');

plot([ceil(numsteps/2) ceil(numsteps/2)],[0 alpha],'g-');

%set(gca, 'YScale','log');

