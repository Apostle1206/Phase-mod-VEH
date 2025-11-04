% parameters
bom = 10; %\Omega
g = 1;
t = 0:0.002:(20*pi/bom); 
%som= \omega

% Calling the periodic_cosine_wave function
r = g*cos(bom*t);
% Initializing the functions
h = zeros(size(t));
b = zeros(size(t));
P = g * (mod(bom * t / pi, 2) - 1);

% Defining the piecewise functions with periodicity
for i = 1:length(t)
    % Calculating the periodic time variable
    t_mod = mod(t(i), 2*pi/bom);
    
    if t_mod < pi/bom
        h(i) = g;
        b(i) = g * (2 * bom * t_mod / pi - 1);
    else
        h(i) = -g;
        b(i) = -g * (2 * bom * t_mod / pi - 3);
    end
end



%subplots
figure;

% Plot of periodic cosine waveform
subplot(2, 2, 1);
plot(t, r, 'b','LineWidth', 2);
xlabel('t');
ylabel('P(t)');
title('Periodic Cosine Waveform');
ylim([-2 2]);
xlim([0 4]);


% Plot of Square waveform
subplot(2, 2, 2);

plot(t, h, 'b', 'LineWidth', 2);
xlabel('t');
ylabel('F(t)');
title('Square Waveform');
ylim([-2 2]);
xlim([0 4]);


% Plot of symmetric sawtooth waveform

subplot(2, 2, 3);
plot(t, b, 'b','LineWidth', 2);
xlabel('t');
ylabel('H(t)');
title('Symmetric Sawtooth Waveform');
ylim([-2 2]);
xlim([0 4]);


% Plot of asymmetric sawtooth waveform
subplot(2, 2, 4);
plot(t, P, 'b','LineWidth', 2);
xlabel('t');
ylabel('K(t)');
title('Asymmetric Sawtooth Waveform');
ylim([-2 2]);
xlim([0 4]);

