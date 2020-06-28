%% Neurons Spiking Visual
close all;clear all;
load('Poisson_Spikes_Observation.txt');
load('Poisson_Spike_Time.txt');
load('Poisson_Spike_Count.txt');
load('Poisson_Neurons.txt');
load('Poisson_Spikes.txt');
for i = 1 : Poisson_Spike_Count
    for j = 1 : Poisson_Neurons
        y(i,j) = Poisson_Spikes_Observation(Poisson_Neurons*(i-1)+j);
        if (y(i,j) == 1)
            plot(Poisson_Spike_Time(i),j,'r.');
            hold on;
        end
    end
end
xlim([0 10]);
ylim([1 Poisson_Neurons]);

%% All-To-All Neuron Gain Curve Fitting
close all;clear all;
load('Poisson_Spikes_Observation.txt');
load('Poisson_Spike_Time.txt');
load('Poisson_Spike_Count.txt');
load('Poisson_Neurons.txt');
load('Poisson_Spikes.txt');
Avg = zeros(1,51);
for i = 1: 51
    for j = 1:50
        Avg(i) = Avg(i)+Poisson_Spikes((j-1)*51+i);
    end
    Avg(i) = Avg(i)/(j*Poisson_Neurons*10);
end

for i = 1:51
    plot(0.5+(i-1)*0.02,Avg(i),'r.');
    hold on;
end
