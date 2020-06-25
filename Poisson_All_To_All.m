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
ylim([1 100]);
