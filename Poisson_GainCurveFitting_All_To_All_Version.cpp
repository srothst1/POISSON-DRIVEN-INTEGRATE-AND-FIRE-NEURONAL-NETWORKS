//
//  main.cpp
//  Poisson_Driven
//
//  Created by Sam and Alex on 6/28/20.
//  All rights reserved.
//

// #include "/Users/haoyixuan/Desktop/cpp/test/stdc++.h"
// #include <bits/stdc++.h>
#include <iostream>
#include <vector>
#include <random>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <queue>
using namespace std;

//Coupling strength
double S = 2;
//Number of neurons
int N = 100;
//Firing at spike time?
//bool observation[100000][10000];
//Spike time
double spike_time[100000];
//Spike count at a given time
int spike_count = 0;
//Individual spike count
int individual_spike_count[10000];
//Total spikes
int spikes = 0;
//Starting voltages
double initial_voltage[10000];
//voltages
double voltage[10000];
//Starting time
double t0 = 0;
//Ending time
double tF = 11;
//Current time
double current_time;

// int N = 1000; double sn = 0.0001; double t0 = 0.0; double tf = 10; double tau = 1;
// int A[1005][1005]; double fp[1000000]; double v[1000000]; double t;
// double observation[10000][10000]; bool spike[100000][100000]; int j; int spike_cnt;
// int spikes; double tsp; int array_sparsity = 10; double alpha = 6.50;
// int output_from[10000], output_to[10000]; int edge = 0; int key_neuron = 0;
// double Vr = 0.0; double Vt = 1.0;

int main(){
    srand (time(NULL));
    // srand (1);
    //Initialize voltages in voltage
    //ofstream Initial_Voltage;
    //Initial_Voltage.open ("Poisson_Initial_Voltage.txt");
    //for (int i=0;i<N;i++){
    //    initial_voltage[i] = ((double) rand() / (RAND_MAX)); //Value between 0 and 1
    //    Initial_Voltage << initial_voltage[i] << endl;
    //}
    //Initial_Voltage.close();
    // ofstream Spikes_Observation;
    // Spikes_Observation.open ("Poisson_Spikes_Observation.txt");
    ofstream Spike_Count;
    ofstream Spikes;
    ofstream Spike_Time;
    Spike_Count.open ("Poisson_Spike_Count.txt");
    Spikes.open ("Poisson_Spikes.txt");
    Spike_Time.open ("Poisson_Spike_Time.txt");
    for (int BIG_LOOP=0; BIG_LOOP<100;BIG_LOOP++){
    for (int loop_fv=0;loop_fv<=50;loop_fv++){
        //Voltage jump caused by outside delta function spikes
        double f = 0.01;
        //Velocity of spikes
        double v = (0.5+(loop_fv)*0.02) / f;
        for (int i=0;i<N;i++){
            initial_voltage[i] = ((double) rand() / (RAND_MAX)); //Value between 0 and 1
            voltage[i] = initial_voltage[i];
            individual_spike_count[i] = 0;
        }
        for (int i=1;i<=spike_count;i++) spike_time[i] = 0;
        current_time = t0;
        priority_queue< pair<double,int>, vector<pair<double,int>>,
        greater<pair<double,int>> > pq;
        spikes = 0; spike_count = 0;
        for (int i=0;i<N;i++) {
            double tem = ((double) rand() / (RAND_MAX));
            // next_spike_time[i] = log (tem) * (-1 * (1/v));
            tem = log (tem) * (-1 * (1/v));
            pq.push(make_pair(tem, i));
        }
        while(current_time < tF){
            //construct N-1 array
            // long double next_spike_time[N];
            
            //each N_i is a randomly generated time (between 0 and inf)
            
            /*
            int smallest_index = -1;
            long double smallest_value = 100;
            //find "closest time" -> smallest N_i
            //Determine the next spike time
            
            for (int i = 0; i < N; i++){
                if (next_spike_time[i] < smallest_value){
                    smallest_index = i;
                    smallest_value = next_spike_time[i];
                }
            }
            current_time += smallest_value;
            voltage[smallest_index] += f;
            */
            //deduct the smallest value from all other N_k
            pair<double, int> top = pq.top();
            pq.pop();
            for (int i = 0; i < N; i++){
                voltage[i] = voltage[i] * exp(-top.first+current_time);
            }
            voltage[top.second] += f;
            current_time = top.first;
            if (voltage[top.second] > 1){
                if (current_time>1) spike_count++;
                spike_time[spike_count] = current_time;
                bool lock[10000]; // Lock the ones that have already spiked;
                for (int i=0;i<N;i++) lock[i] = false;
                vector<int> waitlist; // This contains neurons that spike at the same time and haven't been dealt with.
                waitlist.push_back(top.second);
                while (!waitlist.empty()){
                    int crt = waitlist.back();
                    waitlist.pop_back();
                    if (current_time>1) spikes++;
                    individual_spike_count[crt]++;
                    lock[crt] = true;
                    voltage[crt] = 0;
                    for (int i=0;i<N;i++){
                        if (i != crt && !lock[i]) {
                            voltage[i] = voltage[i] + ((double) S/ (double) N);
                            if (voltage[i] >= 1) {
                                waitlist.push_back(i);
                                lock[i] = true;
                            }
                        }
                    }
                }
                //Loop through
                //for (int i = 0; i < N; i++){
                //    if (lock[i]){
                //        observation[spike_count][i] = true;
                //    }
                //}
            }
            double tem = ((double) rand() / (RAND_MAX));
            tem = log (tem) * (-1 * (1/v));
            pq.push(make_pair(top.first + tem, top.second));
            
            //update whole system -> if fire -> go into "fire" loop
            //give a new random variable and add N_i to current_time
            //repeat
            
            // double tem = ((double) rand() / (RAND_MAX));
            // next_spike_time[smallest_index] = log (tem) * (-1 * (1/v));
        }
        cout << spike_count << " " << spikes << endl;
        Spike_Count << spike_count << endl;
        Spikes << spikes << endl;
        for (int i=1;i<=spike_count;i++){
            Spike_Time << spike_time[i] << endl;
        }
    }
    }
    //for (int i=1;i<=spike_count;i++)
    //    for (int j=0;j<N;j++) {
    //        Spikes_Observation << observation[i][j] << endl;
    //    }
    //Spikes_Observation.close();
    Spike_Count.close();
    Spikes.close();
    ofstream Neurons;
    Neurons.open ("Poisson_Neurons.txt");
    Neurons << N;
    Neurons.close();
    Spike_Time.close();
    return 0;
}
