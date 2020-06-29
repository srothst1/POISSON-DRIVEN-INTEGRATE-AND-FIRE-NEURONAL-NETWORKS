//
//  main.cpp
//  Poisson_Gain_Curve_Fitting_Reconstruction FOR TESTING PURPOSE
//
//  Created by Sam and Alex on 6/29/20.
//  Copyright © 2020 Haoyi Xuan. All rights reserved.
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

//Voltage jump caused by outside delta function spikes
int a[10005][10005];
int b[1005][10005];
double f;
//Coupling strength
double S = 2;
//Number of neurons
int N = 1000;
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
double voltage[10000];
//Starting time
double t0 = 0;
//Ending time
double tF = 10;
//Current time
double current_time;
// Image In
double arr[10005];
// \mu value
double v[10000];
int Input_Edges_num,Edges_num;


// int N = 1000; double sn = 0.0001; double t0 = 0.0; double tf = 10; double tau = 1;
// int A[1005][1005]; double fp[1000000]; double v[1000000]; double t;
// double observation[10000][10000]; bool spike[100000][100000]; int j; int spike_cnt;
// int spikes; double tsp; int array_sparsity = 10; double alpha = 6.50;
// int output_from[10000], output_to[10000]; int edge = 0; int key_neuron = 0;
// double Vr = 0.0; double Vt = 1.0;

int main(){
    srand (time(NULL));
    // srand (1);
    //Velocity of spikes
    
    string line; ifstream img("/Users/haoyixuan/Desktop/Poisson_Driven/Rec_Im.txt");
    if (img.is_open()){ // if there is an input file
    // cout << "Opened a file -> generating input based on image" << endl;
        int indx = 0;
        while (getline (img,line)){
            //IM.push_back(atof(line.c_str()));
            arr[indx] = atof(line.c_str());
            indx++;
        }
        img.close();
    }
    ifstream Input_f("/Users/haoyixuan/Desktop/Poisson_Driven/Poisson_f.txt");
    Input_f >> f;
    Input_f.close();
    ifstream Input_E_N("/Users/haoyixuan/Desktop/Poisson_Driven/Poisson_Input_Edges_num.txt");
    Input_E_N >> Input_Edges_num;
    Input_E_N.close();
    cout << Input_Edges_num << endl;
    for (int i=0;i<N;i++)
        for (int j=0;j<10*N;j++) {
            b[i][j] = 0;
        }
    ifstream Input_E("/Users/haoyixuan/Desktop/Poisson_Driven/Poisson_Input_Edges.txt");
    for (int i=0;i<Input_Edges_num;i++){
        int j,k;
        Input_E >> j >> k;
        b[j-1][k-1] = 1;
    }
       
    ifstream E_N("/Users/haoyixuan/Desktop/Poisson_Driven/Poisson_Edges_num.txt");
    E_N >> Edges_num;
    E_N.close();
    cout << Edges_num << endl;
    for (int i=0;i<N;i++)
        for (int j=0;j<N;j++) {
            a[i][j] = 0;
        }
    ifstream E("/Users/haoyixuan/Desktop/Poisson_Driven/Poisson_Edges.txt");
    for (int i=0;i<Edges_num;i++){
        int j,k;
        E >> j >> k;
        a[j-1][k-1] = 1;
    }
    ifstream Initial_V("/Users/haoyixuan/Desktop/Poisson_Driven/Poisson_Initial_volt.txt");
    for (int i=0;i<N;i++){
        Initial_V >> voltage[i];
    }
    
    for (int i=0;i<N;i++){
        for (int j=0;j<10*N;j++){
            v[i] += (b[i][j] * arr[j]);
        }
        v[i] = v[i] / 10;
        v[i] = v[i] / f;
        // fp[i] *= (b*0.05); // Big For Loop
    }
    
    //ofstream Spikes_Observation;
    //Spikes_Observation.open ("Poisson_Spikes_Observation.txt");
    ofstream Individual_Spike_Count;
    Individual_Spike_Count.open("Rec_Individual_Spike_Count.txt");
    
    priority_queue< pair<double,int>, vector<pair<double,int>>,
        greater<pair<double,int>> > pq;
    
    current_time = t0;
    for (int i=0;i<N;i++) {
        double tem = ((double) rand() / (RAND_MAX));
        // next_spike_time[i] = log (tem) * (-1 * (1/v));
        tem = log (tem) * (-1 * (1/v[i]));
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
            spike_count++;
            spike_time[spike_count] = current_time;
            bool lock[10000]; // Lock the ones that have already spiked;
            for (int i=0;i<N;i++) lock[i] = false;
            vector<int> waitlist; // This contains neurons that spike at the same time and haven't been dealt with.
            waitlist.push_back(top.second);
            while (!waitlist.empty()){
                int crt = waitlist.back();
                waitlist.pop_back();
                spikes++;
                individual_spike_count[crt]++;
                lock[crt] = true;
                voltage[crt] = 0;
                for (int i=0;i<N;i++){
                    if (i != crt && a[crt][i] && !lock[i]) {
                        voltage[i] = voltage[i] + ((double) S/ (double) N);
                        if (voltage[i] >= 1) {
                            waitlist.push_back(i);
                            lock[i] = true;
                        }
                    }
                }
            }
            //Loop through
            /*
            for (int i = 0; i < N; i++){
                if (lock[i]){
                    observation[spike_count][i] = true;
                }
            }
            */
        }
        double tem = ((double) rand() / (RAND_MAX));
        tem = log (tem) * (-1 * (1/v[top.second]));
        pq.push(make_pair(top.first + tem, top.second));
        
        //update whole system -> if fire -> go into "fire" loop
        //give a new random variable and add N_i to current_time
        //repeat
        
        // double tem = ((double) rand() / (RAND_MAX));
        // next_spike_time[smallest_index] = log (tem) * (-1 * (1/v));
    }
    cout << spike_count << " " << spikes << endl;
    /*
    for (int i=1;i<=spike_count;i++)
        for (int j=0;j<N;j++) {
            Spikes_Observation << observation[i][j] << endl;
        }
    */
    for (int j=0;j<N;j++) {
        Individual_Spike_Count << individual_spike_count[j] << endl;
    }
    Individual_Spike_Count.close();
    // Spikes_Observation.close();
    ofstream Spike_Count;
    Spike_Count.open ("Poisson_Spike_Count.txt");
    Spike_Count << spike_count;
    Spike_Count.close();
    ofstream Spikes;
    Spikes.open ("Poisson_Spikes.txt");
    Spikes << spikes;
    Spikes.close();
    ofstream Neurons;
    Neurons.open ("Poisson_Neurons.txt");
    Neurons << N;
    Neurons.close();
    ofstream Spike_Time;
    Spike_Time.open ("Poisson_Spike_Time.txt");
    for (int i=1;i<=spike_count;i++){
        Spike_Time << spike_time[i] << endl;
    }
    Spike_Time.close();
    return 0;
}
