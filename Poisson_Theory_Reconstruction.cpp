//
//  main.cpp
//  Poisson_Theory_Reconstruction(Small f, (1,1/2), treat it as constant inputs for each neuron)
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
double f = 0.02; //TODO: Vary f values between .500 and .001
//Coupling strength
double S = 2;  //TODO: Keep S the same? -> yes
//Number of neurons
int N = 1000;
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
double tF = 10; //TODO: change this (3/4 runs per tF)
//Current time
double current_time;
// Image In
double arr[10005];
// \mu value
double v[10000];

int main(){
    srand (time(NULL));
    //Velocity of spikes
    string line; ifstream img("/Users/samuelrothstein/Desktop/Poisson_Driven/Im.txt");
    if (img.is_open()){ // if there is an input file
        int indx = 0;
        while (getline (img,line)){
            //IM.push_back(atof(line.c_str()));
            arr[indx] = atof(line.c_str());
            indx++;
        }
        img.close();
    }

    // Initialize Matrix B with 1/1000 chance of connection
    for (int i=0;i<N;i++)
        for (int j=0;j<10*N;j++){
            int rnd = rand() % 1000 + 1;
            if (rnd == 1000) b[i][j] = 1;
            else b[i][j]=0;
        }
    int Input_output_from[100000], Input_output_to[100000];
    int Input_edge = 0;
    for (int i=0;i<N;i++){
        for (int j=0;j<10*N;j++) {
            if (b[i][j] == 1) {
                Input_edge++;
                Input_output_from[Input_edge] = i+1;
                Input_output_to[Input_edge] = j+1;
            }
        }
    }
    ofstream Input_Edges;
    Input_Edges.open ("Poisson_Input_Edges.txt");
    for (int i=1;i<=Input_edge;i++) {
        Input_Edges << Input_output_from[i] << " " << Input_output_to[i];
        if (i != Input_edge) Input_Edges << endl;
    }
    Input_Edges.close();
    ofstream Input_Edges_num;
    Input_Edges_num.open ("Poisson_Input_Edges_num.txt");
    Input_Edges_num << Input_edge;
    Input_Edges_num.close();

    for (int i=0;i<N;i++){
        for (int j=0;j<10*N;j++){
            v[i] += (b[i][j] * arr[j]);
        }
        v[i] = v[i] / 10;
        v[i] = v[i] / f;
        // fp[i] *= (b*0.05); // Big For Loop
    }

    // Initialize Matrix A with 99% 0 and 1% 1;
    for (int i=0;i<N;i++)
        for (int j=0;j<N;j++){
            if (i == j) a[i][j] = 0;
            else {
                int rnd = rand() % 100 + 1;
                if (rnd == 100) a[i][j] = 1;
                else a[i][j]=0;
            }
        }
    int output_from[100000], output_to[100000];
    int edge = 0;
    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++) {
            if (a[i][j] == 1) {
                edge++;
                output_from[edge] = i+1;
                output_to[edge] = j+1;
            }
        }
    }
    ofstream Edges;
    Edges.open ("Poisson_Edges.txt");
    for (int i=1;i<=edge;i++) {
        Edges << output_from[i] << " " << output_to[i];
        if (i != edge) Edges << endl;
    }
    Edges.close();
    ofstream Edges_num;
    Edges_num.open ("Poisson_Edges_num.txt");
    Edges_num << edge;
    Edges_num.close();

    //Initialize voltages in voltage
    ofstream Initial_Voltage;
    Initial_Voltage.open ("Poisson_Initial_Voltage.txt");
    for (int i=0;i<N;i++){
        voltage[i] = ((double) rand() / (RAND_MAX)); //Value between 0 and 1
        Initial_Voltage << voltage[i] << endl;
    }
    Initial_Voltage.close();
    //ofstream Spikes_Observation;
    //Spikes_Observation.open ("Poisson_Spikes_Observation.txt");
    ofstream Individual_Spike_Count;
    Individual_Spike_Count.open("Poisson_Individual_Spike_Count.txt");

    priority_queue< pair<double,int>, vector<pair<double,int> >,
        greater<pair<double,int> > > pq;

    current_time = t0;
    for (int i=0;i<N;i++) {
        double tem = ((double) rand() / (RAND_MAX));
        // next_spike_time[i] = log (tem) * (-1 * (1/v));
        tem = log (tem) * (-1 * (1/v[i]));
        pq.push(make_pair(tem, i));
    }
    while(current_time < tF){

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
        }
        double tem = ((double) rand() / (RAND_MAX));
        tem = log (tem) * (-1 * (1/v[top.second]));
        pq.push(make_pair(top.first + tem, top.second));

    }
    cout << spike_count << " " << spikes << endl;

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
