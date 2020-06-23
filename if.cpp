//Sam Rothstein
//May 2020

#include <iostream>
#include <vector>
#include <random>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

//TODO: #3 clean up all code and clean up file directory (see //NEW: tags)

//helper function header
void initializeArray(int arr[][1005], int n);
void initializeInputImage(int n,double arr[], double lower, double upper, double alpha);
void initializeVoltage(int n,double arr[]);
void initializeOutput(int output_to[], int output_from[], int edge, int N);
void spikeCout(int spike_cnt, int N);
void generalInfo(double tao, double Vr, double Vt, double sn);

int N = 1000; double sn = 0.0001; double t0 = 0.0; double tf = 10; double tau = 1;
int A[1005][1005]; double fp[1000000]; double v[1000000]; double t;
double observation[10000][10000]; bool spike[100000][100000]; int j; int spike_cnt;
int spikes; double tsp; int array_sparsity = 10; double alpha = 6.50;
int output_from[10000], output_to[10000]; int edge = 0; int key_neuron = 0;
double Vr = 0.0; double Vt = 1.0;

int main() {
    //update general information about the model
    generalInfo(tau, Vr, Vt,sn);
    double N_spike_rate[1000] = {};
    //initialize random
    srand (time(NULL));
    //initialize the connectivity matrix with semi-random values
    initializeArray(A, N);
    //Handle graph edges for demo
    initializeOutput(output_to, output_from, edge, N);
    // Initialize constant input vector for each neuron between 0.8 and 1.2;
    initializeInputImage(N, fp, 1.7, 2.3, alpha);
    //Step 1 "Randomly initialize processing neurons’ voltages"
    initializeVoltage(N,v);
    t = t0;
    //Open spike time file
    ofstream new_file_5;
    new_file_5.open("Spike_time.txt");
    //Open key neuron file
    ofstream new_file_6;
    new_file_6.open("Neuron_I.txt");
    //Open key neuron time file
    ofstream new_file_7;
    new_file_7.open("Neuron_I_T.txt");
    while (t < tf) {
      //Step 2 "Compute next spike time, for which vi = VT , using (15) for
      //some i = 1,...,m"
      //Step 3 "Determine the minimum spike time, tsp, and the spiking neuron
      //index, j."
      tsp = 1000; j = -1;
      for (int i = 0;i < N;i ++){
        double tsp_temp = - ( tau * log((Vt - fp[i]) / (v[i] - fp[i])));
        if (tsp_temp < tsp && tsp_temp > 0) {
          tsp = tsp_temp;
          j = i;
        }
      }
      //edge case with no spikes
      if (tsp == 1000 || j < 0){
        break;
      }
      //updating current time
      t += tsp;
      //Step 4 "Evaluate the voltage of all neurons at tsp, vi(tsp), for i ̸= j,
      //and set the voltage of the spiking neuron to the reset voltage,
      //vj(tsp) = VR."
      for (int i=0;i<N;i++){
        v[i] = (v[i] * exp(-(tsp/tau))) + (fp[i] * (1-exp(-(tsp/tau))));
        //write out the voltage of the key neuron
        if (i == key_neuron){
          new_file_6 << v[i] << endl;
          new_file_7 << t << endl;
        }
      }
      //Step 5 Add voltage jump (S/NA) to all neurons post-connected to
      //neuron j. If the voltage of any additional neuron reaches threshold, VT
      //repeat the step 4.
      //lock neurons that have already spiked
      bool lock[1000];
      for (int i=0;i<N;i++){
        lock[i] = false;
      }
      vector<int> waitlist; // This contains neurons that spike at the same time and haven't been dealt with.
      if (t != 10) {
        waitlist.push_back(j);
        spike_cnt++;
        new_file_5 << t << " ";
        //NEW:
        N_spike_rate[j] += 1;
      }
      // cout << waitlist.back() << " " << waitlist.size();
      while (!waitlist.empty()){
        int crt = waitlist.back();
        waitlist.pop_back();
        spikes++;
        lock[crt] = true;
        v[crt] = 0;
        for (int i=0;i<N;i++){
          if (i != crt && A[crt][i] == 1 && !lock[i]) {
            v[i] = v[i] + sn;
              if (v[i] >= 1) {
                waitlist.push_back(i);
                lock[i] = true;
              }
          }
        }
      }
      for (int i=0;i<N;i++){
        observation[spike_cnt][i] = v[i];
        if (lock[i] == true){
          spike[spike_cnt][i] = true;
        }
        else{
          spike[spike_cnt][i] = false;
        }
      }
    }//Step 6 Repeat steps 2-5 until t = tf .
    new_file_5.close();
    new_file_6.close();
    new_file_7.close();
    //update spike time and spike cout
    spikeCout(spike_cnt, N);

    //NEW:
    ofstream new_file_15;
    new_file_15.open("Neuron_Firing_Rates.txt");
    for (int i =0; i < N; i++){
      N_spike_rate[i] = N_spike_rate[i]/(double)tf;
      new_file_15 << N_spike_rate[i] << endl;
    }
    new_file_15.close();

    return 0;
}

//______________________________________________________________________________
//helper function body

void initializeArray(int arr[][1005], int n){
ofstream new_file_0;
new_file_0.open("Initial_Array.txt");
  for (int i=0;i<N;i++) {
    for (int j=0;j<N;j++){
        if (i == j){
          A[i][j] = 0;
        }
        else {
          int rnd = rand() % 100 + 1;
          if (rnd >= 100){
            A[i][j] = 1;
          }
          else{
            A[i][j]=0;
          }
        }
    new_file_0 << A[i][j] << endl;
    }
  }
new_file_0.close();
}


void initializeInputImage(int n,double arr[], double lower, double upper, double alpha){
  string line; ifstream myfile("Im.txt");
  if (myfile.is_open()){ // if there is an input file
    cout << "Opened a file -> generating input based on image" << endl;
    int index = 0;
    while (getline (myfile,line)){
      //IM.push_back(atof(line.c_str()));
      arr[index] = atof(line.c_str());
      index++;
    }
    myfile.close();
  }
  else{ //no input file
    cout << "Unable to open file -> generating input based on parameters" << endl;
    for (int i=0;i<n;i++){
        double r = ((double) rand() / (RAND_MAX));
        r = r * lower + upper; // [0.8,1.2];
        r *= alpha;
        arr[i] = r;
    }
  }
  //write input out to text document
  ofstream new_file_3;
  new_file_3.open("Initial_Input.txt");
  for (int i=0;i<n;i++){
      new_file_3 << arr[i] << endl;
  }
  new_file_3.close();
}

void initializeVoltage(int n,double arr[]){
  ofstream new_file_2;
  new_file_2.open ("Initial_Voltage.txt");
  for (int i=0;i<n;i++){
      arr[i] = ((double) rand() / (RAND_MAX));
      new_file_2 << arr[i] << endl;
  }
  new_file_2.close();
}

void initializeOutput(int output_to[], int output_from[], int edge, int N){
  for (int i=0;i<N;i++){
      for (int j=0;j<N;j++) {
          if (A[i][j] == 1) {
              edge++;
              output_from[edge] = i+1;
              output_to[edge] = j+1;
          }
      }
  }
  ofstream myfile;
  myfile.open ("Edges.txt");
  for (int i=1;i<=edge;i++) {
      myfile << output_from[i] << " " << output_to[i];
      if (i != edge) myfile << endl;
  }
  myfile.close();
  ofstream myfile_1;
  myfile_1.open ("Edges_num.txt");
  myfile_1 << edge;
  myfile_1.close();

  //writing the total number of nodes to a file
  ofstream new_file_1;
  new_file_1.open ("Total_node.txt");
  new_file_1 << N;
  new_file_1.close();
}

void spikeCout(int spike_cnt, int N){
  ofstream myfile_3;
  myfile_3.open ("Spike_count.txt");
  myfile_3 << spike_cnt;
  myfile_3.close();

  ofstream myfile_4;
  myfile_4.open ("Spikes.txt");
  for (int i=1;i<=spike_cnt;i++) {
      for (int j=0;j<N;j++){
          myfile_4 << spike[i][j] << " ";
      }
      if (i != spike_cnt) myfile_4 << endl;
  }
  myfile_4.close();

  ofstream new_file_6;
  new_file_6.open("FiringRate.txt");
  double average_num_firing_events = (double)spike_cnt / (double)N;
  double firing_rate;
  firing_rate = average_num_firing_events/tf;
  new_file_6 << "Alpha: " << alpha << endl;
  new_file_6 << firing_rate << endl;
  new_file_6.close();
}
void generalInfo(double tao, double Vr, double Vt, double sn){
  ofstream new_file_10;
  new_file_10.open("General_Inf.txt");
  new_file_10 << tau << endl;
  new_file_10 << Vr << endl;
  new_file_10 << Vt << endl;
  new_file_10 << sn << endl;
  new_file_10.close();
}
