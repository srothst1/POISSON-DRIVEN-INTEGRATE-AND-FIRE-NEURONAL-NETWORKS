%% Image Pre-Processing 1
% load in image file
close all; clear all;
load Orca.jpg;
B = imread('Orca.jpg');
% convert to useful form
if(ndims(B)==3)
    A = rgb2gray(B);
else
    A = double(B);
end
S = size(A);
Divide_cnt = 35; % normalize the grey-scale picture. 
imshow(A,[0 255]);
% reshape from matrix into vector, stacking columns under one another from
% left to right
k=1;  
C = zeros(S(1,1)*S(1,2),1);
for i = 1:S(1,2)
    for j = 1:S(1,1)
        C(k,1)=A(i,j) / Divide_cnt;
        k=k+1;
    end
end
save('/Users/samuelrothstein/Desktop/Poisson_Driven/Im.txt', 'C', '-ascii');

%% Fourier Basis
for i = 1 : 100
    for j = 1 : 100
        if i == 1
            tem = (i-1)*(2*j-1)/200;
            D(i,j) = sqrt(1/100) * cos(tem*pi);
        else
            tem = (i-1) * (2*j-1)/200;
            D(i,j) = sqrt(1/50) * cos(tem*pi);
        end
    end
end
Cross_D = [];
for i = 1 : 100
    InterM = [];
    for j = 1 : 100
        InterM = [InterM D(i,j)*D];
    end
    if (i == 1)
        Cross_D = InterM;
    else
        Cross_D = [Cross_D; InterM];
    end
end
Cross_D_inv = inv(Cross_D);
%mtx = dctmtx(100);
%Cross_D_inv =kron(mtx', mtx');

%% Gain Curve Fitting
load('Poisson_Input_Edges.txt');
load('Poisson_Input_Edges_num.txt');
load('Poisson_Edges.txt');
load('Poisson_Edges_num.txt');
load('Poisson_Individual_Spike_Count.txt');
load('Poisson_Neurons.txt');

for i = 1 : Poisson_Neurons
    for j = 1 : 11
        ISK(i,j) = Poisson_Individual_Spike_Count(1000*(j-1)+i) / 10;
    end
end

for i = 1:11
    plot((i+19)*0.05,ISK(345,i),'r.');
    hold on;
end
xlabel('Scalar Value','FontSize',16);
ylabel('Spike Count Per Unit Time','FontSize',16);
hold off;

%% Gain Curve Fitting Reconstruction (Base on fbp = slope * mu + bias)
for i = 1 : Poisson_Neurons
    for j = 1 : (Poisson_Neurons*10)
        b(i,j) = 0;
    end
end
for i = 1 : Poisson_Input_Edges_num
    b(Poisson_Input_Edges(i,1),Poisson_Input_Edges(i,2)) = 1;
end
RHS = (1/10) * b * C;
% bias = zeros(21,1000);
slope = zeros(21,1000);

%for loop3 = 1 : 21 % Big Loop Only 
%for i = 1 : Poisson_Neurons
%noise_level = 0.02 * (loop3-1);
for i = 1 : 1000
    for j = 1 : 11
        RHS_Fitting(j) = RHS(i) * (j+19) * 0.05;
    end
    % For individual neuron, Scaler * fbp = k * mu + b;
    % Start Fitting
    for k = 1 : 11
        Fitting_A(k,1) = ISK(i,k);
        % Fitting_A(k,1) = Fitting_A(k,1) * (1 + noise_level * randn);
        Fitting_A(k,2) = 1;
    end
    for k = 1 : 11
        Fitting_b(k) = RHS_Fitting(k);
    end
    if det(Fitting_A' * Fitting_A) == 0
        x = [0;mean(RHS_Fitting)];
    else
        x = inv(Fitting_A' * Fitting_A) * Fitting_A' * Fitting_b';
    end
    % slope(loop3,i) = x(1);
    % bias(loop3,i) = x(2);
    slope(i) = x(1);
    bias(i) = x(2);
end
%end % Big Loop Only 

%% Reconstuction of a new picture initialization 
load cameraman100.jpg;
Rec_B = imread('cameraman100.jpg');
% convert to useful form
if(ndims(Rec_B)==3)
    Rec_A = rgb2gray(Rec_B);
else
    Rec_A = double(Rec_B);
end
Rec_S = size(Rec_A);
imshow(Rec_A,[0 255]);
% reshape from matrix into vector, stacking columns under one another from
% left to right
k=1;  
Rec_C = zeros(Rec_S(1,1)*Rec_S(1,2),1);
for i = 1:Rec_S(1,2)
    for j = 1:Rec_S(1,1)
        Rec_C(k,1)=Rec_A(i,j) / Divide_cnt;
        k=k+1;
    end
end
save('/Users/samuelrothstein/Desktop/Poisson_Driven/Im.txt', 'Rec_C', '-ascii');

%% Reconstuction of a new picture (Base on fbp = slope * mu + bias)

load('Poisson_Input_Edges.txt');
load('Poisson_Input_Edges_num.txt');
load('Poisson_Edges.txt');
load('Poisson_Edges_num.txt');
load('Poisson_Individual_Spike_Count.txt');
load('Poisson_Neurons.txt');

% Reconstruct Matrix A
for i = 1 : Poisson_Neurons
    for j = 1 : Poisson_Neurons
        a(i,j) = 0;
    end
end
for i = 1 : Poisson_Edges_num
    a(Poisson_Edges(i,2),Poisson_Edges(i,1)) = 1;
end

% Reconstruct Matrix B
for i = 1 : Poisson_Neurons
    for j = 1 : (Poisson_Neurons*10)
        b(i,j) = 0;
    end
end
for i = 1 : Poisson_Input_Edges_num
    b(Poisson_Input_Edges(i,1),Poisson_Input_Edges(i,2)) = 1;
end
Poisson_Individual_Spike_Count = Poisson_Individual_Spike_Count / 10;

% error = zeros(21,4);
% for loop1 = 1:21 % Big Loop Only
% for loop2 = 1:4 % Big Loop Only
% noise_level = (loop1-1) * 0.02; % Big Loop Only
% Matrix R = (1/10) * Matrix b * Matrix Cross_D_inv
R = b * Cross_D_inv;
R = 0.1 * R; 
RHS = zeros(Poisson_Neurons, 1);
for i = 1 : Poisson_Neurons
    % Tem_Rec_Individual_spike_cnt(i) = Rec_Individual_Spike_Count(i) * (1 + noise_level * randn);
    % RHS(i) = Tem_Rec_Individual_spike_cnt(i) * slope(loop1,i) + bias(loop1,i);
    RHS(i) = Poisson_Individual_Spike_Count(i) * slope(i) + bias(i);
end

r1 = rand(1000,1); 
r1 = r1 .* 1; %TODO: update noise here
r1 = ones(1000,1) + r1;
RHS = RHS .* r1;

P_hat = zeros(10000,1);
cnt = 0;
while cnt < 250
    max_index = -1;
    max_value = -100;
    for i = 1 : 10000
        l = (R(:,i)' * RHS) / norm(R(:,i));
        if ~isnan(l)
            if abs(l) > max_value
                max_value = abs(l);
                max_index = i;
            end
        end
    end
    P_hat(max_index) = ((R(:,max_index)' * RHS) / norm(R(:,max_index))) / norm(R(:,max_index));
    RHS = RHS - P_hat(max_index)*R(:,max_index);
    cnt = cnt + 1;
end
image_theory = Cross_D_inv * P_hat;
image_theory = image_theory * Divide_cnt;
image_display = zeros(100,100);
for i = 1 : 100
    for j = 1 : 100
        image_display(i,j) = image_theory(100*(i-1)+j);
    end
end
%imshow(image_display,[0 255]);
Image_error = norm(image_display-A,'fro')/norm(A,'fro')


%Image_error = norm(image_display-Rec_A,'fro')/norm(Rec_A,'fro');


% error(loop1,loop2) = Image_error; % Big Loop Only
% end % Big Loop Only
% end % Big Loop Only