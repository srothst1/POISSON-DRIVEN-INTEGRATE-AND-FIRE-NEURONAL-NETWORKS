% Image Pre-Processing
% load in image file
close all; clear all;
load cameraman100.jpg;
B = imread('cameraman100.jpg');
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
            tem = (i-1)*(2*j-1)/200
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

%% Reconstruction
% f is effective 1/(10*40) = 1/400; 1/40 already factored in image p;
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

% Matrix R = (1/10) * Matrix b * Matrix Cross_D_inv
R = b * Cross_D_inv;
R = 0.1 * R;
Poisson_Individual_Spike_Count = Poisson_Individual_Spike_Count / 10;
RHS = Poisson_Individual_Spike_Count + (1/2) * ones(Poisson_Neurons,1) - (1/Poisson_Neurons) * (a * Poisson_Individual_Spike_Count);


%TODO: add noise level
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
        if abs(l) > max_value
            max_value = abs(l);
            max_index = i;
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

