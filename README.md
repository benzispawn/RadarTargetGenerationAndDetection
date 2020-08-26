#### 1. FMCW Waveform Design
Here we start do define some of the values as velocity, initial range, range resolution, max tange, speed of light and 
the FMCW waveform(Bandwidth, Chirp time, slope).

```Matlab

c = 3e8; %Speed of light
v = 20; %Velocity
si = 0; %Initial Range
RR = 1; %Range Resolution 
Rm = 200; % Max Range


%% FMCW Waveform Generation
B = c / 2*RR; %Bandwidth
Tchirp = 5.5*2*Rm/c; % chirp
slope = B/Tchirp; % Slope
```

#### 2. Target Generation and Detection
Set initial position and velocity. A loop to simulate the movement and we allocate all the values in matrix.

```Matlab
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples

Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));

for i=1:length(t)         
    r_t(i) = si + v*t(i);
    td(i) = 2*r_t(i)/c; 
    Tx(i) = cos(2*pi*(fc*t(i)+ slope*t(i)^2/2));
    Rx(i) = cos(2*pi*(fc*(t(i) - td(i)) + slope*(t(i) - td(i))^2/2));
    Mix(i) = Tx(i) * Rx(i);
end
```

#### 3. FFT Operation - 1D

We implemented de 1D FFT.

```Matlab
sig_fft = fft(Mix,Nr)./Nr;

sig_fft = abs(sig_fft);

sig_fft = sig_fft(1:Nr/2);

figure ('Name','Range from First FFT');
subplot(2,1,1)
plot(sig_fft);
title('FFT 1D');
axis ([0 200 0 1]);
xlabel('Range');
```
<img width="480" src="/home/benzi/radar_detect/fig_1.png">

#### 4. FFT Operation - 2D

```Matlab
Mix=reshape(Mix,[Nr,Nd]);

sig_fft2 = fft2(Mix,Nr,Nd);
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure('Name', 'Range from First FFT'),surf(doppler_axis,range_axis,RDM);
```
<img width="480" src="/home/benzi/radar_detect/fig_2.png">

#### 5. 2D CFAR
Setting the training and guard cell of the grid. 

```Matlab
Tcr = 10;
Tcd = 4;

Gcr = 5;
Gcd = 2;

offset = 1.4;
```

Declare the noise level, size of the grid, number of training cells and implementing the CFAR.

```Matlab
noise_level = zeros(Nr/2-2*(Tcd+Gcd), Nd-2*(Tcr+Gcr));
gridSize = (2*Tcr+2*Gcr+1)*(2*Tcd+2*Gcd+1);
numTrainCells = gridSize-(2*Gcr+1)*(2*Gcd+1);

CFAR_sig = zeros(size(RDM));
```

We loop with the grid values stablished above and we some the average values converting logarithmic to liner using the function db2pow(I am using Octave and I had to install a module for this step to work). At the end of every interaction we select specific values using the stablished threshold based on the offset above;

```Matlab
for j=1:Nd-2*(Tcr+Gcr)
    for i=1:Nr/2-2*(Tcd+Gcd)
      trainCells = db2pow(RDM(i:i+2*(Tcd+Gcd),j:j+2*(Gcr+Tcr)));
      trainCells(Tcd+1:end-Tcd,Tcr+1:end-Tcr) = 0;
      
      noise_level(i,j) = pow2db(sum(sum(trainCells))/numTrainCells);
      threshold = noise_level(i,j)*offset;
      if RDM(i+(Tcd+Gcd),j+(Tcd+Gcr))>threshold
          CFAR_sig(i+(Tcd+Gcd),j+(Tcd+Gcr)) = 1;
      else
          CFAR_sig(i+(Tcd+Gcd),j+(Tcd+Gcr)) = 0;
      end
  end
end

Finally we show the result. The values of training cell and guard cells were at the beginning 5 and 2 for Tcr and Tcd, 5 and 2 for Gcr and Gcd, respectively. We choose randomly at first but changing some of the values on training cells was enough to get a homogenous CFAR sign.  
```Matlab
figure('Name', 'CA-CFAR'),surf(doppler_axis,range_axis, CFAR_sig);
colorbar;

```
<img width="480" src="/home/benzi/radar_detect/fig_3.png">