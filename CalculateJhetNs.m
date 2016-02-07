% Calculate Jhet and Ns from input files from Freezing detection software
% input file contains freezing temperature in celcius of droplets of water in oil as
% well as droplets diameter in microns.
% file name of freezing data from freezing detection software
% waterActivity (in numbers from 0 to 1)
% meltingTemperature (K) assumed a single number in which the droplets froze
% coolingRate (in Kelvin/min )
% surfaceArea from BET measurement ( meter^2/gm )
% concentration (% weight)
function CalculateJhetNs(inputFile, waterActivity, meltingTemperature, coolingRate, surfaceAreaBET, concentration)
  %% Read file 
  % This will give T.FreezeTemperature and T.Diameter_microns_
  T = readtable(inputFile,'Format', '%s%f%f', 'Delimiter', 'tab');
  
  % Calculate the Surface Area per droplet in cm^2
  % 1. calculate the concentration of ice nuclei in gm/cm^3
  % 2. calculate the volume of the droplets in cm^3
  % 3. Calculte the amount of Ice Nuclei mass in each droplet in gm
  % 4. Calculate the surface  per droplet using BET (m^2/gm) and Ice Nuclei
  % mass
  
  ConcgmPercm3 = concentration/100;
  volumeDropletcm3 = 4*pi/3*(T.Diameter_microns_/2*1e-4).^3; 
  massINPerDropletgm = volumeDropletcm3*ConcgmPercm3; 
  SurfaceAreaPerDropletcm2 = massINPerDropletgm*surfaceAreaBET*1e4;
  
  % Build the input matrix for Knopf calculation
  arrayLength = length(SurfaceAreaPerDropletcm2);
  A(:,1) = 1:arrayLength;
  A(:,2) = ones(arrayLength,1)*waterActivity; % water activity
  A(:,3) = 0.2675997*ones(arrayLength,1); % knopf correction
  A(:,4) = T.FreezeTemperature + 273.15; % freezing Temperature
  A(:,5) = ones(arrayLength,1) * meltingTemperature; % melting Temperature
  A(:,6) = ones(arrayLength,1) * coolingRate; % cooling rate
  A(:,7) = SurfaceAreaPerDropletcm2;
  
  KnopfHeadFunction(inputFile, A);
end