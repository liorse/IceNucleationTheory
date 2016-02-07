% file name of freezing data from freezing detection software
% waterActivity (in numbers from 0 to 1)
% meltingTemperature ( assumed a single number in which the droplets froze
% coolingRate (in Kelvin/min )
% surfaceArea from BET measurement ( meter^2/gm )
% concentration (% weight)
%CalculateJhetNs('Kaolinite IN- RAW -20d00-60d00--40d00-0d00-0d01-0d00.txt', ...
%                 1, 273.15, 1, 11.8, 1);
%CalculateJhetNs('K-Feldspar- RAW -20d00-60d00--60d00-0d00-0d01-0d00.txt', ...
 %                1, 273.15, 1, 1.86, 1);

%CalculateJhetNs('ATD 1 WT%- RAW -20d00-60d00--60d00-0d00-0d01-0d00.txt', ...
%                 1, 273.15, 1,5, 1);

CalculateJhetNs('0.1 g ncc heterogeneous- RAW -20d00-60d00--36d00-0d00-0d01-0d00.txt', ...
                 1, 273.15, 10, 8, 0.1);
%CalculateJhetNs('0.2 g NCC heterogeneous- RAW -20d00-60d00--36d00-0d00-0d01-0d00.txt', ...
%                 1, 273.15, 10, 8, 0.1);
             
             