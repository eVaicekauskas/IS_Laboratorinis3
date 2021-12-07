clc
clear all
close all
%% Spindulio tipo bazinė funkcija 

%% 1. Duomenys
x = 0.1:1/22:1;
d = (1+ 0.6*sin(2*pi*x/0.7) + 0.3*sin(2*pi*x))/2;
plot(x,d,'ro'); hold on

%% 2. Rankiniu būdu pasirenkamos c ir r vertės

c1 = 0.19;
c2 = 0.91;
r1 = 0.35/2;
r2 = 0.4/2;

%% 3. Atsitiktiniai svoriai w ir bias b

w1 = rand(1);
w2 = rand(1);
b = rand(1);
n = 0.1;
%% 4. Tinklo atsakas

for iteration = 1:100
    for i = 1:20
        
        % spindulio funkcijos
        Phi1 = exp(-(x(i)-c1)^2/(2*r1^2));
        Phi2 = exp(-(x(i)-c2)^2/(2*r2^2));

        y = Phi1*w1 + Phi2*w2 + b;

        %% 5. Klaida

        e = d(i)-y;

        %% 6. Atnaujinami svoriai

        w1 = w1 + n*e*Phi1; 
        w2 = w2 + n*e*Phi2;
        b = b + n*e;

    end
end

%% 7. Tikrinimas

x = 0.1:1/100:1;

for i = 1:length(x)
    
    Phi1 = exp(-(x(i)-c1)^2/(2*r1^2));
    Phi2 = exp(-(x(i)-c2)^2/(2*r2^2));

    y(i) = Phi1*w1 + Phi2*w2 + b;
    
end

plot(x,y,'g.')