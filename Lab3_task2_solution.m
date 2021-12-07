clc
clear all
close all
%% Spindulio tipo bazinė funkcija su c ir r parinkimu

%% 1. Duomenys
x = 0.1:1/22:1;
d = (1+ 0.6*sin(2*pi*x/0.7) + 0.3*sin(2*pi*x))/2;
plot(x,d,'ro'); hold on

%% 2. Randamos c vertės ieškant funkcijos maksimumų - skirtumui keičiantis iš teigiamo į neigiamą, funkcija yra maksimume (išvestinės)
k = 0;
for j = 1:length(x)-1
    diff(j) = d(j+1)-d(j);
    
    if j >1
        
        if diff(j-1) >0 & diff(j) <0
            k = k+1;
            maxima(k) = j;
        end
 
    end
    
end

c1 = x(maxima(1));
c2 = x(maxima(2));

%% Parenkamos mažiausios r vertės pradžiai, nustatoma maksimalios r vertės

r_min = 0.05; %r1_min = r2_min
r_max = 1; %r1_max = r2_max

%% 3. Atsitiktiniai svoriai w ir bias b

w1 = rand(1);
w2 = rand(1);
b = rand(1);
n = 0.1;
rCount = 0;
%% 4. Tinklo atsakas kintant r reikšmei 
for r_change = r_min:(r_max-r_min)/20:r_max
    rCount = rCount+1;
    r1 = r_change;
    r2 = r_change;
            
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

    x = 0.1:1/22:1;

    for i = 1:length(x)

        Phi1 = exp(-(x(i)-c1)^2/(2*r1^2));
        Phi2 = exp(-(x(i)-c2)^2/(2*r2^2));

        y(i) = Phi1*w1 + Phi2*w2 + b;
        E(i,rCount) = d(i) - y(i);
    end
end

%% Randama r reiksšmė kuriai esant buvo mažiausia klaida
avgErrors = mean(E);
[M,k] = min(avgErrors);

%% Apmokoma iš naujo su optimalia r reikšme
for iteration = 1:100
    for i = 1:20
        
        r1 = ((r_max-r_min)/20)*k;
        r2 = ((r_max-r_min)/20)*k;
         
        % spindulio funkcijos
        Phi1 = exp(-(x(i)-c1)^2/(2*r1^2));
        Phi2 = exp(-(x(i)-c2)^2/(2*r2^2));

        y = Phi1*w1 + Phi2*w2 + b;

        %% Klaida

        e = d(i)-y;

        %% Atnaujinami svoriai

        w1 = w1 + n*e*Phi1; 
        w2 = w2 + n*e*Phi2;
        b = b + n*e;

    end
end

%% Tikrinimas

x = 0.1:1/100:1;

for i = 1:length(x)

    Phi1 = exp(-(x(i)-c1)^2/(2*r1^2));
    Phi2 = exp(-(x(i)-c2)^2/(2*r2^2));

    y(i) = Phi1*w1 + Phi2*w2 + b;

end

plot(x,y,'g.')