function [s1, s2] = trayectoria_solar(t,h0)

t = t./3600;
t = t+h0;
% Vector de tiempo como entrada (en horas)
s1 = zeros(1,length(t));
s2 = zeros(1,length(t));

for i=1:1:length(t)
[s1(i), s2(i)] = posicion_solar(t(i));
end

function [z, h] = posicion_solar(t)
    n = 93; % Numero de día;
    lambda = 40.48*(pi/180.0); % Latitud
    L = -3.68*(pi/180.0); % Longitud
    sn = 12; % Medio día solar
    delta=(23.45*(pi/180))*sin(2*pi*((284+n)/365)); % Declinación solar
    tau=(sn-t)*(15*(pi/180)); % Angulo horario 
    sinh1 =(cos(lambda)*cos(delta)*cos(tau))+(sin(lambda)*sin(delta));
    h=asin(sinh1); % Elevacion
    cosz1=(sin(h)*sin(lambda)-sin(delta))/(cos(h)*cos(lambda));
    z=acos(cosz1); % Acimutal
    if t>sn
        z=-z;
    end
end

end