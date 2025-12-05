% Este codigo simula el funcionamiento de controladores On-Off, PID y MPC
% en la plataforma de pruebas de 1 grado de libertad.

clear
clc
close all

% Configuracion de simulacion
h = 0.1; % Paso de integracion para simulacion
Ts = 1; % Tiempo de muestreo para controlador MPC [s]
T = 9*3600; % Tiempo de simulacion [s]  (9 horas) x (3600s/hora)
pos0 = 73; % Posicion inicial al amanecer [°]. 73° --> acimutal, 25° --> elevacion
h0 = 8; % Hora de inicio de operacion (amanecer) [h]
eje = 0; % Seleccion de eje: 0 --> acimutal, 1 --> Elevacion 
controlador = 1; % Seleccion de controlador 1 --> On-Off, 2 --> PID, 3 --> MPC

% Parametros del sistema 
J = 5.1279e-4; % Momento de inercia
beta = 38.6673e-4; % Coeficiente de friccion viscosa
V_max = 12; % Voltaje maximo para el motor (positivo y negativo)
ka = 0.0382; % Constante de proporcionalidad par/corriente
R = 0.6454; % Resistencia de armadura
r = 6900; % Relacion de transmision

% Definicion del sistema en tiempo continuo
Ac = [0 1;0 -beta/J];
Bc = [0;(ka*V_max)/(J*R)];
Cc = (1/r)*(180/pi)*[1 0];
Dc = 0;
sysc = ss(Ac,Bc,Cc,Dc);

% Definición del sistema en tiempo discreto para simulacion
sys = c2d(sysc,h);
A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;

% Parametros del controlador On-Off
e_umbral = 0.3; % Error de umbral de activacion [°]
u_On = 0.3; % Magnitud de control en estado "On"

% Parametros del controlador PID
kp = 1.05; % Ganancia proporcional 
ki = 0.5; % Ganancia integral
kd = 0.1; % Ganancia derivativa

% Parámetros del controlador MPC
Nh = 10; % Horizonte de prediccion
Nc = 5; % Horizonte de control
rp = 1; % Coeficiente de ponderacion del error
qp = 3000; % Coeficiente de ponderacion del consumo
lb = -1*ones(Nc,1);
ub = ones(Nc,1);
U0 = zeros(Nc,1);
% Definicion del sistema en tiempo discreto para MPC
    sys_MPC = c2d(sysc,Ts);
    A_MPC = sys_MPC.A;
    B_MPC = sys_MPC.B;
    C_MPC = sys_MPC.C;
    D_MPC = sys_MPC.D;

t = 0:h:T; % Creacion del vector de tiempo

x0 = [pos0*r*(pi/180);0]; % Condicion inicial del vector de estado

% Incializacion del vector de estado
    x = zeros(size(A,1),size(t,2));
    x(:,1) = x0; 

y = zeros(size(t,2),1); % Inicializacion del vector de salida

%%

% Creacion del vector de referencia
[s1, s2] = trayectoria_solar(t,h0);
if eje == 0 
    y_ref = s1'*(180/pi);
end
if eje == 1
    y_ref = s2'*(180/pi);
end
% Inicializacion del vector de error
e = zeros(size(t,2)+1,1); 
% Considerar el corrimiento a la derecha del vector, de tal manera
% que e(k-1) (programción) = e[k] (tiempo discreto);

u = zeros(size(t,2),1); % Inicializacion del vector de control

I = 0; % Inicializacion del valor integral

t0 = 0; % Referencia de tiempo para peridos del controlador MPC

uk = 0; % Inicializacion de variable de control en tiempo de controlador
uk_ant = 0; % Variable para registro de la acción de control anterior
k0 = 0; % Referencia para tiempo de muestreo de controladores
K = Ts/h; % Relacion del tiempo de muestreo y el paso de integracion
k_MPC = 0; % Contador de instantes para horizonte de MPC

% Simulacion iterativa del sistema
for k = 1:1:length(t)-Nh*(Ts/h)

    y(k) = C*x(:,k); % Calculo de salida

    e(k+1) = y_ref(k) - y(k); % Calculo de error
    if(floor((k-k0)/K)==1)
        k0 = k;
    switch controlador
        case 1 % Controlador On-Off
            uk = 0;
            if e(k+1) >= e_umbral
                    uk = u_On;
            end
            if e(k+1) <= -e_umbral
                    uk = -u_On;
            end

        case 2 % Controlador PID
            I = I + Ts*((e(k)+e(k+1))/2); % Calculo de termino integral
            D = (e(k+1)-e(k))/Ts; % Calculo de termino derivativo
            uk = kp*e(k+1) + ki*I + kd*D; % Calculo de accion de control
            if k == 1; uk=0; end

        case 3 % Controlador MPC
            k_MPC = k_MPC+1;
            if(1 <= k_MPC && k_MPC <= Nh)
                if(k_MPC == 1)
                    fun = @(U) funcion_costo(A_MPC,B_MPC,C_MPC,y_ref(k+1:Ts/h:k+(Ts/h)*Nh),rp,qp,Nh,Nc,x(:,k),U);
                    U = fmincon(fun,U0,[],[],[],[],lb,ub);
                end

                if(1 <= k_MPC && k_MPC <= Nc)
                    uk = U(k_MPC);
                else
                    uk = 0;
                end
            end
            if(k_MPC == Nh)
                k_MPC = 0;
            end
    end
    else
    uk = uk_ant;
    end

    u(k) = uk;
    x(:,k+1) = A*x(:,k) + B*u(k); % Calculo del vector de estado siguiente
    
    uk_ant = uk;
end
%%
y = y(1:end-Nh*(Ts/h));
y_ref = y_ref(1:end-Nh*(Ts/h));
t = t(1:end-Nh*(Ts/h));
u = u(1:end-Nh*(Ts/h));
e = e(1:end-Nh*(Ts/h));
%%
figure
plot(t',y,"Color",'b','LineStyle','-','LineWidth',2)
hold on
plot(t',y_ref,"Color",'r','LineStyle','--','LineWidth',2)
grid on
legend({'Posición de junta','Referencia'})
ylabel('$q(t)$ [$^\circ$]','Interpreter','latex','FontSize',12)
xlabel('Tiempo [$h$]', 'Interpreter','latex','FontSize',12)
etiquetas_x = linspace(0,t(end),10);
xticks(etiquetas_x)
xticklabels({'8','9','10','11','12','13','14','15','16','17'})
xlim([0 t(end)])
ylim([-100 80])

figure
plot(t',y_ref-y,'Color','b')
grid on
ylabel('Error de seguimiento $[^{\circ}]$','Interpreter','latex','FontSize',12)
xlabel('Tiempo [$h$]', 'Interpreter','latex','FontSize',12)
etiquetas_x = linspace(0,t(end),10);
xticks(etiquetas_x)
xticklabels({'8','9','10','11','12','13','14','15','16','17'})
xlim([0 t(end)])
ylim([-0.4 0.4])

figure
plot(t',u,"Color",'b','LineStyle','-')
grid on
ylabel('$u(t)$','Interpreter','latex','FontSize',12)
xlabel('Tiempo [$h$]', 'Interpreter','latex','FontSize',12)
etiquetas_x = linspace(0,t(end),10);
xticks(etiquetas_x)
xticklabels({'8','9','10','11','12','13','14','15','16','17'})
xlim([0 t(end)])
ylim([-1.1 1.1])

MSE = (e'*e)/length(e)
MAE = sum(abs(e))/length(e)
Energia = sum(abs(u))

function J = funcion_costo(A_MPC,B_MPC,C_MPC,Y_ref,r,q,Nh,Nc,xk,U)
Y = zeros(Nh,1);
for k=1:1:Nh
    if k<=Nc
    xk1 = A_MPC*xk + B_MPC*U(k);
    else
    xk1 = A_MPC*xk;
    end
    Y(k) = C_MPC*xk1;
    xk = xk1;
end
J = r*(Y_ref-Y)'*(Y_ref-Y) + q*(U')*U;
end 