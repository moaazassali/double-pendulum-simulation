% DOUBLE PENDULUM SIMULATION SOLVED USING EULER AND RK4 METHOD
% Moaaz Assali (ma5679@nyu.edu)
% Brian Chesney Quartey (bcq205@nyu.edu)
% ENGR-UH 2017 Numerical Methods with Prof. Saif Eddin Jabari (sej7@nyu.edu)
% New York University - Abu Dhabi

%-------READ ME------------------------------------------------------------
% Only parameters to change to begin simulation: 
% time, timeStep, th1(1), th2(1), w1(1), w2(1)
% Choose whether the simulation should run with Euler's method or 4th-order
% Runge-Kutta method by commenting out the appropriate one
%--------------------------------------------------------------------------

clear; clc;

%Declaring global constants to be used in all functions
global L1; L1 = 1;      
global L2; L2 = 1;
global m1; m1 = 1;
global m2; m2 = 1;
global g; g = 9.81;

%Runtime of the simulation and step size
time = 10;
global timeStep; timeStep = 0.01;
iterations = time/timeStep+1;

%Initializing the vectors of every quantity
t = zeros(iterations,1);
th1 = zeros(iterations,1);
th2 = zeros(iterations,1);
w1 = zeros(iterations,1);
w2 = zeros(iterations,1);
x1 = zeros(iterations,1);
x2 = zeros(iterations,1);
y1 = zeros(iterations,1);
y2 = zeros(iterations,1);
E = zeros(iterations,1);

%Initial conditions of the ODE system
t(1) = 0; %setting element 1 to be t = 0 seconds
th1(1) = 2*pi/3; %initial angle of upper bob
th2(1) = pi; %initial angle of lower bob
w1(1) =  0; %initial angular velocity of upper bob
w2(1) = 0; %initial angular velocity of lower bob

%Calculation of initial quantities
x1(1) = L1*sin(th1(1));
y1(1) = -L1*cos(th1(1));
x2(1) = x1(1) + L2*sin(th2(1));
y2(1) = y1(1) - L2*cos(th2(1));
dx1 = w1(1)*L1*cos(th1(1));
dy1 = w1(1)*L1*sin(th1(1));
dx2 = dx1 + w2(1)*L2*cos(th2(1));
dy2 = dy1 + w2(1)*L2*sin(th2(1));
K = (m1*(dx1^2 + dy1^2))/2 + (m2*(dx2^2 + dy2^2))/2;
U = m1*g*(y1(1)+2) + m2*g*(y2(1)+2);
E(1) = K + U;

global h; h = timeStep;
for i = 2:iterations
    %Euler or RK4 method can be used to solve the system of ODEs (comment
    %out one method only)
    %Euler method
%     [th1(i), th2(i), w1(i), w2(i)] = euler(th1(i-1),th2(i-1),w1(i-1),w2(i-1));
    
    %4th-order runge-Kutta (RK4)
    [th1(i), th2(i), w1(i), w2(i)] = rk4(th1(i-1),th2(i-1),w1(i-1),w2(i-1));
    
    %Calculation of rectangular coordinates of each ball
    x1(i) = L1*sin(th1(i));
    y1(i) = -L1*cos(th1(i));
    x2(i) = x1(i) + L2*sin(th2(i));
    y2(i) = y1(i) - L2*cos(th2(i));
    
    
    t(i) = t(i-1) + h; %increments time by 1 timeStep
    
    %Calculation of energy of the system
    dx1 = w1(i)*L1*cos(th1(i));
    dy1 = w1(i)*L1*sin(th1(i));
    dx2 = dx1 + w2(i)*L2*cos(th2(i));
    dy2 = dy1 + w2(i)*L2*sin(th2(i));
    K = (m1*(dx1^2 + dy1^2))/2 + (m2*(dx2^2 + dy2^2))/2;
    U = m1*g*(y1(i)+2) + m2*g*(y2(i)+2);
    E(i) = K + U;
    
    % Plot of both balls in x-y coordinate space (skips some frames to keep
    % the simulation as close to real-time as possible (doesn't impact any
    % values or plots)
    if (rem(i,0.02/timeStep) == 0) || (i==iterations)
        plot([0, x1(i), x2(i)], [0, y1(i), y2(i)], '-o');
        axis([-2 2 -2 2]);
        xlabel("x-position (m)");
        ylabel("y-positoin (m)");
        title(['Plot of double pendulum system at time = ', num2str(t(i), '% 1.2f'), 's']); %dynamically changes time displayed in the title
        hold on
        plot(x2(1:i), y2(1:i), 'g');
        hold off
        pause(0.00001);
    end
end

figure
plot(t,th1(:),'b',t,th2(:),'r');    %Plot of angles of both pendulum balls
xlabel("Time (s)");
ylabel("Angle (rad)");
legend('\Theta_1','\Theta_2')
title("Plot of angles of each pendulum bob over time");
figure
plot(t,E(:),'b');   %Plot of energy of double pendulum system over time
xlabel("Time (s)");
ylabel("Total Energy (J)");
title("Plot of total energy of the double pendulum system over time");

%First 2nd order ODE but written as 1st order ODE
function [y] = dw1(th1,th2,w1,w2)
    global L1;
    global L2;
    global m1;
    global m2;
    global g;
    
    w1Num_1 = -g*(2*m1+m2)*sin(th1) - m2*g*sin(th1-2*th2);
    w1Num_2 = -2*sin(th1-th2)*m2*((w2.^2)*L2 + (w1.^2)*L1*cos(th1-th2));
    w1Num = w1Num_1 + w1Num_2;
    wDen = L1*(2*m1+m2-m2*cos(2*th1-2*th2));
    y = w1Num / wDen;
end

%Second 2nd order ODE but written as 1st order ODE
function [y] = dw2(th1,th2,w1,w2)
    global L1;
    global L2;
    global m1;
    global m2;
    global g;
    
    w2Num = 2*sin(th1-th2) * ((w1.^2)*L1*(m1+m2) + g*(m1+m2)*cos(th1) + (w2.^2)*L2*m2*cos(th1-th2));
    wDen = L2*(2*m1+m2-m2*cos(2*th1-2*th2));
    y = w2Num / wDen;
end

function [th1new,th2new,w1new,w2new] = euler(th1,th2,w1,w2)
    global h;
    
    th1new = th1 + h*w1;
    th2new = th2 + h*w2;
    w1new = w1 + h*dw1(th1,th2,w1,w2);
    w2new = w2 + h*dw2(th1,th2,w1,w2);
end

function [th1new,th2new,w1new,w2new] = rk4(th1,th2,w1,w2)
    global h;

    K1th1 = w1;
    K1th2 = w2;
    K1w1 = dw1(th1,th2,w1,w2);
    K1w2 = dw2(th1,th2,w1,w2);
    
    K2th1 = w1 + (h*K1w1)/2;
    K2th2 = w2 + (h*K1w2)/2;
    K2w1 = dw1(th1+(h*K1th1)/2, th2+(h*K1th2)/2, w1+(h*K1w1)/2, w2+(h*K1w2)/2);
    K2w2 = dw2(th1+(h*K1th1)/2, th2+(h*K1th2)/2, w1+(h*K1w1)/2, w2+(h*K1w2)/2);
    
    K3th1 = w1 + h*K2w1/2;
    K3th2 = w2 + h*K2w2/2;
    K3w1 = dw1(th1+h*K2th1/2, th2+h*K2th2/2, w1+h*K2w1/2, w2+h*K2w2/2);
    K3w2 = dw2(th1+h*K2th1/2, th2+h*K2th2/2, w1+h*K2w1/2, w2+h*K2w2/2);
    
    K4th1 = w1 + h*K3w1;
    K4th2 = w2 + h*K3w2;
    K4w1 = dw1(th1+h*K3th1, th2+h*K3th2, w1+h*K3w1, w2+h*K3w2);
    K4w2 = dw2(th1+h*K3th1, th2+h*K3th2, w1+h*K3w1, w2+h*K3w2);
    
    th1new = th1 + (h/6)*(K1th1 + 2*K2th1 + 2*K3th1 + K4th1);
    th2new = th2 + (h/6)*(K1th2 + 2*K2th2 + 2*K3th2 + K4th2);
    w1new = w1 + (h/6)*(K1w1 + 2*K2w1 + 2*K3w1 + K4w1);
    w2new = w2 + (h/6)*(K1w2 + 2*K2w2 + 2*K3w2 + K4w2);
end