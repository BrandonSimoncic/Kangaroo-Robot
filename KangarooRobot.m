%By Brandon Simoncic and Eddie Doemer

%==========================================================================
%Control System Design and Implementation
%==========================================================================

syms xt xb xl

%Set up our variables
Ib = (1/3)*3.45*.10^2
It = ((.1/3)+0.52)*.4318^2
Il = (.25*.2921^2)/3
g = 3.711
Pb = 3.45*.10*g
Pt = (.1 + .52)*.213*g
Pl = Il*g

Ts = .01

%Continuous state space using our small angle approximations.
%Even though small angle approximations will not work with the expected 20
%degree sweep, we assumed we could design a controller that would still
%work with this
%Because we approximate the potential energy of the leg to be minimal under
%small angles, it's multiplied by .1
A = [0 1 0 0 0 0; -Pb/Ib 0 0 0 0 0; 0 0 0 1 0 0; 0 Pt/It 0 0 0 0; 0 0 0 0 0 1; 0 0 0 0 (.1)*Pl/Il 0]
B = [0 0; 1/Ib 1/Ib; 0 0; 1/Il 0; 0 0; 0 1/It]
C = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];
D = 0;

ro = .001

%Q and R matricies are chosen to hold body still
%The Q matricie puts a large penelty on the body not being at it's setpoint
%It is penalized the least for the motion of the leg and tail, since they
%are intended to be moved around.
Q = [100000000 0 0 0 0 0; 0 10000 0 0 0 0; 0 0 1 0 0 0; 0 0 0 .1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1];

%The R matrix penalizes the controller for trying to move the leg, since
%the leg should only be moved by our predetermined input.
%since ro is very small, the penalty for trying to move the tail is small
%We have alot of energy available for the tail, but almost none
%available for the legs
R = ro*[1 0; 0 100000];





%Create a state space system and use C2D with zoh to discretize it
sys = ss(A,B,C,D)
sysz = c2d(sys, Ts, 'zoh');

%Arbitrary stable poles were placed to test the system, but commented out
%to use LQR
%k = place (sysz.A, sysz.B, [.9 .88 .87 .86 .85 .84]);

%Get our gain matrix from our digital system and our Q and R matricies.
%This k matrix should be very close to zero for the bottom row.
%We assume this makes it similar to a feedfoward input
k = dlqr(sysz.A, sysz.B, Q, R)

%Get our closed loop feedback with controller applied using our A and B
%matricies along with the gain matrix found with dlqr()
Gz = sysz.A - sysz.B*k;


%Create some labels for all of our states, that way it is easier to view
states = {'Body Angle', 'Body Omega', 'Leg Angle', 'Leg Omega', 'Tail Angle', 'Tail Omega'}
inputs = {'Tl', 'Tt'}
outputs = {'Body Orientation', 'Leg Orientation', 'Tail Orientation'}
%outputs = {'Body Orientation', 'Tail Orientation'}

%apply these labels
sys = ss(Gz, sysz.B, sysz.C, sysz.D, .01, 'statename', states, 'inputname', inputs, 'outputname', outputs)

%Setup a time vector and force vector for the leg
t = 0:.01:20;
%We wanted to have it stabalize a leg with a frequency of 1/0.89 hz


Legmove = 22*cos(.89*pi()*2*t);
Txx = Legmove';

%It would be better to simulate this with an external force being applied
%to the body, but matlab didn't like to do that with state space, thought
%it would be simulated in simulink
%This input vector is passed to it. 
u = [zeros(length(t), 1) Txx];


[x1 t1 xs] = lsim(sys, u, t);

%Create the plots 
subplot(2,1,1)
plot(t1, x1)
legend({'y = Body', 'y = Leg', 'y = Tail'})
title('Cosine T-leg Input Response')
xlabel('Time')
ylabel('Angle in Degrees')

%unco = length(A) - rank(ctrb(A,B))


%==========================================================================
%Simulation Start
%==========================================================================

% define the x- and y-data for the original line we would like to rotate
xT = [0 -1 -2];
yT = [0 0 0];
xB = [0 1 1 0 0];
yB = [0 0 1 1 0];
xL = [0 -0.5 0];
yL = [0 -1 -2];

% create a matrix of these points, which will be useful in future calculations
vT = [xT;yT];
vB = [xB;yB];
vL = [xL;yL];

% choose a points which will be the center of rotation
xT_center = xT(1);
yT_center = yT(1);
centerT = repmat([xT_center; yT_center], 1, length(xT));
xB_center = xB(1);
yB_center = yB(1);
centerB = repmat([xB_center; yB_center], 1, length(xB));
xL_center = xL(1);
yL_center = yL(1);
centerL = repmat([xL_center; yL_center], 1, length(xL));
subplot(2,1,2)
title('Simulation')
legend('Body', 'Tail', 'Leg')


for i = 1:16:max(size(t))
    
vect = [x1(i,1),x1(i,2),x1(i,3)];

thetaBody = vect(1)*pi()/180;
thetaLeg = vect(2)*pi()/180;
thetaTail = vect(3)*pi()/180;


RT = [cos(thetaTail) -sin(thetaTail); sin(thetaTail) cos(thetaTail)];
vT2 = RT*(vT - centerT) + centerT;
% pick out the vectors of rotated x- and y-data
xT = vT2(1,:);
yT = vT2(2,:);

RB = [cos(thetaBody) -sin(thetaBody); sin(thetaBody) cos(thetaBody)];
vB2 = RB*(vB - centerB) + centerB;
% pick out the vectors of rotated x- and y-data
xB = vB2(1,:);
yB = vB2(2,:);

RL = [cos(thetaLeg) -sin(thetaLeg); sin(thetaLeg) cos(thetaLeg)];
vL2 = RL*(vL - centerL) + centerL;
% pick out the vectors of rotated x- and y-data
xL = vL2(1,:);
yL = vL2(2,:);

% Iterate through the plot

%%%%%Annimated line (VERY SLOW)
% addpoints(anB, t1(i), x1(i, 1));
% drawnow
% addpoints(anL, t1(i), x1(i, 2));
% drawnow
% addpoints(anT, t1(i), x1(i, 3));
% drawnow
% legend({'y = Body', 'y = Leg', 'y = Tail'})


plot(xT, yT, 'b-', xB,yB, 'g-',xL,yL, 'r-', xT_center, yT_center, 'bo');
axis equal
axis([-3 3 -3 3]);
pause(.1)  % Shows results at each time interval

end

%By Brandon Simoncic and Eddie Doemer