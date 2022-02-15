
%% Simulation edit section
Ts = 1/20; %sampling frequency 60Hz

t = 0:Ts:30; % simulation time step
l = 0.0765 ; % length robot center to wheel center
d = 0.0695; % width robot center to middle of wheel
r_w = 0.024; %radius of wheel

wMax = pi/2; %maximum body angular velocity
vMax = 0.377;  %maximum body velocity
noise = 0.00; % Gaussian noise of location tracking (e.g. 0.001)


%state feedback controller
desPoles = [-1-0.5i; -1+0.5i]; % pole placement
%desPoles = [1+0.5i; 1-0.5i];
%%

%kinematic transformation matrix
J = [1 -1 -(l+d); 1 1 -(l+d); 1 -1 l+d; 1 1 l+d]; %(body to motors)
J_plus = inv((J')*J)*(J'); %pseudo-inverse kinematic, motors to body

%reference trajectory x(t), y(t)(circle)
freq = 2*pi/30;
xRef = freq*t;
yRef = sin(freq*t);
dxRef = freq*ones([1 length(t)]);
dyRef = freq*cos(freq*t);
ddxRef = zeros([1 length(t)]);
ddyRef = -freq^2*sin(freq*t);
qRef = [xRef; yRef; atan2(dyRef, dxRef)]; %target location and orientation
uRef = [ddxRef; ddyRef]; %target acceleration

%recorded state
q = [0;0; 0]; %initial locaiton
z1 = [q(1); dxRef(1)]; % x, dx state
z2 = [q(2); dyRef(1)]; % y, dy state
v = sqrt(z1(2)^2+z2(2)^2); % body velocity

%matrix of linearisation
A = [0 1; 0 0];
B = [0; 1];
C = [1 0];
K = place(A, B, desPoles); %pole placement

%plot variables
v_r = zeros([4 length(t)]); %wheel velocity
w_m = zeros([4 length(t)]); %angular velocity wheel
q_r = zeros([3 length(t)]); %position
w_r = zeros([1 length(t)]); %angular velocity
dq_r = zeros([3 length(t)]); %body motion (velocity x,y, angular velocity)
D_r = zeros([1 length(t)]); %distance from body to destination point

for k = 1:length(t)
    %reference state 
    zRef1 = [xRef(k); dxRef(k)];
    zRef2 = [yRef(k); dyRef(k)];
    
    %error and control
    ez1 = zRef1 - z1; %z1:  x position and x speed
    ez2 = zRef2- z2; %z2:  y position and y speed
    uu = [ddxRef(k); ddyRef(k)] + [K*ez1; K*ez2];
    D = norm([z1(1)-zRef1(1) z2(1)-zRef2(1)]);
    %Compute inputs to the robot
    F = [cos(q(3)), -v*sin(q(3)); sin(q(3)), v*cos(q(3))];
    vv = F\uu; %translational acceleration and angular velocity
    v = v + Ts*vv(1);
    u = [v; vv(2)];

    %constraint
    if abs(u(2))>wMax, u(2) = wMax*sign(u(2)); end
    if abs(v)>vMax, v = vMax*sign(v); end
    
    % input to robot
    dq = [u(1)*cos(q(3)); u(1)*sin(q(3)); u(2)] ; %body velocity output
    v_m = (J*dq); %setpoint motor velocity    
    w_m(:,k) = v_m./(2*pi*r_w);

    %simulation
    q = q + Ts*dq + randn(3,1)*noise; % Euler integration, coordinate received from camera
    q(3) = wrapToPi(q(3));
    v_r(:,k) = v_m;
    q_r(:,k) = q;
    dq_r(:,k) = dq;
    D_r(:,k) = D;
    %update current state
    z1 = [q(1); u(1)*cos(q(3))];
    z2 = [q(2); u(1)*sin(q(3))];
    
end
tiledlayout(2,2)
nexttile;
plot(q_r(1,:), q_r(2,:));
hold on
plot(xRef, yRef, '--','Color','red'); 
grid on;
title("xy plane charting");
legend('robot path', 'reference track');

nexttile;
plot(t, v_r);
grid on;
title("Individual wheel input velocity");
legend('wheel 1', 'wheel 2', 'wheel 3', 'wheel 4');

nexttile([1 2]);
plot(t, D_r);
grid on;
title("Error distance from body to reference point");