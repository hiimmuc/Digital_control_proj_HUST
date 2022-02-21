
%% Simulation edit section
Ts = 1/20; %sampling frequency 20Hz

%reference trajectory x(t), y(t)(sine)
T_cycle = 35; % time taken to traverse 1 cycle of sine trajectory (second)
t = 0:Ts:T_cycle; % point step of trajectory 
freq = 2*pi/T_cycle; %(NEW)
xRef = freq*t;
yRef = sin(freq*t);
dxRef = freq*ones([1 length(t)]);
dyRef = freq*cos(freq*t);
ddxRef = zeros([1 length(t)]);
ddyRef = -freq^2*sin(freq*t);
qRef = [xRef; yRef; atan2(dyRef, dxRef)]; %target location and orientation
uRef = [ddxRef; ddyRef]; %target acceleration

%INPUT state (q, w)
q = [0;0;0]; %INPUT X,Y LOCATION, PHI ANGLE (FROM CAMERA)
z1 = [q(1); dxRef(1)]; % initial x, dx state 
z2 = [q(2); dyRef(1)]; % initial y, dy state 
v = sqrt(z1(2)^2+z2(2)^2); % body velocity
dq = [dxRef(1); dyRef(1); atan2(dyRef(1), dxRef(1))]; %inital body velocity
w = J*dq; %initial INPUT wheel velocity (FROM ENCODER) (NEW)

%controller function
function [volt] = controller(w, q, xRef(k), yRef(k), dxRef(k), dyRef(k), ddxRef(k), ddyRef(k))
    %control parameter:
    l = 0.0765 ; % length robot center to wheel center
    d = 0.0695; % width robot center to middle of wheel
    r_w = 0.024; %radius of wheel
    
    %kinematic transformation matrix
    J = [1 -1 -(l+d); 1 1 -(l+d); 1 -1 l+d; 1 1 l+d]; %(body to motors)
    J_plus = inv((J')*J)*(J'); %pseudo-inverse kinematic, motors to body

    wMax = 0.2582; %maximum body angular velocity
    vMax = 0.377;  %maximum body velocity

    %state feedback controller
    desPoles = [-1-0.5i; -1+0.5i]; % pole placement
    %matrix of linearisation
    A = [0 1; 0 0];
    B = [0; 1];
    C = [1 0];
    K = place(A, B, desPoles); %pole placement
    %PID parameter
    Kp = [1;1;1;1];
    Ki = [10;10;10;10];
    Kd = [0;0;0;0];
    
    %update current state (NEW)
    dq = J_plus*w; %current body velocity 
    z1 = [q(1); dq(1)]; 
    z2 = [q(2); dq(2)];

    %reference state grouping
    zRef1 = [xRef(k); dxRef(k)];
    zRef2 = [yRef(k); dyRef(k)];

    %error and control
    ez1 = zRef1 - z1; %z1: INPUT x position and x speed
    ez2 = zRef2- z2; %z2: INPUT y position and y speed
    uu = [ddxRef(k); ddyRef(k)] + [K*ez1; K*ez2];
    D = norm([z1(1)-zRef1(1) z2(1)-zRef2(1)]);

    %Compute reference robot velocity
    F = [cos(q(3)), -v*sin(q(3)); sin(q(3)), v*cos(q(3))];
    vv = F\uu; %translational acceleration and angular velocity
    v = v + Ts*vv(1);
    u = [v; vv(2)];

    %constraint
    if abs(u(2))>wMax, u(2) = wMax*sign(u(2)); end %rotation speed constraint
    if abs(v)>vMax, v = vMax*sign(v); end %long speed constraint

    % Calculate OUTPUT voltage for 4 motors
    dqRef = [u(1)*cos(q(3)); u(1)*sin(q(3)); u(2)] ; %body velocity output
    vRef = (J*dqRef); % motor velocity    
    wRef = vRef.*60./(2*pi*r_w); %setpoint wheel rpm
    % PID controller (NEW)
    e = wRef - w;
    volt = Kp.*e + Ki.*e.*Ts + Kd.*(e-e_last)./Ts; %OUTPUT voltage 
    e_last = e;
end
