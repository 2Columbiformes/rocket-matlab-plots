%
inertiaM4 = sym([1 0 0 0;
                0 0.023 0 0;
                0 0 0.023 0;
                0 0 0 0.005]);
inertiaM = sym([0.023 0 0;
                0 0.023 0;
                0 0 0.005]);
mass = 1;
g = 9.81;
%Isp = c/g;
kP = 1;
kD = 0.5;
P2 = 3;
D2 = 50;
Leverarm = 0.5; %(m leverarm of nozzle)
Thrust = 20; %(N force of thruster)



function ANS = sym(M)
    ANS = 0.5*(M + transpose(M));
end
function ANS = asym(M)
    ANS = 0.5*(M - transpose(M));
end
function ANS = KE(vec,inertiaM)
    ANS = 0.5*dot(vec,inertiaM*vec);
end
function ANS = lnq(q)
    q2 = [q(2,:);q(3,:);q(4,:)];
    q0 = [q(1,:)];
    ANS = [1/2*log(dot(q,q));normalize(q2).*atan2(sqrt(dot(q2,q2)),q0)];
end
function ANS = expq(q)
    q2 = [q(2,:);q(3,:);q(4,:)];
    q0 = [q(1,:)];
    M = sqrt(dot(q2,q2));
    ANS = exp(q0).*[cos(M);normalize(q2).*sin(M)];
end
function ANS = pow(q,p)
    ANS = expq(p*lnq(q));
end
function ANS = magn(q)
    ANS = zeros(1,size(q,2));
    for i = 1:size(q,2)
        ANS(1,i) = sqrt(dot(q(:,i),q(:,i)));
    end
end

function ANS = Rotate(vec,quat)
    quatv = [quat(2,:);quat(3,:);quat(4,:)];
    quatw = [quat(1,:)];
    ANS = zeros(3,size(quat, 2));
    if size(vec,2) == 1
        vec = repmat(vec, 1, size(quat, 2));
    end
    for i = 1:size(quat, 2)
        ANS(:,i) = vec(:,i) + 2*cross(quatv(:,i),quatw(:,i)*vec(:,i) + cross(quatv(:,i),vec(:,i)))./(dot(quat(:,i),quat(:,i)));
    end
 
end
function ANS = torquetoq(vec)
    ANS = zeros(4,size(vec, 2));
    for i = 1:size(vec, 2)
        ANS(:,i) = 1/2*[sqrt(1-magn(vec(:,i)))+sqrt(1+magn(vec(:,i)));(sqrt(1+magn(vec(:,i)))-sqrt(1-magn(vec(:,i))))*normalize(vec(:,i))];
    end
end
function ANS = vtoq(vec1,vec2)
    ANS = zeros(4,size(vec1, 2));
    for i = 1:(size(vec1, 2)+size(vec2, 2))/2
        ANS(:,i) = normalize(pow([dot(vec1(:,i),vec2(:,i));cross(vec1(:,i),vec2(:,i))],1/2));
    end
end
function ANS = angfrom0(quat)
    ANS = zeros(1,size(quat, 2));
    for i = 1:size(quat,2)
        Q(:,i) = vtoq([0;0;1],Rotate([0;0;1],quat(:,i)));
        ANS(1,i) = 2*acos(Q(1,i));
    end
end
function ANS = clamp(quat)
    A = angfrom0(quat);
    if A > 5*pi/180
       ANS = pow(quat,pi/36/A);
    else
       ANS = quat;
    end
end
function ANS = clampv(vec)
    ANS = zeros(3,size(vec, 2));
    V = vec;
    C = magn(cross([0;0;1],Rotate([0;0;1],pow([0;1;0;0],5/180))));
    for i = 1:size(vec,2)
        V(3,i) = 0;
        if magn(V(:,i)) > C
            ANS(:,i) = C*normalize(V(:,i));
        else
           ANS(:,i) = V(:,i);
        end
    end
end


function ANS = Rotate2(vec,quat)
    vec1 = [vec(1,:);vec(3,:);vec(5,:)];
    vec2 = [vec(2,:);vec(4,:);vec(6,:)];
    ANS1 = zeros(3,size(vec,2));
    ANS2 = zeros(3,size(vec,2));
    for i = 1:size(vec,2)
        ANS1(:,i) = Rotate(vec1(:,i),quat);
        ANS2(:,i) = Rotate(vec2(:,i),quat);
    end
    ANS = [ANS1(1,:);ANS2(1,:);ANS1(2,:);ANS2(2,:);ANS1(3,:);ANS2(3,:)];
end
function ANS = normalize(q)
    ANS = zeros(size(q,1),size(q,2));
    for i = 1:size(q,2)
        D = dot(q(:,i),q(:,i));
        if D > 10^(-8)
            ANS(:,i) = q(:,i)/sqrt(D); 
        else
            if size(q,1) == 4
                ANS(:,i) = [1;0;0;0];
            else
                ANS(:,i) = zeros(size(q,1),1);
            end

        end
    end
end
function ANS = d2(q)
    ANS=dot(q,q);
end
function ANS = Biv(q)
    ANS = [q(1,:),q(2,:),q(3,:)];
end

function ANS = qinv(q)
    ANS=[q(1,:);-q(2,:);-q(3,:);-q(4,:)]./dot(q,q);
end
function ANS = conj(q)
    ANS=[q(1,:);-q(2,:);-q(3,:);-q(4,:)];
end
% takes in column vector, spits out matrix with time in 3d dim.
function ANS = qmat(q)
    ANS = zeros(4, 4, size(q, 2));
    for i = 1:size(q, 2)
        ANS(:,:,i) = [ q(1,i) -q(2,i) -q(3,i) -q(4,i) 
                       q(2,i)  q(1,i) -q(4,i)  q(3,i)
                       q(3,i)  q(4,i)  q(1,i) -q(2,i) 
                       q(4,i) -q(3,i)  q(2,i)  q(1,i) ];
    end
end
%takes in column vector q, takes in column vector qdot spits out row vector
function ANS = Pconj(q,qdot,inertiaM4)
    Qmat = qmat(q);
    ANS = zeros(4,size(q, 2));
    for i = 1:size(q, 2)
        ANS(:,i) = 4*Qmat(:,:,i)*inertiaM4*(Qmat(:,:,i)')*qdot(:,i)/(dot(q(:,i),q(:,i))^2);
    end
end
function ANS = multq(q,p)
    ANS = zeros(4, size(q, 2));
    Qmat = qmat(q);
    for i = 1:size(q, 2)
        ANS(:,i) = Qmat(:,:,i)*p(:,i);
    end
end
%{
function ANS = QC(q,p)
    Dcontrib = (pow(p,0.1*sqrt(dot(p,p))));
    qid = [1;0;0;0];
    %ANS = normalize(multq(Dcontrib,pow(q,0.05)));
    ANS = multq(expq(0.5*multq(conj(q),p)),pow(q,0.2));
end
%}
function ANS = QC(q,p,kP,kD)
    qid = [1;0;0;0];
    inertiaM4 = sym([1 0 0 0;
                0 0.023 0 0;
                0 0 0.023 0;
                0 0 0 0.005]);
    %ANS = normalize(multq(Dcontrib,pow(q,0.05)));
    ANS = multq(expq(kD*inertiaM4\multq(conj(q),p)),pow(q,kP));
end
function ANS = TOR(q,Amn,kP,kD,Leverarm,Thrust,inertiaM)
    ANS = zeros(3, size(q, 2));
    for i = 1:size(q,2)    %real(sqrt(1-q(4,:)^2))
        ANS(:,i) = -Leverarm*Thrust*clampv(((2*kP-0.5*d2(inertiaM\Amn))*inertiaM*(([q(2,i);q(3,i);q(4,i)])/(q(1,i)) + kD*Amn(:,i) - 0*cross(Amn(:,i),inertiaM\Amn(:,i)))));
        %ANS(:,i) = -Leverarm*Thrust*clampv(((2*kP-0.5*d2(inertiaM\Amn))*inertiaM*(([q(2,i);q(3,i);q(4,i)])*real(sqrt(1-q(4,:)^2))/(q(1,i)) + kD*Amn(:,i) - 0*cross(Amn(:,i),inertiaM\Amn(:,i)))));
        %ANS(:,i) = -Leverarm*Thrust*clampv((2*kP)*inertiaM*([q(2,i);q(3,i);q(4,i)]) + kD*Amn(:,i));
        
        %ANS(:,i) = -Leverarm*Thrust*clampv((Pro*inertiaM*2*normalize([q(2,i);q(3,i);q(4,i)]).*atan2(sqrt(d2([q(2,i);q(3,i);q(4,i)])),q(1,i)) + Der*p(:,i)));
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ODE initial conditions
opts = odeset('Reltol',1e-6,'AbsTol',1e-7,'Stats','on');
tspan = [0  30];
Quat0 = normalize([1 1 0 0]')';
Ang0 = [0 1 2];
Angv0 = [0 Ang0];
Pos0 = [0 0 0];
Vel0 = [0 0 0];

Pcon0 = Pconj(Quat0',0.5*qmat(Quat0')*(Angv0'),inertiaM4)';
v0 = [Quat0 Ang0 Pcon0 Pos0 Vel0];
[t,v] = ode113(@(t, v) odefun(t,v,inertiaM,inertiaM4,mass,Leverarm,Thrust,g,kP,kD,P2,D2),tspan,v0,opts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ODE function system of differntial equations 

function dv_dt = odefun(t,v,inertiaM,inertiaM4,mass,Leverarm,Thrust,g,kP,kD,kP2,kD2)
    Ang = [v(5);v(6);v(7)];
    Anm = inertiaM*Ang;
    Angv = [0;v(5);v(6);v(7)];
    Quat = [v(1);v(2);v(3);v(4)];
    Quac = conj(Quat);
    Qmat = qmat(Quat);
    Qmai = qmat(qinv(Quat));
    Qdot = 0.5*Qmat*Angv;
    %Pcon = Pconj(Quat,Qdot,inertiaM4);
    Pcon = [v(8);v(9);v(10);v(11)];
    %{
    Qcontrol = clamp(QC(Quat,Pcon,kP,kD));
    Tor = -Leverarm*Thrust*(cross([0;0;1],Rotate([0;0;1],Qcontrol)));
    %}
    Tor = TOR(Quat, Anm, kP2, kD2, Leverarm, Thrust, inertiaM);
    Tor(3) = 0;
    Qcontrol = Leverarm*Thrust*torquetoq(Tor/Leverarm/Thrust);
    Torq =[0;Tor];
    
    Vel = [v(15);v(16);v(17)];
    Acc = Thrust/mass*Rotate(Rotate([0;0;1],Qcontrol),conj(Quat))-g*[0;0;1];
    dv_dt = [Qdot; 
             inertiaM\(Tor+cross(Anm,Ang));
             0.5*qmat(Pcon)*Angv+(2*Qmat*Torq);  
             Vel; 
             Acc];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Quat = [v(:,1) v(:,2) v(:,3) v(:,4)]';
LEN = size(Quat, 2);
Ang = [v(:,5) v(:,6) v(:,7)]';
Anm = inertiaM*Ang;
Angv = [zeros(LEN,1) v(:,5) v(:,6) v(:,7)]';
Angm = inertiaM4*Angv;

Qmat = qmat(Quat);
Qmai = qmat(qinv(Quat));
Qdot = zeros(4,LEN); %Initialization
Pcon = [v(:,8) v(:,9) v(:,10) v(:,11)]'; 
Pmat = qmat(Pcon);
Pdot = zeros(4,LEN); %Initialization
E2 = zeros(1,LEN); %Initialization
Quac = conj(Quat);
Qcontrol = zeros(4,LEN); %Initialization
Tor = zeros(3,LEN); %Initialization
%Qcontrol = [Quat(1,:);-Quat(2,:);-Quat(3,:);zeros(1,LEN)];

for i = 1:LEN
    Qdot(:,i) = 1/2*Qmat(:,:,i)*Angv(:,i);
    %{
    Qcontrol(:,i) = QC(Quat(:,i),Pcon(:,i),P,D);
    Tor(:,i) = -Leverarm*Thrust*clampv(cross([0;0;1],Rotate([0;0;1],(Qcontrol(:,i)))));
    %}
    Tor(:,i) = TOR(Quat(:,i), Anm(:,i), P2, D2, Leverarm, Thrust, inertiaM);
    Tor(3,i) = 0;
    Qcontrol(:,i) = Leverarm*Thrust*torquetoq(Tor(:,i)/Leverarm/Thrust);
    
    
end
Torq = [zeros(1,LEN);Tor];
Euler = zeros(3,LEN);
Dcontrib = zeros(4,LEN);
for i = 1:LEN
    Pdot(:,i) = 0.5*qmat(Pcon(:,i))*Angv(:,i)+2*Qmat(:,:,i)*Torq(:,i);
    E2(:,i) = 1/2*dot(Qmat(:,:,i)*inertiaM4*Qmai(:,:,i)*Pcon(:,i),Pcon(:,i));
    Euler(1,i) = atan2(2*(Quat(1,i)*Quat(4,i) + Quat(2,i)*Quat(3,i)), 1 - 2*(Quat(3,i)^2 + Quat(4,i)^2));
    Euler(3,i) = asin(2*(Quat(1,i)*Quat(3,i) - Quat(4,i)*Quat(2,i)));
    Euler(2,i) = atan2(2*(Quat(1,i)*Quat(2,i) + Quat(3,i)*Quat(4,i)), 1 - 2*(Quat(2,i)^2 + Quat(3,i)^2));
end
E2 = 1/2*dot(Qdot,Pcon); 
E1 = 1/2*dot(Ang,inertiaM*Ang);
CH1qca = angfrom0(Qcontrol);
t_uniform = min(t):mean(diff(t))/2:max(t);
y_uniform = interp1(t, Euler(2,:), t_uniform,"spline","extrap");
% Perform FFT
Fs = 1 / mean(diff(t_uniform)); % Sampling frequency
N = length(y_uniform); % Number of data points
Y = fft(y_uniform);
P2 = abs(Y/N); % Two-sided spectrum
P1 = P2(1:floor(N/2+1)); % Single-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(N/2))/N; % Frequency vector
% Find the dominant frequency
[~, idx] = max(P1);
dominant_freq = f(idx);
disp(['Dominant Frequency: ', num2str(dominant_freq), ' Hz']);
% Fit an exponential decay model
Pos(3,end)
% Range of s values for Laplace transform
s_max = 10; % Inverse of the smallest time step
s_min = 0; % Smallest decay rate
s_values = linspace(s_min, s_max, 2000); % 100 points uniformly spaced
%{
% Initialize Laplace transform result
L_integral = zeros(size(s_values));
% Numerical Laplace Transform using MATLAB's integral function
for l = 1:length(s_values)
    s = s_values(l);
    % Numerical integration using integral function
    L_integral(l) = integral(@(t_int) interp1(t_uniform, y_uniform, t_int, "spline") .* exp(-s * t_int), min(t_uniform), max(t_uniform));
end

[~, ind] = max(abs(L_integral)); % Find the index of the maximum magnitude in the Laplace transform
decay_rate = s_values(ind); % The corresponding s value is the decay rate
disp(['Exponential Decay Rate: ', num2str(decay_rate)]);
%}
%{
fit_result = fit(t_uniform', y_uniform', 'exp1');
decay_rate = fit_result.a;
disp(['Decay Rate: ', num2str(decay_rate)]);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pos = [v(:,12) v(:,13) v(:,14)]';
Vel = [v(:,15) v(:,16) v(:,17)]';
Acc = Thrust/mass*Rotate(Rotate([0;0;1],Qcontrol),conj(Quat))-g*[0;0;1];
Basis3 = Rotate(normalize([0;0;1]),Quat);

%%%%%%%%%%%%%%%%%%%%%
figure;

[X1, Y1, Z1] = cylinder(2, 10); % Cylinder with radius 1 and 20 faces
Z1 = 10*(Z1-0.5);
[X2, Y2, Z2] = cylinder([0.5,2], 30);
Z2 = -20*(Z2);
h1 = surf(X1, Y1, Z1); % Initial plot
hold on;
h2 = surf(X2, Y2, Z2); % Additional plot
plot3(Pos(1,:),Pos(2,:),Pos(3,:),"r")
%axis equal
daspect([1 1 1]);

title('Rotating Cylinder with Quaternion Animation'); 
xlabel("x (m)")
ylabel("y (m)")
zlabel("z (m)")
%timeText = text(100, 100, 12, 'Time: 0', 'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'right'); % Initial text% Example quaternion data (replace with your actual data)
 % 4x100 array of quaternions

% Create an interpolated time vector for smoother animation
t_int = linspace(t(1), t(end), 20*round(t(end))); % 500 steps for smoothness
q_int = interp1(t, Quat', t_int, 'spline')';
q_intn = interp1(t, (Qcontrol)', t_int, 'spline')';
r_int = interp1(t, Pos', t_int, 'spline')';

% Iterate through interpolated time steps

L = 20;
for k = 1:size(t_int,2)
    %disp('C2'); 
    % Get the current quaternion
    Q_k = q_int(:, k);
    Qn_k = q_intn(:, k);
    r_k = r_int(:,k);
    % Apply rotation to cylinder vertices
    nr1 = Rotate2([X1;Y1;Z1],Q_k)+[r_k(1);r_k(1);r_k(2);r_k(2);r_k(3);r_k(3)];
    nr2 = -2.5*Rotate2([0;0;0;0;1;1],Q_k)+Rotate2([X2;Y2;Z2],multq(Q_k,Qn_k))+[r_k(1);r_k(1);r_k(2);r_k(2);r_k(3);r_k(3)];    
    xlim([-L+r_k(1) r_k(1)+L])
    ylim([-L+r_k(2) r_k(2)+L])
    zlim([-L+r_k(3) r_k(3)+L])
    %disp(int2str(size(t_int,2))); 
    % Update the surface plot
    set(h1, 'XData', [nr1(1,:);nr1(2,:)], 'YData', [nr1(3,:);nr1(4,:)], 'ZData', [nr1(5,:);nr1(6,:)]);
    set(h2, 'XData', [nr2(1,:);nr2(2,:)], 'YData', [nr2(3,:);nr2(4,:)], 'ZData', [nr2(5,:);nr2(6,:)]);
    %set(timeText, 'String', ['Time: ' num2str(t_interp(k), '%.2f')]);
     
    drawnow limitrate nocallbacks;
    
    pause(0.05);
    
end

hold off;

%}
