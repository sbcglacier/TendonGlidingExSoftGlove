% The following code is used to calculate the scale factors whic needs to be input into OpenSim for the index finger.
% The measurements made for a finger is transformed to the co-ordinate frame in OpenSim using affine transforms.

% Sung-bok Chung
% 21/06/2022

% notes on abbreviations:
% MC or CMC i- metacarpal
% PP - proximial phalange 
% MP - intermediate phalange
% DP - distal phalange

% MCP = metacarpophalangeal joint
% PIP = proximal interphalangeal joint
% DIP = distal interphalangeal joint

% In OpenSim 1-breadth(x-axis)  2-length(y-axis)  3-depth(z-axis)

%% Import data from .txt files containing finger bone mesh data from OpenSim
% Change the directory "C:\Users\18507522\Downloads\FingerMeshData" to the folder containing the F2_MC.txt file
%-------------------------------------------------------------------------------
% Import index finger MC bone
opts = delimitedTextImportOptions("NumVariables", 3);
opts.DataLines = [1, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["x", "y", "z"];
opts.VariableTypes = ["double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
F2MC = readtable("C:\Users\18507522\Downloads\FingerMeshData\F2_MC.txt", opts);
clear opts

%Import PP index finger
opts = delimitedTextImportOptions("NumVariables", 3);
opts.DataLines = [1, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["x", "y", "z"];
opts.VariableTypes = ["double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
F2PP = readtable("C:\Users\18507522\Downloads\FingerMeshData\F2_PP.txt", opts);
clear opts

%import MP index finger
opts = delimitedTextImportOptions("NumVariables", 3);
opts.DataLines = [1, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["x", "y", "z"];
opts.VariableTypes = ["double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
F2MP = readtable("C:\Users\18507522\Downloads\FingerMeshData\F2_MP.txt", opts);
clear opts

% Import index finger DP bone
opts = delimitedTextImportOptions("NumVariables", 3);
opts.DataLines = [1, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["x", "y", "z"];
opts.VariableTypes = ["double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
F2DP = readtable("C:\Users\18507522\Downloads\FingerMeshData\F2_DP.txt", opts);
clear opts

%% 
% Variable to select the parameters in the CCD design
Ds= 0.7071%-0.7071
Ls= 0.7071%-0.7071
Bs= -0.7071%-0.7071

% Variables rotate the bone about joints. The goal is to align the joint centers to be in a line.
angscalar_cmc = 0;
angscalar_pp = -2;
angscalar_mp = -6.56;
angscalar_dp = 0;

% Enable this flag to plot the parent axis system
p_flag = 0;

% Enable this flag to plot the child axis system
c_flag = 0;

% Enable this flag to plot the used reference frame axis system (default = enabled)
n_flag = 1;

% Enable this flag to plot the CMC axis system (default = enabled)
cmc_flag = 1;

figure(1);
hold on;
grid on;
axis equal;
hold off;

ax = gca;
ax.Clipping = 'off'; 

M_T = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];

%========================================================
% Import the mesh bone data to be used
%========================================================

% Assign the MC bone mesh to matrix A
mesh = F2MC;
A1(1,1:length(mesh.x)) = mesh.x(:).';
A1(2,1:length(mesh.y)) = mesh.y(:).';
A1(3,1:length(mesh.z)) = mesh.z(:).';
A1(4,1:length(mesh.z)) = 1;

hold on;
S = 0.01;

if cmc_flag == 1
plot3([0 1*S],[0 0],[0 0],'color','r','LineWidth',1);
plot3([0 0],[0 1*S],[0 0],'color','g','LineWidth',1);
plot3([0 0],[0 0],[0 1*S],'color','b','LineWidth',1);
end

%===================================
%Rotation vectors for MC
% ==================================
angscalar = angscalar_cmc * (2*pi)/360;
% Variable rot_v takes its vector from the ARMS OpenSim model (variable: 2MCP rotation2)
rot_v = [0.978577 0.0260223 -0.197329];

rot1 = angscalar*rot_v;
rot1_x =  rot1(1);% unit: radians
rot1_y =  rot1(2);
rot1_z = rot1(3);
% matrices to rotate the bone mesh
Rx1 = [1 0 0 0; 0 cos(rot1_x) -sin(rot1_x) 0; 0 sin(rot1_x) cos(rot1_x) 0; 0 0 0 1];
Ry1 = [cos(rot1_y) 0 sin(rot1_y) 0; 0 1 0 0; -sin(rot1_y) 0 cos(rot1_y) 0; 0 0 0 1];
Rz1 = [cos(rot1_z) -sin(rot1_z) 0 0; sin(rot1_z) cos(rot1_z) 0 0; 0 0 1 0; 0 0 0 1];

% Variable Txyz takes its vector from the ARMS OpenSim model (variable: 2MCP secondmc_offset)
Txyz = [0.003501 -0.028243 0.00362]; % offset from CMC to MCP
T_m1 = [1 0 0 Txyz(1); 0 1 0 Txyz(2); 0 0 1 Txyz(3); 0 0 0 1];

grid on;
axis equal;
ax = gca;
ax.Clipping = 'off';  

%========================================================
% Rotation vectors for PP
%========================================================
angscalar = angscalar_pp * (2*pi)/360; % sets flexion angle
rot_v = [0.978577 0.0260223 -0.197329]; % (OpenSim variable: 3MCP rotation2)
Txyz = [0.003501 -0.028243 0.00362];% (OpenSim variable: thirdmc_offset)

clear mesh;
clear A;
mesh = F2PP;
A2(1,1:length(mesh.x)) = mesh.x(:).';
A2(2,1:length(mesh.y)) = mesh.y(:).';
A2(3,1:length(mesh.z)) = mesh.z(:).';
A2(4,1:length(mesh.z)) = 1;
hold on;

rot2 = angscalar*rot_v;
rot2_x =  rot2(1);% this is in radians
rot2_y =  rot2(2);%
rot2_z = rot2(3);
Rx2 = [1 0 0 0; 0 cos(rot2_x) -sin(rot2_x) 0; 0 sin(rot2_x) cos(rot2_x) 0; 0 0 0 1];
Ry2 = [cos(rot2_y) 0 sin(rot2_y) 0; 0 1 0 0; -sin(rot2_y) 0 cos(rot2_y) 0; 0 0 0 1];
Rz2 = [cos(rot2_z) -sin(rot2_z) 0 0; sin(rot2_z) cos(rot2_z) 0 0; 0 0 1 0; 0 0 0 1];

% Translation matrix:
T_m2 = [1 0 0 Txyz(1); 0 1 0 Txyz(2); 0 0 1 Txyz(3); 0 0 0 1];

grid on;
axis equal;
ax = gca;
ax.Clipping = 'off';  

%========================================================
% Rotation vectors for MP
%========================================================
angscalar = angscalar_mp * (2*pi)/360; % MP flexion
rot_v = [0.976336 0.126712 -0.175394];% (OpenSim variable: 2prox_midph_b rotation1)
Txyz = [0.007517 -0.041273 0.007326]; % (OpenSim variable: 2proxph_offset)

clear mesh;
clear A;
mesh = F2MP;
A3(1,1:length(mesh.x)) = mesh.x(:).';
A3(2,1:length(mesh.y)) = mesh.y(:).';
A3(3,1:length(mesh.z)) = mesh.z(:).';
A3(4,1:length(mesh.z)) = 1;
hold on;

rot3 = angscalar*rot_v;
rot3_x =  rot3(1);% this is in radians
rot3_y =  rot3(2);%
rot3_z = rot3(3);
Rx3 = [1 0 0 0; 0 cos(rot3_x) -sin(rot3_x) 0; 0 sin(rot3_x) cos(rot3_x) 0; 0 0 0 1];
Ry3 = [cos(rot3_y) 0 sin(rot3_y) 0; 0 1 0 0; -sin(rot3_y) 0 cos(rot3_y) 0; 0 0 0 1];
Rz3 = [cos(rot3_z) -sin(rot3_z) 0 0; sin(rot3_z) cos(rot3_z) 0 0; 0 0 1 0; 0 0 0 1];

% Translation matrix:
T_m3 = [1 0 0 Txyz(1); 0 1 0 Txyz(2); 0 0 1 Txyz(3); 0 0 0 1];

grid on;
axis equal;
ax = gca;
ax.Clipping = 'off';  

%========================================================
% Rotation vectors for DP
%========================================================
angscalar = angscalar_dp * (2*pi)/360; % DP flexion
rot_v = [0.986691 0.147151 0.0749456]; % (OpenSim variable: rotation1 2mid_disph)
Txyz = [0.003333 -0.024888 0.00122]; % (OpenSim variable: 2midph_offset)

clear mesh;
clear A;
mesh = F2DP;
A4(1,1:length(mesh.x)) = mesh.x(:).';
A4(2,1:length(mesh.y)) = mesh.y(:).';
A4(3,1:length(mesh.z)) = mesh.z(:).';
A4(4,1:length(mesh.z)) = 1;

hold on;
S = 0.01;
rot4 = angscalar*rot_v;
rot4_x =  rot4(1);% this is in radians
rot4_y =  rot4(2);%
rot4_z = rot4(3);
Rx4 = [1 0 0 0; 0 cos(rot4_x) -sin(rot4_x) 0; 0 sin(rot4_x) cos(rot4_x) 0; 0 0 0 1];
Ry4 = [cos(rot4_y) 0 sin(rot4_y) 0; 0 1 0 0; -sin(rot4_y) 0 cos(rot4_y) 0; 0 0 0 1];
Rz4 = [cos(rot4_z) -sin(rot4_z) 0 0; sin(rot4_z) cos(rot4_z) 0 0; 0 0 1 0; 0 0 0 1];

% Translation matrix:
T_m4 = [1 0 0 Txyz(1); 0 1 0 Txyz(2); 0 0 1 Txyz(3); 0 0 0 1];

grid on;
axis equal;
xlabel('Global X-axis',FontWeight='bold');% (R)
ylabel('Global Y-axis',FontWeight='bold');% (G)
zlabel('Global Z-axis',FontWeight='bold');% (B)
ax = gca;
ax.Clipping = 'off';  

%===============================================
% rotate the default bones to be aligned to the desired orientation (vector along length of finger aligned with Y-axis)
%===============================================
% Adjust for yaw pitch and roll
vec1 = Rz1*Ry1*Rx1*T_m2*[0; 0; 0; 1];
vec2 = Rz1*Ry1*Rx1*T_m2*Rz2*Ry2*Rx2*T_m3*Rz3*Ry3*Rx3*T_m4*[0; 0; 0; 1];
vec_r = vec2-vec1;
rot_v = [0.978577 0.0260223 -0.197329]'; % from 2MCP rotation2
cross_prod = cross(vec_r(1:3),rot_v(1:3));

% angle between two vectors:
%---------------------------------
A = vec_r(1:3);
B =cross_prod;
ang = rad2deg(acos( dot(A,B)/(norm(A)*norm(B))));
%---------------------------------

% rotating to be in line with the x-axis
p1 = ([vec_r(1) vec_r(2) vec_r(3)]).';
size_1 = size(p1);
p1(4,1:size_1(2)) = 1;
yaw = atan2(p1(1),-p1(2));
pitch = -atan2(p1(3),sqrt(p1(2)^2+p1(1)^2));

yaw_M = [cos(yaw) -sin(yaw) 0 0; sin(yaw) cos(yaw) 0 0; 0 0 1 0; 0 0 0 1];
pitch_M = [1 0 0 0; 0 cos(pitch) -sin(pitch) 0; 0 sin(pitch) cos(pitch) 0; 0 0 0 1];
Rot_M_pitch =  pitch_M*yaw_M;
Rot_MT_pitch = Rot_M_pitch\M_T; 
p2 = Rot_MT_pitch*p1;

cross_prod2 = cross_prod;
size_1 = size(cross_prod2);
cross_prod2(4,1:size_1(2)) = 1;
p3=Rot_MT_pitch*cross_prod2;

roll = atan2(p3(1),p3(3));
Rot_M_roll = [cos(roll) 0 sin(roll) 0; 0 1 0 0; -sin(roll) 0 cos(roll) 0; 0 0 0 1]; % SB1
Rot_MT_roll = Rot_M_roll\M_T;

% Correct for yaw
p3 = Rot_M_roll\p3;

% final transform for roll and pitch
M_T = Rot_MT_roll*Rot_MT_pitch*M_T;

% =====================================================
%                    Plotting the mesh
%======================================================

% this offset_w and offset_h is used to readjust the position of the MCP such that it aligns with the origin
offset_ref =  M_T*Rz1*Ry1*Rx1*T_m2*[0 1 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 0 1; 1/S 1/S 1/S 1/S 1/S 1/S]*S; % taken from PP
offset_w = 0.00113882;% change this value such that the mcp aligns with zero
offset_h = 0.002391;%
T_offset = [1 0 0 offset_w; 0 1 0 0; 0 0 1 offset_h; 0 0 0 1];
boneCMC = T_offset*M_T*Rz1*Ry1*Rx1*A1;
plot3(boneCMC(1,:),boneCMC(2,:),boneCMC(3,:),'o','MarkerSize',1,'Color',[0 0 1]);

% =======================================
% Initialize the variables

% literature depth measurements mean values
D2DP=	13.3e-3 +Ds*1.55e-3;%1.55
D2DIP=	14.7e-3 +Ds*1.65e-3;%14.7	1.65
D2MP=	16.6e-3 +Ds*1.86e-3;%16.6	1.86
D2PIP=	19.0e-3 +Ds*2.13e-3;%19.0	2.13
D2PP=	21.4e-3 +Ds*2.39e-3;%21.4	2.39
D2MCP = 25.7e-3 +Ds*1.84e-3;%	1.84

% Bias is B and O indicates the measured values in the OpenSim model
% note: direction for depth vector goes into the direction of the palm of the hand
D2DP_B = 0.186;
D2DP_O = 0.01043;   % mult by 3.3538 ch
D2DIP_B = 0.118;
D2DIP_O = 0.0123;   % mult by 2.8703 ch
D2MP_B = 0.038;
D2MP_O = 0.0112*1.3;   % mult by 2.652 ch *1.3
D2PIP_B = -0.008;
D2PIP_O = 0.01487;  % mult by 2.3757 ch
D2PP_B = -0.031;
D2PP_O = 0.01902;   % mult by 3.1552

D2MCP_B = 0.067;
D2MCP_O = 0.0219507;        % mult by 2.0794

% literature mean values for breadth
B2DP=	16.7e-3+Bs*1.84e-3;%16.7	1.84
B2DIP=	17.5e-3+Bs*1.72e-3; %17.5	1.72
B2MP=	18.8e-3+Bs*1.82e-3; %18.8	1.82
B2PIP=	20.6e-3+Bs*1.92e-3; %20.6	1.92
B2PP=	22.4e-3+Bs*2.16e-3; %22.4	2.16

% Bias is B and O indicates the measured values in the OpenSim model
% OpenSim Breadth values: breadth vector goes in direction of joint rotation axis
B2DP_B = -0.055;
B2DP_O = 0.01326;%  multiply by 3.3318 ch
B2DIP_B = -0.123;
B2DIP_O = 0.01521;  % mult by 1.96556 ch
B2MP_B = -0.077;
B2MP_O = 0.01203*1.3; % mult by 1.9494 ch * 1.3
B2PIP_B = -0.039;
B2PIP_O = 0.01646;   % mult by 1.743 ch
B2PP_B = -0.023;
B2PP_O = 0.0171;    % mult by 1.8702 ch

% literature mean values for length
L2DP=	26.5e-3+Ls*2.18e-3;%24.5	2.18
L2MP=	30.6e-3+Ls*2.60e-3;%25.6	2.60
L2PP=	50.6e-3+Ls*3.96e-3;%45.3*	3.96
L2CMC = 80.9e-3+Ls*7.12e-3;%7.12

% OpenSim model mean values
L2DP_O = 0.02139; % ch multiply length (DIP to bone tip) value by 1.191
L2MP_O = 0.02513; % ch
L2PP_O = 0.0425847; % ch
L2CMC_O = 0.05987;% ch

% this cmc offset value used to adjust Y direction marker points to align with joints
CMC_OFFSET = 7.2736e-3+0.0212924;%0.0262665 % important to reposition - from midpoint of cmc  only changes along length of finger


%~~~~~~~~~~~~~~~~~~~~~~~ MC SCALE FACTOR ~~~~~~~~~~~~~~~~~~~~~~~~~~
boneCMC = T_offset*M_T*Rz1*Ry1*Rx1*A1;
plot3(boneCMC(1,:),boneCMC(2,:),boneCMC(3,:),'o','MarkerSize',1,'Color',[0 0 1]);
OS = 0;
Rev = T_offset*M_T*Rz1*Ry1*Rx1;

V_MC = Rev\[0,0,0,0,0,0;...
            -CMC_OFFSET+OS,-CMC_OFFSET+OS,-CMC_OFFSET+OS,-CMC_OFFSET+OS+L2CMC_O,-CMC_OFFSET+OS+L2CMC_O,-CMC_OFFSET+OS+L2CMC_O;...
            topbias(D2MCP_B,D2MCP_O),-botbias(D2MCP_B,D2MCP_O),0,topbias(D2MCP_B,D2MCP_O),-botbias(D2MCP_B,D2MCP_O),0;...
            1,1,1,1,1,1];

E_MC = Rev\[0,0,0,0,0,0;...
            -CMC_OFFSET+OS,-CMC_OFFSET+OS,-CMC_OFFSET+OS,-CMC_OFFSET+OS+L2CMC,-CMC_OFFSET+OS+L2CMC,-CMC_OFFSET+OS+L2CMC;...
            topbias(D2MCP_B,D2MCP),-botbias(D2MCP_B,D2MCP),0,topbias(D2MCP_B,D2MCP),-botbias(D2MCP_B,D2MCP),0;...
            1,1,1,1,1,1];

TP = [0,0,0,0,0,0;...
            -CMC_OFFSET+OS,-CMC_OFFSET+OS,-CMC_OFFSET+OS,-CMC_OFFSET+OS+L2CMC_O,-CMC_OFFSET+OS+L2CMC_O,-CMC_OFFSET+OS+L2CMC_O;...
            topbias(D2MCP_B,D2MCP_O),-botbias(D2MCP_B,D2MCP_O),0,topbias(D2MCP_B,D2MCP_O),-botbias(D2MCP_B,D2MCP_O),0;...
            1,1,1,1,1,1];

plot3(TP(1,:),TP(2,:),TP(3,:),'o','color','g');

marker_cnt = length(V_MC(1,:));
marker_cnt_x = 0;
marker_cnt_y = 0;
marker_cnt_z = 0;
s_totx = 0;
s_toty = 0;
s_totz = 0;

for i = 1:marker_cnt
if V_MC(1,i) ~= 0
if E_MC(1,i)/V_MC(1,i) <3
if E_MC(1,i)/V_MC(1,i) >0
s_totx = s_totx + E_MC(1,i)/V_MC(1,i);
marker_cnt_x = marker_cnt_x +1;
end
end
end

if V_MC(2,i) ~= 0 
if E_MC(2,i)/V_MC(2,i) <3
if E_MC(2,i)/V_MC(2,i) >0
s_toty = s_toty + E_MC(2,i)/V_MC(2,i);
marker_cnt_y = marker_cnt_y +1;
end
end
end

if V_MC(3,i) ~= 0
if E_MC(3,i)/V_MC(3,i) <3
if E_MC(3,i)/V_MC(3,i) >0
s_totz = s_totz + E_MC(3,i)/V_MC(3,i);
marker_cnt_z = marker_cnt_z +1;
end
end
end
end

s_totx = s_totx/marker_cnt_x;
s_toty = s_toty/marker_cnt_y;
s_totz = s_totz/marker_cnt_z;
Scale_MC = [s_totx; s_toty; s_totz]

% ~~~~~~~~~~~~~~~~~~~ PP SCALE FACTOR ~~~~~~~~~~~~~~~
OS = -CMC_OFFSET;
Rev = T_offset*M_T*Rz1*Ry1*Rx1*T_m1*Rz2*Ry2*Rx2;
E_PP = Rev\[0,0,topbias(B2PP_B,B2PP),-botbias(B2PP_B,B2PP),0,0,0,topbias(B2PIP_B,B2PIP),-botbias(B2PIP_B,B2PIP),0;...
        -L2PP/2+OS,-L2PP/2+OS,-L2PP/2+OS,-L2PP/2+OS,-L2PP/2+OS,-L2PP+OS,-L2PP+OS,-L2PP+OS,-L2PP+OS,-L2PP+OS;...
        topbias(D2PP_B,D2PP),-botbias(D2PP_B,D2PP),0,0,0,topbias(D2PIP_B,D2PIP),-botbias(D2PIP_B,D2PIP),0,0,0;...
        1,1,1,1,1,1,1,1,1,1];

OS = -CMC_OFFSET;
Rev = T_offset*M_T*Rz1*Ry1*Rx1*T_m1*Rz2*Ry2*Rx2;
V_PP = Rev\[0,0,topbias(B2PP_B,B2PP_O),-botbias(B2PP_B,B2PP_O),0,0,0,topbias(B2PIP_B,B2PIP_O),-botbias(B2PIP_B,B2PIP_O),0;...
        -L2PP_O/2+OS,-L2PP_O/2+OS,-L2PP_O/2+OS,-L2PP_O/2+OS,-L2PP_O/2+OS,-L2PP_O+OS,-L2PP_O+OS,-L2PP_O+OS,-L2PP_O+OS,-L2PP_O+OS;...
        topbias(D2PP_B,D2PP_O),-botbias(D2PP_B,D2PP_O),0,0,0,topbias(D2PIP_B,D2PIP_O),-botbias(D2PIP_B,D2PIP_O),0,0,0;...
        1,1,1,1,1,1,1,1,1,1];

TP = [0,0,topbias(B2PP_B,B2PP_O),-botbias(B2PP_B,B2PP_O),0,0,0,topbias(B2PIP_B,B2PIP_O),-botbias(B2PIP_B,B2PIP_O),0;...
        -L2PP_O/2+OS,-L2PP_O/2+OS,-L2PP_O/2+OS,-L2PP_O/2+OS,-L2PP_O/2+OS,-L2PP_O+OS,-L2PP_O+OS,-L2PP_O+OS,-L2PP_O+OS,-L2PP_O+OS;...
        topbias(D2PP_B,D2PP_O),-botbias(D2PP_B,D2PP_O),0,0,0,topbias(D2PIP_B,D2PIP_O),-botbias(D2PIP_B,D2PIP_O),0,0,0;...
        1,1,1,1,1,1,1,1,1,1];

plot3(TP(1,:),TP(2,:),TP(3,:),'o','color','r');
marker_cnt = length(V_PP(1,:));
marker_cnt_x = 0;
marker_cnt_y = 0;
marker_cnt_z = 0;
s_totx = 0;
s_toty = 0;
s_totz = 0;

for i = 1:marker_cnt
if V_PP(1,i) ~= 0
if E_PP(1,i)/V_PP(1,i) <3
if E_PP(1,i)/V_PP(1,i) >0
s_totx = s_totx + E_PP(1,i)/V_PP(1,i);
marker_cnt_x = marker_cnt_x +1;
end
end
end


if V_PP(2,i) ~= 0 
if E_PP(2,i)/V_PP(2,i) <3
if E_PP(2,i)/V_PP(2,i) >0
s_toty = s_toty + E_PP(2,i)/V_PP(2,i);
marker_cnt_y = marker_cnt_y +1;
end
end
end

if V_PP(3,i) ~= 0
if E_PP(3,i)/V_PP(3,i) <3
if E_PP(3,i)/V_PP(3,i) >0
s_totz = s_totz + E_PP(3,i)/V_PP(3,i);
marker_cnt_z = marker_cnt_z +1;
end
end
end
end

s_totx = s_totx/marker_cnt_x;
s_toty = s_toty/marker_cnt_y;
s_totz = s_totz/marker_cnt_z;
Scale_PP = [s_totx; s_toty; s_totz]

%Parent reference
if p_flag == 1
a_v =  T_offset*M_T*Rz1*Ry1*Rx1*T_m2*[0 1 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 0 1; 1/S 1/S 1/S 1/S 1/S 1/S]*S;
plot3(a_v(1,1:2),a_v(2,1:2),a_v(3,1:2),'--','color','r','LineWidth',1);%r
plot3(a_v(1,3:4),a_v(2,3:4),a_v(3,3:4),'--','color','g','LineWidth',1);%g
plot3(a_v(1,5:6),a_v(2,5:6),a_v(3,5:6),'--','color','b','LineWidth',1);%b
end

%Child reference
if c_flag == 1
a_v = T_offset*M_T*Rz1*Ry1*Rx1*T_m2*Rz2*Ry2*Rx2*[0 1 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 0 1; 1/S 1/S 1/S 1/S 1/S 1/S]*S;
plot3(a_v(1,1:2),a_v(2,1:2),a_v(3,1:2),'color','r','LineWidth',1);
plot3(a_v(1,3:4),a_v(2,3:4),a_v(3,3:4),'color','g','LineWidth',1);
plot3(a_v(1,5:6),a_v(2,5:6),a_v(3,5:6),'color','b','LineWidth',1);
end

if n_flag ==1
a_v = T_offset*M_T*Rz1*Ry1*Rx1*T_m2*[0; 0; 0; 1/S]*S;
plot3([a_v(1) a_v(1)+S],[a_v(2) a_v(2)],[a_v(3) a_v(3)],'color','r','LineWidth',1);
plot3([a_v(1) a_v(1)],[a_v(2) a_v(2)+S],[a_v(3) a_v(3)],'g','LineWidth',1);
plot3([a_v(1) a_v(1)],[a_v(2) a_v(2)],[a_v(3) a_v(3)+S],'color','b','LineWidth',1);
end

bonePP = T_offset*M_T*Rz1*Ry1*Rx1*T_m1*Rz2*Ry2*Rx2*A2;
plot3(bonePP(1,:),bonePP(2,:),bonePP(3,:),'o','MarkerSize',1,'Color',[0 0.4 0.7]);

% ~~~~~~~~~~~~~~~~~~~ MP SCALE FACTOR ~~~~~~~~~~~~~~~
OS = -CMC_OFFSET-L2PP_O;
Rev = T_offset*M_T*Rz1*Ry1*Rx1*T_m2*Rz2*Ry2*Rx2*T_m3*Rz3*Ry3*Rx3;
E_MP = Rev\[0,0,topbias(B2MP_B,B2MP),-botbias(B2MP_B,B2MP),0,0,0,topbias(B2DIP_B,B2DIP),-botbias(B2DIP_B,B2DIP),0;...
        -L2MP/2+OS,-L2MP/2+OS,-L2MP/2+OS,-L2MP/2+OS,-L2MP/2+OS,-L2MP+OS,-L2MP+OS,-L2MP+OS,-L2MP+OS,-L2MP+OS;...
        topbias(D2MP_B,D2MP),-botbias(D2MP_B,D2MP),0,0,0,topbias(D2DIP_B,D2DIP),-botbias(D2DIP_B,D2DIP),0,0,0;...
        1,1,1,1,1,1,1,1,1,1];

OS = -CMC_OFFSET-L2PP_O; 
Rev = T_offset*M_T*Rz1*Ry1*Rx1*T_m2*Rz2*Ry2*Rx2*T_m3*Rz3*Ry3*Rx3;
V_MP = Rev\[0,0,topbias(B2MP_B,B2MP_O),-botbias(B2MP_B,B2MP_O),0,0,0,topbias(B2DIP_B,B2DIP_O),-botbias(B2DIP_B,B2DIP_O),0;...
        -L2MP_O/2+OS,-L2MP_O/2+OS,-L2MP_O/2+OS,-L2MP_O/2+OS,-L2MP_O/2+OS,-L2MP_O+OS,-L2MP_O+OS,-L2MP_O+OS,-L2MP_O+OS,-L2MP_O+OS;...
        topbias(D2MP_B,D2MP_O),-botbias(D2MP_B,D2MP_O),0,0,0,topbias(D2DIP_B,D2DIP_O),-botbias(D2DIP_B,D2DIP_O),0,0,0;...
        1,1,1,1,1,1,1,1,1,1];

TP2 = [0,0,topbias(B2MP_B,B2MP_O),-botbias(B2MP_B,B2MP_O),0,0,0,topbias(B2DIP_B,B2DIP_O),-botbias(B2DIP_B,B2DIP_O),0;...
        -L2MP_O/2+OS,-L2MP_O/2+OS,-L2MP_O/2+OS,-L2MP_O/2+OS,-L2MP_O/2+OS,-L2MP_O+OS,-L2MP_O+OS,-L2MP_O+OS,-L2MP_O+OS,-L2MP_O+OS;...
        topbias(D2MP_B,D2MP_O),-botbias(D2MP_B,D2MP_O),0,0,0,topbias(D2DIP_B,D2DIP_O),-botbias(D2DIP_B,D2DIP_O),0,0,0;...
        1,1,1,1,1,1,1,1,1,1];

plot3(TP2(1,:),TP2(2,:),TP2(3,:),'o','color','b');

marker_cnt = length(V_MP(1,:));
s_totx = 0;
s_toty = 0;
s_totz = 0;
marker_cnt_x = 0;
marker_cnt_y = 0;
marker_cnt_z = 0;

for i = 1:marker_cnt
if V_MP(1,i) ~= 0
if E_MP(1,i)/V_MP(1,i) <3
if E_MP(1,i)/V_MP(1,i) >0
s_totx = s_totx + E_MP(1,i)/V_MP(1,i);
marker_cnt_x = marker_cnt_x +1;
end
end
end

if V_MP(2,i) ~= 0 
if E_MP(2,i)/V_MP(2,i) <3
if E_MP(2,i)/V_MP(2,i) >0
s_toty = s_toty + E_MP(2,i)/V_MP(2,i);
marker_cnt_y = marker_cnt_y +1;
end
end
end

if V_MP(3,i) ~= 0
if E_MP(3,i)/V_MP(3,i) <3
if E_MP(3,i)/V_MP(3,i) >0
s_totz = s_totz + E_MP(3,i)/V_MP(3,i);
marker_cnt_z = marker_cnt_z +1;
end
end
end
end

s_totx = s_totx/marker_cnt_x;
s_toty = s_toty/marker_cnt_y;
s_totz = s_totz/marker_cnt_z;

% Scale_MP contains the MP scale factors
fprintf('\nScale value for MP:\n')
Scale_MP = [s_totx; s_toty; s_totz]

%Parent reference
if p_flag == 1
a_v = T_offset*M_T*Rz1*Ry1*Rx1*T_m2*Rz2*Ry2*Rx2*T_m3*[0 1 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 0 1; 1/S 1/S 1/S 1/S 1/S 1/S]*S;
plot3(a_v(1,1:2),a_v(2,1:2),a_v(3,1:2),'--','color','r','LineWidth',1);
plot3(a_v(1,3:4),a_v(2,3:4),a_v(3,3:4),'--','color','g','LineWidth',1);
plot3(a_v(1,5:6),a_v(2,5:6),a_v(3,5:6),'--','color','b','LineWidth',1);
end
%Child reference
if c_flag == 1
a_v = T_offset*M_T*Rz1*Ry1*Rx1*T_m2*Rz2*Ry2*Rx2*T_m3*Rx3*Rz3*Ry3*[0 1 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 0 1; 1/S 1/S 1/S 1/S 1/S 1/S]*S;
plot3(a_v(1,1:2),a_v(2,1:2),a_v(3,1:2),'color','r','LineWidth',1);
plot3(a_v(1,3:4),a_v(2,3:4),a_v(3,3:4),'color','g','LineWidth',1);
plot3(a_v(1,5:6),a_v(2,5:6),a_v(3,5:6),'color','b','LineWidth',1);
end

if n_flag ==1
a_v = T_offset*M_T*Rz1*Ry1*Rx1*T_m2*Rz2*Ry2*Rx2*T_m3*[0; 0; 0; 1/S]*S;
plot3([a_v(1) a_v(1)+S],[a_v(2) a_v(2)],[a_v(3) a_v(3)],'color','r','LineWidth',1);
plot3([a_v(1) a_v(1)],[a_v(2) a_v(2)+S],[a_v(3) a_v(3)],'g','LineWidth',1);
plot3([a_v(1) a_v(1)],[a_v(2) a_v(2)],[a_v(3) a_v(3)+S],'color','b','LineWidth',1);
end

boneMP = T_offset*M_T*Rz1*Ry1*Rx1*T_m2*Rz2*Ry2*Rx2*T_m3*Rz3*Ry3*Rx3*A3;
plot3(boneMP(1,:),boneMP(2,:),boneMP(3,:),'o','MarkerSize',1,'Color','b');

% ~~~~~~~~~~~~~~~~~~~ DP SCALE FACTOR ~~~~~~~~~~~~~~~

OS = -CMC_OFFSET-L2MP_O-L2PP_O; 
Rev = (T_offset*M_T*Rz1*Ry1*Rx1*T_m2*Rz2*Ry2*Rx2*T_m3*Rz3*Ry3*Rx3*T_m4*Rz4*Ry4*Rx4);
E_DP =  Rev\[0,0,topbias(B2DP_B,B2DP),-botbias(B2DP_B,B2DP),0,0;...
        -L2DP/2+OS,-L2DP/2+OS,-L2DP/2+OS,-L2DP/2+OS,-L2DP/2+OS,-L2DP+OS;...
        topbias(D2DP_B,D2DP),-botbias(D2DP_B,D2DP),0,0,0,0;...
        1,1,1,1,1,1];
hold on;
%
OS = -CMC_OFFSET-L2MP_O-L2PP_O; % This is same because starting point is relative to axis !!!!
Rev = (T_offset*M_T*Rz1*Ry1*Rx1*T_m2*Rz2*Ry2*Rx2*T_m3*Rz3*Ry3*Rx3*T_m4*Rz4*Ry4*Rx4);
V_DP = Rev\[0,0,topbias(B2DP_B,B2DP_O),-botbias(B2DP_B,B2DP_O),0,0;...
        -L2DP_O/2+OS,-L2DP_O/2+OS,-L2DP_O/2+OS,-L2DP_O/2+OS,-L2DP_O/2+OS,-L2DP_O+OS;...
        topbias(D2DP_B,D2DP_O),-botbias(D2DP_B,D2DP_O),0,0,0,0;...
        1,1,1,1,1,1];

TP = [0,0,topbias(B2DP_B,B2DP_O),-botbias(B2DP_B,B2DP_O),0,0;...
        -L2DP_O/2+OS,-L2DP_O/2+OS,-L2DP_O/2+OS,-L2DP_O/2+OS,-L2DP_O/2+OS,-L2DP_O+OS;...
        topbias(D2DP_B,D2DP_O),-botbias(D2DP_B,D2DP_O),0,0,0,0;...
        1,1,1,1,1,1];
plot3(TP(1,:),TP(2,:),TP(3,:),'o','color','r');

marker_cnt = length(V_DP(1,:));
s_totx = 0;
s_toty = 0;
s_totz = 0;
marker_cnt_x = 0;
marker_cnt_y = 0;
marker_cnt_z = 0;

for i = 1:marker_cnt
if V_DP(1,i) ~= 0
if E_DP(1,i)/V_DP(1,i) <3
if E_DP(1,i)/V_DP(1,i) >0
s_totx = s_totx + E_DP(1,i)/V_DP(1,i);
marker_cnt_x = marker_cnt_x +1;
end
end
end

if V_DP(2,i) ~= 0 
if E_DP(2,i)/V_DP(2,i) <3
if E_DP(2,i)/V_DP(2,i) >0
s_toty = s_toty + E_DP(2,i)/V_DP(2,i);
marker_cnt_y = marker_cnt_y +1;
end
end
end

if V_DP(3,i) ~= 0
if E_DP(3,i)/V_DP(3,i) <3
if E_DP(3,i)/V_DP(3,i) >0
s_totz = s_totz + E_DP(3,i)/V_DP(3,i);
marker_cnt_z = marker_cnt_z +1;
end
end
end
end
s_totx = s_totx/marker_cnt_x
s_toty = s_toty/marker_cnt_y
s_totz = s_totz/marker_cnt_z
% Scale_DP presents the DP scale factors
Scale_DP = [s_totx; s_toty; s_totz]

%Parent reference
if p_flag == 1
a_v = T_offset*M_T*Rz1*Ry1*Rx1*T_m2*Rz2*Ry2*Rx2*T_m3*Rz3*Ry3*Rx3*T_m4*[0 1 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 0 1; 1/S 1/S 1/S 1/S 1/S 1/S]*S;
plot3(a_v(1,1:2),a_v(2,1:2),a_v(3,1:2),'--','color','r','LineWidth',1);
plot3(a_v(1,3:4),a_v(2,3:4),a_v(3,3:4),'--','color','g','LineWidth',1);
plot3(a_v(1,5:6),a_v(2,5:6),a_v(3,5:6),'--','color','b','LineWidth',1);
end

%Child reference
if c_flag ==1
a_v = T_offset*M_T*Rz1*Ry1*Rx1*T_m2*Rz2*Ry2*Rx2*T_m3*Rx3*Rz3*Ry3*T_m4*Rx4*Rz4*Ry4*[0 1 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 0 1; 1/S 1/S 1/S 1/S 1/S 1/S]*S;
plot3(a_v(1,1:2),a_v(2,1:2),a_v(3,1:2),'color','r','LineWidth',1);
plot3(a_v(1,3:4),a_v(2,3:4),a_v(3,3:4),'color','g','LineWidth',1);
plot3(a_v(1,5:6),a_v(2,5:6),a_v(3,5:6),'color','b','LineWidth',1);
end

% New reference
if n_flag ==1
a_v = T_offset*M_T*Rz1*Ry1*Rx1*T_m2*Rz2*Ry2*Rx2*T_m3*Rz3*Ry3*Rx3*T_m4*[0; 0; 0; 1/S]*S;
plot3([a_v(1) a_v(1)+S],[a_v(2) a_v(2)],[a_v(3) a_v(3)],'color','r','LineWidth',1);
plot3([a_v(1) a_v(1)],[a_v(2) a_v(2)+S],[a_v(3) a_v(3)],'g','LineWidth',1);
plot3([a_v(1) a_v(1)],[a_v(2) a_v(2)],[a_v(3) a_v(3)+S],'color','b','LineWidth',1);
end

boneDP = T_offset*M_T*Rz1*Ry1*Rx1*T_m2*Rz2*Ry2*Rx2*T_m3*Rz3*Ry3*Rx3*T_m4*Rz4*Ry4*Rx4*A4;
plot3(boneDP(1,:),boneDP(2,:),boneDP(3,:),'o','MarkerSize',1,'Color',[0 0.4 0.7]);

xlim([-0.02 0.02])
zlim([-0.02 0.02])
%view(0,90) % view top
view(-90,0) % view side

fprintf('Scale factors:\n')
Scale_MC

Scale_PP

Scale_MP

Scale_DP
%% Code which plots the Bone with markers translated back to the OpenSim frame
% This code allows you to visualize how the markers have been relocated to the axes in which it will be scaled in OpenSim

mesh = F2PP; % F2PP F2MP F2DP
An(1,1:length(mesh.x)) = mesh.x(:).';
An(2,1:length(mesh.y)) = mesh.y(:).';
An(3,1:length(mesh.z)) = mesh.z(:).';
An(4,1:length(mesh.z)) = 1;

figure(14)
hold on
plot3(An(1,:),An(2,:),An(3,:),'o','MarkerSize',1,'Color',[0 0 1]);
plot3(V_PP(1,:),V_PP(2,:),V_PP(3,:),'o','color','b')
plot3(E_PP(1,:),E_PP(2,:),E_PP(3,:),'o','color','r')
view(-90,0)
grid on;
hold off
axis equal

figure(15)
hold on
plot3(An(1,:),An(2,:),An(3,:),'o','MarkerSize',1,'Color',[0 0 1]);
plot3(V_PP(1,:),V_PP(2,:),V_PP(3,:),'o','color','b')
plot3(E_PP(1,:),E_PP(2,:),E_PP(3,:),'o','color','r')
view(0,90)
grid on;
hold off
axis equal


%% Top view
view ([-90 90])

%% Side view
view ([-90 0])

%% Orthogonal view
view([-44 38])
