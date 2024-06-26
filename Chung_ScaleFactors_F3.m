
% The following code is used to calculate the scale factors whic needs to be input into OpenSim for the middle finger.
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
% Change the directory "C:\Users\18507522\Downloads\FingerMeshData" to the folder containing the F3_MC.txt file
%-------------------------------------------------------------------------------
% Import middle finger MC bone
opts = delimitedTextImportOptions("NumVariables", 3);
opts.DataLines = [1, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["x", "y", "z"];
opts.VariableTypes = ["double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
F3MC = readtable("C:\Users\18507522\Downloads\FingerMeshData\F3_MC.txt", opts);
clear opts

%Import PP middle finger
opts = delimitedTextImportOptions("NumVariables", 3);
opts.DataLines = [1, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["x", "y", "z"];
opts.VariableTypes = ["double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
F3PP = readtable("C:\Users\18507522\Downloads\FingerMeshData\F3_PP.txt", opts);
clear opts

%import MP middle finger
opts = delimitedTextImportOptions("NumVariables", 3);
opts.DataLines = [1, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["x", "y", "z"];
opts.VariableTypes = ["double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
F3MP = readtable("C:\Users\18507522\Downloads\FingerMeshData\F3_MP.txt", opts);
clear opts

% Import middle finger DP bone
opts = delimitedTextImportOptions("NumVariables", 3);
opts.DataLines = [1, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["x", "y", "z"];
opts.VariableTypes = ["double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
F3DP = readtable("C:\Users\18507522\Downloads\FingerMeshData\F3_DP.txt", opts);
clear opts

%% 
% Variable to select the parameters in the CCD design
Ds= 0%0.7071%-0.7071
Ls= 0%0.7071%-0.7071
Bs= 0%-0.7071%-0.7071

% Variables rotate the bone about joints. The goal is to align the joint centers to be in a line.
angscalar_cmc = 0;
angscalar_pp = 0;
angscalar_mp = -4.56;
angscalar_dp = 0;

% Enable this flag to plot the parent axis system
p_flag = 0;

% Enable this flag to plot the child axis system
c_flag = 0;

% Enable this flag to plot the used reference frame axis system (default = enabled)
n_flag = 1;

% Enable this flag to plot the CMC axis system (default = enabled)
cmc_flag = 1;

figure(1)
clf(1)
hold on;
grid on;
axis equal

hold off;

ax = gca;
ax.Clipping = 'off'; 

M_T = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];

%========================================================
% Import the mesh bone data to be used
%========================================================

% Assign the MC bone mesh to matrix A
mesh = F3MC;
A1(1,1:length(mesh.x)) = mesh.x(:).';
A1(2,1:length(mesh.y)) = mesh.y(:).';
A1(3,1:length(mesh.z)) = mesh.z(:).';
A1(4,1:length(mesh.z)) = 1;

hold on;
S = 0.01;

if cmc_flag == 1
plot3([0 1*S],[0 0],[0 0],'color','r','LineWidth',1)
plot3([0 0],[0 1*S],[0 0],'color','g','LineWidth',1)
plot3([0 0],[0 0],[0 1*S],'color','b','LineWidth',1)
end

%===================================
%Rotation vectors for MC
% ==================================
angscalar = angscalar_cmc * (2*pi)/360;
% Variable rot_v takes its vector from the ARMS OpenSim model (variable: 5MCP rotation2)
rot_v = [0.998349 0.0473371 -0.0104129];

rot1 = angscalar*rot_v;
rot1_x =  rot1(1);% unit: radians
rot1_y =  rot1(2);%
rot1_z = rot1(3);
% matrices to rotate the bone mesh
Rx1 = [1 0 0 0; 0 cos(rot1_x) -sin(rot1_x) 0; 0 sin(rot1_x) cos(rot1_x) 0; 0 0 0 1];
Ry1 = [cos(rot1_y) 0 sin(rot1_y) 0; 0 1 0 0; -sin(rot1_y) 0 cos(rot1_y) 0; 0 0 0 1];
Rz1 = [cos(rot1_z) -sin(rot1_z) 0 0; sin(rot1_z) cos(rot1_z) 0 0; 0 0 1 0; 0 0 0 1];

% Variable Txyz takes its vector from the ARMS OpenSim model (variable: 5MCP fourthmc_offset)
Txyz = [0.00024 -0.02629 0.001778]; % offset from CMC to MCP
T_m1 = [1 0 0 Txyz(1); 0 1 0 Txyz(2); 0 0 1 Txyz(3); 0 0 0 1];

grid on;
axis equal
ax = gca;
ax.Clipping = 'off';  

%========================================================
% Rotation vectors for PP
%========================================================
angscalar = angscalar_pp * (2*pi)/360; % sets flexion angle
rot_v = [0.998349 0.0473371 -0.0104129]; % (OpenSim variable: 3MCP rotation2)
Txyz = [0.00024 -0.02629 0.001778];% (OpenSim variable: fourthmc_offset) moves joint to MCP

clear mesh
clear A
mesh = F3PP;
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
axis equal
ax = gca;
ax.Clipping = 'off';  


%========================================================
% Rotation vectors for MP
%========================================================
angscalar = angscalar_mp * (2*pi)/360; % MP flexion
rot_v = [0.992829 -0.062002 0.102203];% (OpenSim variable: 5prox_midph_b rotation1) % PIP rot
Txyz = [0.00165 -0.044211 0.00623]; % (OpenSim variable: 5proxph_offset) offsets from MCP to PIP

clear mesh
clear A
mesh = F3MP;
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
axis equal
ax = gca;
ax.Clipping = 'off';  

%========================================================
% Rotation vectors for DP
%========================================================
angscalar = angscalar_dp * (2*pi)/360; % DP flexion
rot_v = [0.994589 0.036103 0.097409]; % (OpenSim variable: rotation1 5mid_disph)
Txyz = [0.001365 -0.029047 0.001954]; % (OpenSim variable: 5midph_offset)

clear mesh
clear A
mesh = F3DP;
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
axis equal
xlabel('Global X-axis',FontWeight='bold')% (R)
ylabel('Global Y-axis',FontWeight='bold')% (G)
zlabel('Global Z-axis',FontWeight='bold')% (B)
ax = gca;
ax.Clipping = 'off';  

%===============================================
% rotate the default bones to be aligned to the desired orientation (vector along length of finger aligned with Y-axis)
%===============================================
% Adjust for yaw pitch and roll
vec1 = Rz1*Ry1*Rx1*T_m2*[0; 0; 0; 1];
vec2 = Rz1*Ry1*Rx1*T_m2*Rz2*Ry2*Rx2*T_m3*Rz3*Ry3*Rx3*T_m4*[0; 0; 0; 1];
vec_r = vec2-vec1;
rot_v = [0.998349 0.0473371 -0.0104129]'; % from 3MCP rotation2
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
Rot_M_pitch =  pitch_M*yaw_M;% *yaw_M
Rot_MT_pitch = Rot_M_pitch\M_T; 
p2 = Rot_MT_pitch*p1;

cross_prod2 = cross_prod;
size_1 = size(cross_prod2);
cross_prod2(4,1:size_1(2)) = 1;
p3=Rot_MT_pitch*cross_prod2;

roll = atan2(p3(1),p3(3));
Rot_M_roll = [cos(roll) 0 sin(roll) 0; 0 1 0 0; -sin(roll) 0 cos(roll) 0; 0 0 0 1]; % SB1
Rot_MT_roll = Rot_M_roll\M_T

% Correct for yaw
p3 = Rot_M_roll\p3;

% final transform for roll and pitch
M_T = Rot_MT_roll*Rot_MT_pitch*M_T

% =====================================================
%                    Plotting the mesh
%======================================================

% this offset_w and offset_h is used to readjust the position of the MCP such that it aligns with the origin
offset_ref =  M_T*Rz1*Ry1*Rx1*T_m2*[0 1 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 0 1; 1/S 1/S 1/S 1/S 1/S 1/S]*S; % taken from PP
offset_w = 0% change this value such that the mcp aligns with zero X-axis
offset_h = 0.00196916% change in Z-axis
T_offset = [1 0 0 offset_w; 0 1 0 0; 0 0 1 offset_h; 0 0 0 1]


% =======================================
% Initialize the variables

% literature depth measurements mean values
D3DP=	11.8e-3+Ds*1.52e-3;%c	
D3DIP=	12.7e-3+Ds*1.42e-3;
D3MP=	14.1e-3+Ds*1.71e-3;
D3PIP=	16.6e-3+Ds*1.66e-3;
D3PP=	17.3e-3+Ds*2.00e-3;

D3MCP = 25.5e-3 +Ds*2.85e-3;%

% Bias is B and O indicates the measured values in the OpenSim model
% note: direction for depth vector goes into the direction of the palm of the hand
D3DP_B = 0.186;
D3DP_O = 0.01147;
D3DIP_B = 0.118;
D3DIP_O = 0.0124;
D3MP_B = 0.038;
D3MP_O = 0.01387;
D3PIP_B = -0.008;
D3PIP_O = 0.01587;
D3PP_B = -0.031;
D3PP_O = 0.01707;

D3MCP_B = 0.067;
D3MCP_O = 0.02064;        % mult by 2.0794

% literature mean values for breadth
B3DP=	16.5e-3+Bs*1.83e-3;
B3DIP=	16.8e-3+Bs*1.62e-3;
B3MP=	17.6e-3+Bs*1.99e-3;
B3PIP=	19.4e-3+Bs*2.02e-3;
B3PP=	18.6e-3+Bs*2.08e-3;

% Bias is B and O indicates the measured values in the OpenSim model
% OpenSim Breadth values: breadth vector goes in direction of joint rotation axis
B3DP_B = -0.055;
B3DP_O = 0.01486;
B3DIP_B = -0.123;
B3DIP_O = 0.01541;
B3MP_B = -0.077;
B3MP_O = 0.01581;
B3PIP_B = -0.039;
B3PIP_O = 0.01743;
B3PP_B = -0.023;
B3PP_O = 0.01743;

% literature mean values for length
L3DP=	26.5e-3+Ls*2.47e-3;
L3MP=	30.6e-3+Ls*2.74e-3;
L3PP=	50.6e-3+Ls*3.79e-3;
L3CMC = 76.8e-3+Ls*6.69e-3;%66.5	6.57 

% OpenSim model mean values
L3DP_O = 0.02216;
L3MP_O = 0.0291;
L3PP_O = 0.0447;
L3CMC_O = 0.05725;

% this cmc offset value used to adjust Y direction marker points to align with joints
CMC_OFFSET = 0.0262665 % important to reposition - from midpoint of cmc
PP_O = -0.0709446
DP_O = -0.100089
DP_End_O = 0.118714

%~~~~~~~~~~~~~~~~~~~~~~~ MC SCALE FACTOR ~~~~~~~~~~~~~~~~~~~~~~~~~~
boneCMC = (T_offset*M_T*Rz1*Ry1*Rx1)*A1;
plot3(boneCMC(1,:),boneCMC(2,:),boneCMC(3,:),'o','MarkerSize',1,'Color',[0 0 1]);
OS = 0;
Rev = T_offset*M_T*Rz1*Ry1*Rx1;

V_MC = Rev\[0,0,0,0,0,0;...
            -CMC_OFFSET+OS,-CMC_OFFSET+OS,-CMC_OFFSET+OS,-CMC_OFFSET+OS+L3CMC_O,-CMC_OFFSET+OS+L3CMC_O,-CMC_OFFSET+OS+L3CMC_O;...
            topbias(D3MCP_B,D3MCP_O),-botbias(D3MCP_B,D3MCP_O),0,topbias(D3MCP_B,D3MCP_O),-botbias(D3MCP_B,D3MCP_O),0;...
            1,1,1,1,1,1];

E_MC = Rev\[0,0,0,0,0,0;...
            -CMC_OFFSET+OS,-CMC_OFFSET+OS,-CMC_OFFSET+OS,-CMC_OFFSET+OS+L3CMC,-CMC_OFFSET+OS+L3CMC,-CMC_OFFSET+OS+L3CMC;...
            topbias(D3MCP_B,D3MCP),-botbias(D3MCP_B,D3MCP),0,topbias(D3MCP_B,D3MCP),-botbias(D3MCP_B,D3MCP),0;...
            1,1,1,1,1,1];

TP = [0,0,0,0,0,0;...
            -CMC_OFFSET+OS,-CMC_OFFSET+OS,-CMC_OFFSET+OS,-CMC_OFFSET+OS+L3CMC_O,-CMC_OFFSET+OS+L3CMC_O,-CMC_OFFSET+OS+L3CMC_O;...
            topbias(D3MCP_B,D3MCP_O),-botbias(D3MCP_B,D3MCP_O),0,topbias(D3MCP_B,D3MCP_O),-botbias(D3MCP_B,D3MCP_O),0;...
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
    fprintf('\n')
    fprintf('a')
    E_MC(2,i)
    V_MC(2,i)
    ratio_d = E_MC(2,i)/V_MC(2,i)
    fprintf('b')
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
OS = -CMC_OFFSET
Rev = T_offset*M_T*Rz1*Ry1*Rx1*T_m1*Rz2*Ry2*Rx2
E_PP = Rev\[0,0,topbias(B3PP_B,B3PP),-botbias(B3PP_B,B3PP),0,0,0,topbias(B3PIP_B,B3PIP),-botbias(B3PIP_B,B3PIP),0;...
        -L3PP/2+OS,-L3PP/2+OS,-L3PP/2+OS,-L3PP/2+OS,-L3PP/2+OS,-L3PP+OS,-L3PP+OS,-L3PP+OS,-L3PP+OS,-L3PP+OS;...
        topbias(D3PP_B,D3PP),-botbias(D3PP_B,D3PP),0,0,0,topbias(D3PIP_B,D3PIP),-botbias(D3PIP_B,D3PIP),0,0,0;...
        1,1,1,1,1,1,1,1,1,1];

OS = -CMC_OFFSET
Rev = T_offset*M_T*Rz1*Ry1*Rx1*T_m1*Rz2*Ry2*Rx2
V_PP = Rev\[0,0,topbias(B3PP_B,B3PP_O),-botbias(B3PP_B,B3PP_O),0,0,0,topbias(B3PIP_B,B3PIP_O),-botbias(B3PIP_B,B3PIP_O),0;...
        -L3PP_O/2+OS,-L3PP_O/2+OS,-L3PP_O/2+OS,-L3PP_O/2+OS,-L3PP_O/2+OS,-L3PP_O+OS,-L3PP_O+OS,-L3PP_O+OS,-L3PP_O+OS,-L3PP_O+OS;...
        topbias(D3PP_B,D3PP_O),-botbias(D3PP_B,D3PP_O),0,0,0,topbias(D3PIP_B,D3PIP_O),-botbias(D3PIP_B,D3PIP_O),0,0,0;...
        1,1,1,1,1,1,1,1,1,1];

TP = [0,0,topbias(B3PP_B,B3PP_O),-botbias(B3PP_B,B3PP_O),0,0,0,topbias(B3PIP_B,B3PIP_O),-botbias(B3PIP_B,B3PIP_O),0;...
        -L3PP_O/2+OS,-L3PP_O/2+OS,-L3PP_O/2+OS,-L3PP_O/2+OS,-L3PP_O/2+OS,-L3PP_O+OS,-L3PP_O+OS,-L3PP_O+OS,-L3PP_O+OS,-L3PP_O+OS;...
        topbias(D3PP_B,D3PP_O),-botbias(D3PP_B,D3PP_O),0,0,0,topbias(D3PIP_B,D3PIP_O),-botbias(D3PIP_B,D3PIP_O),0,0,0;...
        1,1,1,1,1,1,1,1,1,1];

plot3(TP(1,:),TP(2,:),TP(3,:),'o','color','r');
% calculates scale factor
marker_cnt = length(V_PP(1,:))
s_totx = 0;
s_toty = 0;
s_totz = 0;
marker_cnt_x = 0;
marker_cnt_y = 0;
marker_cnt_z = 0;

for i = 1:marker_cnt
if V_PP(1,i) ~= 0
s_totx = s_totx + E_PP(1,i)/V_PP(1,i);
marker_cnt_x = marker_cnt_x +1;
end

if V_PP(2,i) ~= 0 
s_toty = s_toty + E_PP(2,i)/V_PP(2,i);
marker_cnt_y = marker_cnt_y +1;
end

if V_PP(3,i) ~= 0
s_totz = s_totz + E_PP(3,i)/V_PP(3,i);
marker_cnt_z = marker_cnt_z +1;
end
end

s_totx = s_totx/marker_cnt
s_toty = s_toty/marker_cnt
s_totz = s_totz/marker_cnt
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
OS = -CMC_OFFSET-L3PP_O
Rev = T_offset*M_T*Rz1*Ry1*Rx1*T_m2*Rz2*Ry2*Rx2*T_m3*Rz3*Ry3*Rx3
E_MP = Rev\[0,0,topbias(B3MP_B,B3MP),-botbias(B3MP_B,B3MP),0,0,0,topbias(B3DIP_B,B3DIP),-botbias(B3DIP_B,B3DIP),0;...
        -L3MP/2+OS,-L3MP/2+OS,-L3MP/2+OS,-L3MP/2+OS,-L3MP/2+OS,-L3MP+OS,-L3MP+OS,-L3MP+OS,-L3MP+OS,-L3MP+OS;...
        topbias(D3MP_B,D3MP),-botbias(D3MP_B,D3MP),0,0,0,topbias(D3DIP_B,D3DIP),-botbias(D3DIP_B,D3DIP),0,0,0;...
        1,1,1,1,1,1,1,1,1,1];

OS = -CMC_OFFSET-L3PP_O
Rev = T_offset*M_T*Rz1*Ry1*Rx1*T_m2*Rz2*Ry2*Rx2*T_m3*Rz3*Ry3*Rx3
V_MP = Rev\[0,0,topbias(B3MP_B,B3MP_O),-botbias(B3MP_B,B3MP_O),0,0,0,topbias(B3DIP_B,B3DIP_O),-botbias(B3DIP_B,B3DIP_O),0;...
        -L3MP_O/2+OS,-L3MP_O/2+OS,-L3MP_O/2+OS,-L3MP_O/2+OS,-L3MP_O/2+OS,-L3MP_O+OS,-L3MP_O+OS,-L3MP_O+OS,-L3MP_O+OS,-L3MP_O+OS;...
        topbias(D3MP_B,D3MP_O),-botbias(D3MP_B,D3MP_O),0,0,0,topbias(D3DIP_B,D3DIP_O),-botbias(D3DIP_B,D3DIP_O),0,0,0;...
        1,1,1,1,1,1,1,1,1,1];

TP2 = [0,0,topbias(B3MP_B,B3MP_O),-botbias(B3MP_B,B3MP_O),0,0,0,topbias(B3DIP_B,B3DIP_O),-botbias(B3DIP_B,B3DIP_O),0;...
        -L3MP_O/2+OS,-L3MP_O/2+OS,-L3MP_O/2+OS,-L3MP_O/2+OS,-L3MP_O/2+OS,-L3MP_O+OS,-L3MP_O+OS,-L3MP_O+OS,-L3MP_O+OS,-L3MP_O+OS;...
        topbias(D3MP_B,D3MP_O),-botbias(D3MP_B,D3MP_O),0,0,0,topbias(D3DIP_B,D3DIP_O),-botbias(D3DIP_B,D3DIP_O),0,0,0;...
        1,1,1,1,1,1,1,1,1,1];

plot3(TP2(1,:),TP2(2,:),TP2(3,:),'o','color','b');
% Final part to calculate the scale factors
marker_cnt = length(V_MP(1,:))
s_totx = 0;
s_toty = 0;
s_totz = 0;
marker_cnt_x = 0;
marker_cnt_y = 0;
marker_cnt_z = 0;

for i = 1:marker_cnt
if V_MP(1,i) ~= 0
s_totx = s_totx + E_MP(1,i)/V_MP(1,i);
marker_cnt_x = marker_cnt_x +1;
end

if V_MP(2,i) ~= 0 
s_toty = s_toty + E_MP(2,i)/V_MP(2,i);
marker_cnt_y = marker_cnt_y +1;
end

if V_MP(3,i) ~= 0
s_totz = s_totz + E_MP(3,i)/V_MP(3,i);
marker_cnt_z = marker_cnt_z +1;
end
end

s_totx = s_totx/marker_cnt
s_toty = s_toty/marker_cnt
s_totz = s_totz/marker_cnt
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

% new reference
if n_flag ==1
a_v = T_offset*M_T*Rz1*Ry1*Rx1*T_m2*Rz2*Ry2*Rx2*T_m3*[0; 0; 0; 1/S]*S;
plot3([a_v(1) a_v(1)+S],[a_v(2) a_v(2)],[a_v(3) a_v(3)],'color','r','LineWidth',1);
plot3([a_v(1) a_v(1)],[a_v(2) a_v(2)+S],[a_v(3) a_v(3)],'g','LineWidth',1);
plot3([a_v(1) a_v(1)],[a_v(2) a_v(2)],[a_v(3) a_v(3)+S],'color','b','LineWidth',1);
end


boneMP = T_offset*M_T*Rz1*Ry1*Rx1*T_m2*Rz2*Ry2*Rx2*T_m3*Rz3*Ry3*Rx3*A3;
plot3(boneMP(1,:),boneMP(2,:),boneMP(3,:),'o','MarkerSize',1,'Color','b');

% ~~~~~~~~~~~~~~~~~~~ DP SCALE FACTOR ~~~~~~~~~~~~~~~
OS = -CMC_OFFSET-L3MP_O-L3PP_O 
Rev = (T_offset*M_T*Rz1*Ry1*Rx1*T_m2*Rz2*Ry2*Rx2*T_m3*Rz3*Ry3*Rx3*T_m4*Rz4*Ry4*Rx4)
E_DP =  Rev\[0,0,topbias(B3DP_B,B3DP),-botbias(B3DP_B,B3DP),0,0;...
        -L3DP/2+OS,-L3DP/2+OS,-L3DP/2+OS,-L3DP/2+OS,-L3DP/2+OS,-L3DP+OS;...
        topbias(D3DP_B,D3DP),-botbias(D3DP_B,D3DP),0,0,0,0;...
        1,1,1,1,1,1];

OS = -CMC_OFFSET-L3MP_O-L3PP_O % This is same because starting point is relative to axis !!!!
Rev = (T_offset*M_T*Rz1*Ry1*Rx1*T_m2*Rz2*Ry2*Rx2*T_m3*Rz3*Ry3*Rx3*T_m4*Rz4*Ry4*Rx4)
V_DP = Rev\[0,0,topbias(B3DP_B,B3DP_O),-botbias(B3DP_B,B3DP_O),0,0;...
        -L3DP_O/2+OS,-L3DP_O/2+OS,-L3DP_O/2+OS,-L3DP_O/2+OS,-L3DP_O/2+OS,-L3DP_O+OS;...
        topbias(D3DP_B,D3DP_O),-botbias(D3DP_B,D3DP_O),0,0,0,0;...
        1,1,1,1,1,1];

TP = [0,0,topbias(B3DP_B,B3DP_O),-botbias(B3DP_B,B3DP_O),0,0;...
        -L3DP_O/2+OS,-L3DP_O/2+OS,-L3DP_O/2+OS,-L3DP_O/2+OS,-L3DP_O/2+OS,-L3DP_O+OS;...
        topbias(D3DP_B,D3DP_O),-botbias(D3DP_B,D3DP_O),0,0,0,0;...
        1,1,1,1,1,1];
plot3(TP(1,:),TP(2,:),TP(3,:),'o','color','r');
marker_cnt = length(V_DP(1,:))
s_totx = 0;
s_toty = 0;
s_totz = 0;
marker_cnt_x = 0;
marker_cnt_y = 0;
marker_cnt_z = 0;

for i = 1:marker_cnt
if V_DP(1,i) ~= 0
s_totx = s_totx + E_DP(1,i)/V_DP(1,i);
marker_cnt_x = marker_cnt_x +1;
end

if V_DP(2,i) ~= 0 
s_toty = s_toty + E_DP(2,i)/V_DP(2,i);
marker_cnt_y = marker_cnt_y +1;
end

if V_DP(3,i) ~= 0
s_totz = s_totz + E_DP(3,i)/V_DP(3,i);
marker_cnt_z = marker_cnt_z +1;
end
end
s_totx = s_totx/marker_cnt
s_toty = s_toty/marker_cnt
s_totz = s_totz/marker_cnt
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

% new reference
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

view(-90,0) % view side

fprintf('Scale factors:\n')
Scale_MC

Scale_PP

Scale_MP

Scale_DP
%% Code which plots the Bone with markers translated back to the OpenSim frame
% This code allows you to visualize how the markers have been relocated to the axes in which it will be scaled in OpenSim

FD = F5PP;
V_FD = V_PP;
E_FD = E_PP;

mesh = FD; % F5PP F5MP F5DP
An(1,1:length(mesh.x)) = mesh.x(:).';
An(2,1:length(mesh.y)) = mesh.y(:).';
An(3,1:length(mesh.z)) = mesh.z(:).';
An(4,1:length(mesh.z)) = 1;

figure(14)
ax = gca;
ax.Clipping = 'off';  

hold on
plot3(An(1,:),An(2,:),An(3,:),'o','MarkerSize',1,'Color',[0 0 1]);
plot3(V_FD(1,:),V_FD(2,:),V_FD(3,:),'o','color','b')
plot3(E_FD(1,:),E_FD(2,:),E_FD(3,:),'o','color','r')
view(-90,0)
ax = gca;
ax.Clipping = 'off';  

grid on;
hold off
axis equal

figure(15)
ax = gca;
ax.Clipping = 'off';  

hold on
plot3(An(1,:),An(2,:),An(3,:),'o','MarkerSize',1,'Color',[0 0 1]);
plot3(V_FD(1,:),V_FD(2,:),V_FD(3,:),'o','color','b')
plot3(E_FD(1,:),E_FD(2,:),E_FD(3,:),'o','color','r')
view(0,90)
ax = gca;
ax.Clipping = 'off';  

grid on;
hold off
axis equal



%% Top view
view ([-90 90])

%% Side view
view ([-90 0])

%% Orthogonal view
view([-44 38])
