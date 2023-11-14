%%=======================================================================================
% The anisotropic grain size effect on the mechanical response of polycrystals:
% The role of columnar grain morphology in additively manufactured metals
% S. Amir H. Motaman, Dilay Kibaroglu
%
% Matlab Script for Axial Grain Size Measurements and Grain Size Anisotropy Calculations
%
% Note: Importing MTEX toolbox is required (tested with MTEX v5.8.1)
%%=======================================================================================

close all;
clear;
clc;

%% User inputs

RVE_file_name = 'RVE_1145Grains_GR50.txt'; % The path of the RVE's geometry file
PB = true; % false if the RVE has non-periodic boundaries
RVE_edge_length = 200; % The edge length of the RVE in µm
N_r = 128; % Number of axes to be sampled
delta_s = 0.1; % Measurement increment (as a multiple of the edge length of a voxel) for the line parameter s 

%% Uniform sampling from the axis space

% Equispace partitioning of spherical surface (uniform sampling of the axes space)
S2G = equispacedS2Grid('points',N_r,'maxRho',180*degree,'maxTheta',180*degree);

% Extracting the Cartesian and spherical coordinates of the equispaced sampled axes
r_hat = [S2G.x; S2G.y; S2G.z]'; % r_hat(:,1) = r1, r_hat(:,2) = r2, r_hat(:,3) = r3 in line with Eq.(2)
r_hat = r_hat(1:end-1,:);
theta = [S2G.theta]'/degree;
theta = theta(1:end-1,:); % Polar angle
phi = [S2G.rho]'/degree;
phi = phi(1:end-1,:); % Azimuthal angle

% Obtaining the unique axes
for q = 1 : size(theta,1)
    axes(q).theta = theta(q,1);
    axes(q).phi = phi(q,1);
    axes(q).euler_angles(1:3) = [phi(q,1)+90 theta(q,1) 0];
    axes(q).r_hat(1:3) = [r_hat(q,1) r_hat(q,2) r_hat(q,3)]; % r_hat(:,1) = r1, r_hat(:,2) = r2, r_hat(:,3) = r3 in line with Eq.(2)
end

%% Axial grain size measurement

% Importing the RVE
RVE_2D = readmatrix(RVE_file_name);
GR = min(size(RVE_2D)); % Grid resolution of the RVE
RVE = permute(reshape(RVE_2D',[GR,GR,GR]),[1,2,3]); % 3D array of microstructure (RVE)
if PB
    RVE = RVE_reconstruction(RVE);
end
N_G = max(max(max(RVE))); % Number of grains in the RVE
[GR_x, GR_y, GR_z] = size(RVE); % Dimensions of the reconstructed RVE

% Extracting the coordinates of voxel centers
x = []; 
for i = 0.5 : GR_x - 0.5
    for j = 0.5 : GR_y - 0.5
        for k = 0.5 : GR_z - 0.5
            if RVE(floor(i) + 1, floor(j) + 1, floor(k) + 1) ~= 0
                x(:, end+1) = [i; j; k]; % Coordinates of the midpoint of each voxel in the reconstructed RVE
            end
        end
    end
end
x = x';

% Loop over axes and voxels
for q = 1 : size(axes,2)
    for p = 1 : size(x,1)
        
        n = RVE(floor(x(p,1))+1, floor(x(p,2))+1, floor(x(p,3))+1); % n represents the grain ID
        
        n_current = n;        
        s_b = 0;
        pos_count = 0;        
        while n_current == n
            s_b = s_b + delta_s;
            pos_count = pos_count + 1;
            b_1 = axes(q).r_hat(1)*s_b + x(p,1); % Eq.(8)
            b_2 = axes(q).r_hat(2)*s_b + x(p,2);
            b_3 = axes(q).r_hat(3)*s_b + x(p,3);
            if b_1 <= 0 || b_1 >= GR_x || b_2 <= 0 || b_2 >= GR_y || b_3 <= 0 || b_3 >= GR_z
                break;
            end
            n_current = RVE(floor(b_1)+1, floor(b_2)+1, floor(b_3)+1);
        end
        axes(q).positive_count(p,1) = pos_count;        
        
        n_current = n;        
        s_a = 0;
        neg_count = 0;        
        while n_current == n
            s_a = s_a - delta_s;
            neg_count = neg_count + 1;
            a_1 = axes(q).r_hat(1)*s_a + x(p,1); % Eq.(8)
            a_2 = axes(q).r_hat(2)*s_a + x(p,2);
            a_3 = axes(q).r_hat(3)*s_a + x(p,3);
            if a_1 <= 0 || a_1 >= GR_x || a_2 <= 0 || a_2 >= GR_y || a_3 <= 0 || a_3 >= GR_z
                break;
            end
            n_current = RVE(floor(a_1)+1, floor(a_2)+1, floor(a_3)+1);
        end
        axes(q).negative_count(p,1) = neg_count;
        
        axes(q).min_count(p,1) = min([axes(q).positive_count(p,1) axes(q).negative_count(p,1)]);
        axes(q).total_count(p,1) = axes(q).positive_count(p,1) + axes(q).negative_count(p,1);
        
        axes(q).D_m(p,1) = (min([axes(q).positive_count(p,1) axes(q).negative_count(p,1)])) * delta_s / GR * RVE_edge_length; 
        axes(q).D(p,1) = (axes(q).positive_count(p,1) + axes(q).negative_count(p,1)) * delta_s / GR * RVE_edge_length;
    end
    
    axes(q).D_bar = mean(axes(q).D); % Mean axial grain size (normalised by the edge length of the RVE)
    axes(q).D_m_bar = mean(axes(q).D_m); % Mean minimum distance to grain boundary (normalised by the edge length of the RVE)
end

D_bar_bar = mean([axes.D_bar], 'all'); % Effective grain size (Eq.(13))
D_hat = [axes.D_bar] ./ D_bar_bar; % Measured grain size anisotropies (Eq.(15))

% Measured axial grain sizes
D = [];
for q = 1:size(axes,2)
    D = cat(1,D,(axes(q).D));
end

%% Spherical harmonics approximation of the grain size anisotropy function

nodes = [];
for i = 1:length(axes)
    nodes(i,:) = axes(i).r_hat;
end
nodes = [nodes; -nodes];
nodes = vector3d(nodes');
D_hat_prime = [D_hat'; D_hat'];

% Grain size anisotropy function
sF_D_hat = interp(nodes, D_hat_prime, 'harmonicApproximation','bandwidth',10);

% Grain size anisotropy index
A_D_hat = sqrt(norm(sF_D_hat)^2 - 4.0 * pi)

% RMSD between the grain size anisotropy function and the measured grain size anisotropies
rmsd_D = sqrt(mean((eval(sF_D_hat, nodes) - D_hat_prime).^2))

%% Plotting distribution histogram of the measured axial grain size

figure; histogram(D,50,'Normalization','probability');
set(gca,'FontSize',8,'TickLabelInterpreter','latex','LineWidth',0.6)
set(gcf, 'Units', 'Centimeters', 'Position', [10,10,8.08,6.17])
xlim([0 100])
xticks(0:20:100)
ylim([0 0.2])
yticks(0:0.05:0.2)
set(gca, 'Position', [0.1315    0.1712    0.8412    0.7895]);
xlabel('$\rm Axial\ grain\ size\ (\it{D_{ij}}\rm) \left[\mu m\right]$','Interpreter','latex')
ylabel('$\rm Probability\ (\it{P({D}_{ij})}\rm) \left[-\right]$','Interpreter','latex')

%% Streographic projection plot of the discrete measured grain size anisotropies

figure; scatter(nodes, D_hat_prime, 'upper');
colormap jet
CLim(gcm,[0.6,2.4])
colorbar('XTickLabel',{'0.6','0.9','1.2','1.5','1.8','2.1','2.4'},'XTick',0.6:0.3:2.4)

%% Streographic projection plot of the continuous grain size anisotropy function

figure; plot(sF_D_hat,'upper');
colormap jet
CLim(gcm,[0.6,2.4])
colorbar('XTickLabel',{'0.6','0.9','1.2','1.5','1.8','2.1','2.4'},'XTick',0.6:0.3:2.4)

%% Plotting the princiapl sections of the equivalent grain shape

figure; plotSection(sF_D_hat,xvector,'LineWidth',2,'Color','black');
figure; plotSection(sF_D_hat,yvector,'LineWidth',1.5,'color','black');
figure; plotSection(sF_D_hat,zvector,'LineWidth',6,'color','black');

%% Clearing dummy variables

clearvars -except axes D D_bar_bar D_hat D_hat_prime sF_D_hat...
    A_D_hat rmsd_D nodes x delta_s GR N_G N_r...
    PB RVE RVE_edge_length RVE_file_name

%% Function for reconstruction of RVEs with periodic boundaries

function [RVE_recon] = RVE_reconstruction(RVE)

GR = min(size(RVE)); % Grid resolution of the RVE
N_G = max(max(max(RVE))); % Number of grains in the RVE

% Identification of grains cutting the faces of the RVE
 
face_1 = zeros(1,N_G); % Group_1: grains cutting the faces of which i || j || k = 1
face_GR = zeros(1,N_G); % Group_GR: grains cutting the faces of which i || j || k = GR 
for n = 1 : N_G
    idx = find(RVE == n);
    [idx1,idx2,idx3] = ind2sub(size(RVE),idx);
    voxels_n = [idx1 idx2 idx3];
    
    if any(voxels_n(voxels_n == 1))
        face_1(1,n) = n;
    end
    if any(voxels_n(voxels_n == GR))
        face_GR(1,n) = n;
    end    
end

face_grains = intersect(face_1,face_GR)';  face_grains(face_grains==0) = [];

% Identification of clusters of voxels belonging to grain N

for x = 1 : size(face_grains,1)
    N = face_grains(x);
    i = find(RVE == N);
    [i1,i2,i3] = ind2sub(size(RVE),i);
    voxels_n = [i1 i2 i3];

    grain(x).ID = N;
    ID = 0; % Cluster ID
    voxels_checked = zeros(1,3);
    for q = 1 : size(voxels_n,1)
        if ~ismember(voxels_checked,voxels_n(q,:),'rows')
            ID = ID + 1;
            grain(x).cluster(ID).voxels(1,:) = voxels_n(q,:);
            voxels_checked(end+1,:) = voxels_n(q,:);
            for p = 1 : size(voxels_n,1)    
                if all(abs(voxels_n(q,:)-voxels_n(p,:))<=1) && sum(abs(voxels_n(q,:)-voxels_n(p,:))) ~= 0
                    grain(x).cluster(ID).voxels(end+1,:) = voxels_n(p,:);
                    voxels_checked(end+1,:) = voxels_n(p,:);
                end    
            end

            sz_voxels = size(grain(x).cluster(ID).voxels,1);

            p = 2;
            while p <= sz_voxels    
                for d = 1 : size(voxels_n,1)
                    if all(abs(voxels_n(d,:)-grain(x).cluster(ID).voxels(p,:))<=1) && sum(abs(voxels_n(d,:)-grain(x).cluster(ID).voxels(p,:))) ~= 0
                        if ~ismember(grain(x).cluster(ID).voxels,voxels_n(d,:),'rows')
                            grain(x).cluster(ID).voxels(end+1,:) = voxels_n(d,:);
                            sz_voxels = sz_voxels + 1;
                            voxels_checked(end+1,:) = voxels_n(d,:);
                        end
                    end
                end
                p = p + 1;
            end
        end
        
        %%% Definition of face IDs (stated by the assigned sign for each dimension) to find the clusters to be spliced
        % If the cluster is cutting the RVE at i || j || k = 1, then the i || j || k of the sign would be -1
        % If the cluster is cutting the RVE at i || j || k = GR, then the i || j || k of the sign would be 1
        % else the i || j || k of the sign would be 0
        grain(x).cluster(ID).sign = zeros(1,3);
        if any(grain(x).cluster(ID).voxels(:,1) == 1)
            grain(x).cluster(ID).sign(1) = -1;
        elseif any(grain(x).cluster(ID).voxels(:,1) == GR)
            grain(x).cluster(ID).sign(1) = 1;
        end
        if any(grain(x).cluster(ID).voxels(:,2) == 1)
            grain(x).cluster(ID).sign(2) = -1;
        elseif any(grain(x).cluster(ID).voxels(:,2) == GR)
            grain(x).cluster(ID).sign(2) = 1;
        end
        if any(grain(x).cluster(ID).voxels(:,3) == 1)
            grain(x).cluster(ID).sign(3) = -1;
        elseif any(grain(x).cluster(ID).voxels(:,3) == GR)
            grain(x).cluster(ID).sign(3) = 1;
        end
    end
    
    % Assigning priorities depending on the assigned signs 
    for d = 1 : size(grain(x).cluster,2)
        priority1 = 0;
        a = find(grain(x).cluster(d).sign == -1);
        for d2 = 1 : length(a)
            for d3 = 1 : length(grain(x).cluster)
                if grain(x).cluster(d3).sign(a(d2)) == 1
                    priority1 = priority1 + 1;
                end
            end
        end
        grain(x).cluster(d).priority1 = priority1;
        
        priority2 = 0;
        b = find(grain(x).cluster(d).sign == 1);
        for d4 = 1 : length(b)
            for d5 = 1 : size(grain(x).cluster,2)
                if grain(x).cluster(d5).sign(b(d4)) == -1
                    priority2 = priority2 + 1;
                end
            end
        end
        grain(x).cluster(d).priority2 = priority2;
    end
end

% Splicing of the clusters of voxels belonging to a certain grain ID

RVE_recon = RVE;

for x = 1 : size(grain,2)
    while size(grain(x).cluster,2) > 1
        if all([grain(x).cluster.priority1] == 0)
            break;
        end
        ID = find([grain(x).cluster.priority1] == max([grain(x).cluster.priority1]));
        id = find([grain(x).cluster(ID).priority2] == min([grain(x).cluster(ID).priority2]));
        ID = ID(id);
        ID = ID(1);            

        s = find(grain(x).cluster(ID).sign == -1);
        
        % Assigning priorities to identify in which direction the clusters would be spliced
        priority3 = [];
        priority4 = [];
        priority5 = [];
        for p = 1 : length(s)
            a = s(p);
            b = 1 : length(grain(x).cluster(ID).sign);
            b = b(b ~= a);
            c = find(grain(x).cluster(ID).sign(b) ~= 0);
            c = b(c);
            for d = 1 : size(grain(x).cluster,2)
                if grain(x).cluster(d).sign(a) == 1
                    priority3(d,p) = nnz(grain(x).cluster(d).sign(b) == grain(x).cluster(ID).sign(b));
                    priority4(d,p) = nnz(grain(x).cluster(d).sign(c) == grain(x).cluster(ID).sign(c));
                    priority5(d,p) = nnz(grain(x).cluster(d).sign(c) == 0);
                else
                    priority3(d,p) = NaN;
                    priority4(d,p) = NaN;
                    priority5(d,p) = NaN;
                end
            end
        end
        
        p3 = find(priority3 == max(priority3,[],'all'));
        p4 = find(priority4 == max(priority4,[],'all'));
        p5 = find(priority5 == max(priority5,[],'all'));
        [cpart_ID, direction] = ind2sub(size(priority3),intersect(p3,p4));
        direction = s(direction);
        
        if length(cpart_ID) > 1
            if length(intersect(intersect(p3,p4),p5)) == 1
                [cpart_ID, direction] = ind2sub(size(priority3),(intersect(intersect(p3,p4),p5)));
                direction = s(direction);
            else
                cp1 = find([grain(x).cluster(cpart_ID).priority1] == max([grain(x).cluster(cpart_ID).priority1]));
                cp2 = find([grain(x).cluster(cpart_ID(cp1)).priority2] == min([grain(x).cluster(cpart_ID(cp1)).priority2]));
                cp = cpart_ID(cp1(cp2));
                cp = cp(1);            
                direction = direction(cpart_ID == cp);
                cpart_ID = cp;
            end
        end
        
        i0=0;
        for i0 = 1 : size(grain(x).cluster(ID).voxels,1)
            RVE_recon(grain(x).cluster(ID).voxels(i0,1), grain(x).cluster(ID).voxels(i0,2), grain(x).cluster(ID).voxels(i0,3)) = 0;
        end
        
        grain(x).cluster(ID).voxels(:,direction) = grain(x).cluster(ID).voxels(:,direction) + GR;
        grain(x).cluster(cpart_ID).voxels = cat(1,grain(x).cluster(cpart_ID).voxels,grain(x).cluster(ID).voxels);
        
        if any(grain(x).cluster(cpart_ID).voxels(:,1) == 1)
            grain(x).cluster(cpart_ID).sign(1) = -1;
        elseif any(grain(x).cluster(cpart_ID).voxels(:,1) == GR)
            grain(x).cluster(cpart_ID).sign(1) = 1;
        end
        if any(grain(x).cluster(cpart_ID).voxels(:,2) == 1)
            grain(x).cluster(cpart_ID).sign(2) = -1;
        elseif any(grain(x).cluster(cpart_ID).voxels(:,2) == GR)
            grain(x).cluster(cpart_ID).sign(2) = 1;
        end
        if any(grain(x).cluster(cpart_ID).voxels(:,3) == 1)
            grain(x).cluster(cpart_ID).sign(3) = -1;
        elseif any(grain(x).cluster(cpart_ID).voxels(:,3) == GR)
            grain(x).cluster(cpart_ID).sign(3) = 1;
        end
        
        for in = 1 : size(grain(x).cluster(cpart_ID).voxels,1)
            RVE_recon(grain(x).cluster(cpart_ID).voxels(in,1), grain(x).cluster(cpart_ID).voxels(in,2), grain(x).cluster(cpart_ID).voxels(in,3)) = grain(x).ID;
        end
        
        grain(x).cluster(ID) = [];
    end
end
end
