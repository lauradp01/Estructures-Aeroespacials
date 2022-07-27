%-------------------------------------------------------------------------%
% ASSIGNMENT 03 - (A)
%-------------------------------------------------------------------------%
% Date:
% Author/s:
%

clear;
close all;
clc

%% INPUT DATA

% Material properties
E = 85e9; %Pa

% Cross-section parameters
t1 = 1.5e-3; %m
t2 = 4e-3; %m
h1 = 500e-3; %m
h2 = 250e-3; %m
b = 775e-3; %m


% Other data
g = 9.81;
L1 = 5;
L2 = 10;
L = L1+L2;
Me = 2550;
M = 35000;

% Number of elements for each part
hel = [5,2.5,1.25,0.625,0.3125] ;
%nel = L./hel; 
nel = [4,8] ; %Number of beams of each part
Nel = 96; % Number of elements for the "exact" solution

%% PRECOMPUTATIONS

% Compute section: 
% A  - Section area 
% Iz - Section inertia

A = h1*t2 + h2*t2 + (((h1-h2)/2)^2 + (b)^2)^0.5 * t1 * 2 ;
z_cdg = (h2*t2*b + (((h1-h2)/2)^2 + (b)^2)^0.5 * t1 * 2 * b/2) / A ;

Izz1 = 1/12 * t2 * h1^3 ;
Izz2 = 1/12 * t2 * h2^3 ;
theta3 = atan(b / ((h1-h2)/2)) ;
Izz3 = 1/12 * t1 * ((((h1-h2)/2)^2 + (b)^2)^0.5)^3 *sin(theta3)^2 + (((h1-h2)/2)^2 + (b)^2)^0.5 * t1 * ((h1+h2)/4)^2 ;
Izz4 = Izz3 ; 

Iz = Izz1 + Izz2 + Izz3 + Izz4 ; 


% Compute parameter l:
% l - Equilibrium parameter

syms x_int

lambda1 = M / (4*(L1+L2)) + 3*M*(L1-x_int) / (2*L2^2) ;
lambda2 = M / (4*(L1+L2)) ;
q1_no_l = 0.8 - 0.2*cos(pi*x_int/L1) ;
q2_no_l = (1 - (x_int-L1)/L2)*(1 + (x_int-L1)/L2) ;

W = double(int(lambda1*g,[0 L1])) + lambda2*(L-L1) + Me*g ;
q_no_l = double( int(q1_no_l,0,L1) + int(q2_no_l,L1,L) ) ;
l = W / q_no_l ;

% Plot analytical solution
fig = plotBeamsInitialize(L1+L2);

% Loop through each of the number of elements
for k = 1:length(nel)

    %% PREPROCESS
    
    nd = 1 ; 
    nne = 2 ; 
    nnod = sum(nel) + 1 ;
    ni = 2 ;
    ndof = nnod * ni ;

    % Nodal coordinates
    %  x(a,j) = coordinate of node a in the dimension j
    % Complete the coordinates

    x = zeros(nnod,nd) ;
    for i = 1:nel(1)
        x(i) = L1/nel(1) * (i-1) ;
    end
    for i = 1:(nel(2)+1)
        x(i+nel(1)) = L1 +  L2/nel(2) * (i-1) ;
    end
    
    % Nodal connectivities  
    %  Tnod(e,a) = global nodal number associated to node a of element e

    Tnod = zeros(sum(nel),nne) ;
    for i = 1:sum(nel)
        for j = 1:nne
            Tnod(i,j) = i + j - 1 ;
        end
    end

    Td = zeros(sum(nel),nne*ni) ; 
    for i = 1:sum(nel)
        for j = 1:nne
            column = j*ni ;
            valor = ni * Tnod(i,j) ;
            for a = 1:ni
                Td(i,column) = valor ;
                column = column-1 ;
                valor = valor-1 ; 
            end
        end
    end

    
    % Material properties matrix
    %  mat(m,1) = Young modulus of material m
    %  mat(m,2) = Section area of material m
    %  mat(m,3) = Section inertia of material m
    mat = [% Young M.        Section A.    Inertia 
              E,                A,         Iz;  % Material (1)
    ];

    % Material connectivities
    %  Tmat(e) = Row in mat corresponding to the material associated to element e 
    Tmat = ones(sum(nel),1) ;

        
    %% SOLVER
    
    % Compute:
    % u  - Displacements and rotations vector [ndof x 1]
    % pu - Polynomial coefficients for displacements for each element [nel x 4]
    % pt - Polynomial coefficients for rotations for each element [nel x 3]
    % Fy - Internal shear force at each elements's nodes [nel x nne]
    % Mz - Internal bending moment at each elements's nodes [nel x nne]

    u = zeros(ndof,1) ; 
    pu = zeros(sum(nel),4) ; 
    pt = zeros(sum(nel),3) ; 
    Fy = zeros(sum(nel),nne) ;
    Mz = zeros(sum(nel),nne) ;

    Kel = zeros(nne*ni,nne*ni,sum(nel)) ;

    for e = 1:sum(nel)
        x1e = x(Tnod(e,1),1);
        x2e = x(Tnod(e,2),1);
        le = abs(x2e-x1e) ;
        Ke = mat(Tmat(e),3)*mat(Tmat(e),1)/le^3 * [12 6*le -12 6*le ;
            6*le 4*le^2 -6*le 2*le^2 ;
            -12 -6*le 12 -6*le ;
            6*le 2*le^2 -6*le 4*le^2] ;
        for r = 1:(nne*ni)
            for s = 1:(nne*ni)
                Kel(r,s,e) = Ke(r,s) ;
            end
        end
    end
    
    Fel = zeros(nne*ni,sum(nel)) ;

    for e = 1:sum(nel)
        x1e = x(Tnod(e,1),1);
        x2e = x(Tnod(e,2),1);
        le = abs(x2e-x1e) ;
        if x1e < L1
            qe_media = int(q1_no_l*l,x1e,x2e) / le ;
        else
            qe_media = int(q2_no_l*l,x1e,x2e) / le ;
        end
        Fe = qe_media*le/2 * [1 ; le/6 ; 1 ; -le/6] ;
        
        for r = 1:nne*ni
            Fel(r,e) = Fe(r) ;
        end

    end

    Fext = zeros(nnod*ni,1) ;
    KG = zeros(nnod*ni,nnod*ni) ;

    for e = 1:sum(nel)
        for i = 1:nne*ni
            I = Td(e,i) ; 
            Fext(I) = Fext(I) + Fel(i,e) ; 
            for j = 1:nne*ni
                J = Td(e,j) ; 
                KG(I,J) = KG(I,J) + Kel(i,j,e) ; 
            end
        end
    end

    vR = [1 ; 2] ; 
    vL = zeros(ndof-length(vR),1) ;
    for i = 1:length(vL) 
        vL(i) = i+2 ;
    end
    uR = [0 ; 0 ] ; 

    KLL = KG(vL,vL) ; 
    KLR = KG(vL,vR) ; 
    KRL = KG(vR,vL) ; 
    KRR = KG(vR,vR) ; 
    Fext_L = Fext(vL,1) ; 
    Fext_R = Fext(vR,1) ; 

    uL = KLL\(Fext_L-KLR*uR) ; 
    R_R = KRR*uR + KRL*uL - Fext_R ;

    u(vL,1) = uL ; 
    u(vR,1) = uR ; 

    for e = 1:sum(nel)
        x1e = x(Tnod(e,1),1);
        x2e = x(Tnod(e,2),1);
        le = abs(x2e-x1e) ;
        ue = zeros(nne*ni,1) ; 
        for i = 1:nne*ni
            I = Td(e,i) ; 
            ue(i,1) = u(I) ; 
        end

        Fe_int = Kel(:,:,e)*ue ;

        Fy(e,1) = -Fe_int(1) ;
        Fy(e,2) = Fe_int(3) ; 
        Mz(e,1) = -Fe_int(2) ; 
        Mz(e,2) = Fe_int(4) ; 

        abcd = 1/le^3 * [2 le -2 le ;
            -3*le -2*le^2 3*le -le^2 ;
            0 le^3 0 0 ; 
            le^3 0 0 0 ] *ue ;

        pu(e,:) = abcd ;
        pt(e,:) = [3*abcd(1),2*abcd(2),abcd(3)] ;

    end

    
    %% POSTPROCESS
    
    % Number of subdivisions and plots
    nsub = Nel/nel(k);
    plotBeams1D(fig,x,Tnod,nsub,pu,pt,Fy,Mz)
    drawnow;
    
end

% Add figure legends
figure(fig)
legend(strcat('N=',cellstr(string(nel))),'location','northeast');