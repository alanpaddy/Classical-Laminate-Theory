%Developed by Alan Ramirez. January 2017. For SESG6039
disp('Classical Laminate Theory script. ');
disp('-----------------------------------------');
disp('Welcome. Please select a script to run:');
disp(' 1: CLT - outputs ABD matrixes, stress plots in (x,y) & (1,2) coordinates)'); 
disp(' 2: Failure Predictions - Max stress, Max strain, Tsai Wu criterions');
profile = input('Enter your selection: ');

E_11=54e9; E_22=18e9; G_12=6e9; v_12=0.28;%Material inputs in [Pa] 
% theta = [-45 45 90 0 90 0 0 90 0 90 45 -45] %lab report
theta = [0 90 90 0]; %Ply angles in [degrees]
%theta = [0 45 90 90 90 45 0];
t_ply = 0.000125; %Ply thickness in [m]
% t_ply = 0.002759/12; %lab report
%%The above can be commented to uncomment the below:
%  theta = input('Enter fibre angles in degrees in a vector form eg [-45 45]: ');
%  t_ply = input('Enter thicknes of plies in m: ');
%  E_11 = input('Enter E_11 of plies in Pa: ');
%  E_22 = input('Enter E_22 of plies in Pa: ');
%  G_12 = input('Enter G_12 of plies in Pa: ');
%  v_12 = input('Enter v_12 of plies: ');

N_ply = length(theta); %amount of plyes = amount of angles entered
h_f = N_ply*t_ply; %thickness of composite = thicknes of each ply * Num of plyes
theta_rad=deg2rad(theta);



%Q for the sandwich faces
S11 = 1/E_11;
S12 = -v_12/E_11;
S22 = 1/E_22;
S66 = 1/G_12;
S = [S11 S12 0; S12 S22 0; 0 0 S66];
Q = inv(S);

%Parametric thickness of the composite
for i = 1:N_ply;
  z_ply = -h_f/2:t_ply:h_f/2;
end;

%Populate 3 arrays of 3x3 zeros
A = zeros(3,3);
B = zeros(3,3);
D = zeros(3,3);

for i = 1:N_ply; 
    %formulaes from CLT to calculate ABD
    A = A + Q_bar(theta_rad(i)) * (z_ply(i+1)-z_ply(i)); 
    B = B + Q_bar(theta_rad(i)) * 0.5*(z_ply(i+1)^2-z_ply(i)^2);
    D = D + Q_bar(theta_rad(i)) * 1/3*(z_ply(i+1)^3-z_ply(i)^3);
end;

E_x = 1/h_f*(A(1,1)-((A(1,2)^2/A(2,2))));
E_y = 1/h_f*(A(2,2)-((A(1,2)^2)/A(1,1)));
G_xy = 1/h_f*A(3,3);
v_xy = A(1,2)/A(2,2);

ABD = [A B;B D]; %ABD Matrix

%A few calculations to obtain the inverse of ABD. Note inv(ABD) can also be
%used at the price of a longer time. See report.
A_star = inv(A);
B_star = -inv(A)*B;
C_star = B*inv(A); %B/A
D_star = D - B*inv(A)*B; %D-(B/A*B)
A_prime = A_star-B_star/D_star*C_star;
B_prime = B_star/D_star;
C_prime = -inv(D_star)*C_star;
D_prime = inv(D_star);

ABD_inv_star = [A_prime, B_prime; C_prime, D_prime]; %=inv(ABD)

if profile == 1 %if the user has selected CLT mode the this code is executed
   
    Nx=25*10^3; Ny=0; Nxy=0; %Enter Loads here
    Mx=0; My=0; Mxy=0; %Enter Moments here

    N = [Nx ; Ny ; Nxy]; %populate an array
    M = [Mx ; My ; Mxy];

    epsK = ABD_inv_star*[N ; M]; % 1x6 array containing Epsion and K values
    eps_0 = epsK(1:3,1);%splitting Epsilon (strain) values
    K = epsK(4:6, 1); %splitting K (curve) values
    %Individual Epsilon & K values, for de-bugging:
    epsilon_x_0 = epsK(1,1); 
    epsilon_y_0 = epsK(2,1);
    gamma_xy_0 = epsK(3,1);
    K_x = epsK(4,1);
    K_y = epsK(5,1);  
    K_xy = epsK(6, 1);

    stress = []; 
    stress_12 = [];
    strain=[];
    strain_12=[];
    
    layer_Z = zeros(1,N_ply*2); %specifying that z will be double the amount of plies
    
    for i = 0:N_ply-1
        %Zlayer
        layer_Z(2*i+2) = z_ply(i+2);
        layer_Z(2*i+1) = z_ply(i+1);
        %Strain (x,y) coordinates
        tempx_strain = eps_0+z_ply(i+1)*(K);
        tempy_strain = eps_0+z_ply(i+2)*(K);
        strain = [strain, tempx_strain, tempy_strain];
        %Stress (x,y) coordinates
        tempx_stress = Q_bar(theta_rad(i+1))*(eps_0+z_ply(i+1)*(K));
        tempy_stress = Q_bar(theta_rad(i+1))*(eps_0+z_ply(i+2)*(K));
        stress = [stress,tempx_stress,tempy_stress];
        %Strain (1,2) coordinates
        tempx_strain_12 = eps_0+z_ply(i+1)*(K);
        tempy_strain_12 = eps_0+z_ply(i+2)*(K);
        strain_12 = [strain_12, tempx_strain_12, tempy_strain_12];
        %stress (1,2) coordinates
        tempx_stress_12 = T_mat(theta_rad(i+1))*Q_bar(theta_rad(i+1))*(eps_0+z_ply(i+1)*(K));
        tempy_stress_12 = T_mat(theta_rad(i+1))*Q_bar(theta_rad(i+1))*(eps_0+z_ply(i+2)*(K));
        stress_12 = [stress_12,tempx_stress_12,tempy_stress_12];

    end

    sigmaX = stress(1,:);
    sigmaY = stress(2,:);
    Gamma = stress(3,:);
    sigma_1 = stress_12(1,:);
    sigma_2 = stress_12(2,:);
    Gamma_12 = stress_12(3,:);

    %plot stress xy
    grid on
    hold on
    plot(sigmaX*10^-6, layer_Z*10^3,'--','LineWidth',2) 
    plot(sigmaY*10^-6,layer_Z*10^3,'k','LineWidth',2);
    plot(Gamma*10^-6,layer_Z*10^3,'-.','LineWidth',2);
    legend ('\sigma x','\sigma y','\gammaxy');
    xlabel('Stress (MPa)');
    ylabel('Laminate z (mm)');
    title('Stresses on Laminate (x,y) coord'); 

    %plot stress 12
    figure
    grid on
    hold on
    plot(sigma_1*10^-6, layer_Z*10^3,'--','LineWidth',2) 
    plot(sigma_2*10^-6,layer_Z*10^3,'k','LineWidth',2) 
    plot(Gamma_12*10^-6,layer_Z*10^3,'-.','LineWidth',2);
    legend ('\sigma 1','\sigma 2','\gamma 12');
    xlabel('Stress (MPa)');
    ylabel('Laminate z (mm)');
    title('Stresses on Laminate (1,2) coord'); 

    %Display ABD
    disp(' ');
    disp('CLT Script results: '); 
    disp('------------------------------- ');
    disp('A matrix: '); 
    disp(A)
    disp('N/m')
    disp(' ');
    disp('------------------------------- ');
    disp('B matrix: '); 
    disp(B)
    disp('N')
    disp('------------------------------- ');
    disp('D matrix: '); 
    disp(D)
    disp(['GPa mm^3'])
    disp(['------------------------------- ']);
    disp(' ');
    disp(['Exlam = ', num2str(E_x*10^-9),' GPa']);
    disp(['Eylam = ', num2str(E_y*10^-9),' GPa']);
    disp(['Gxylam = ', num2str(G_xy*10^-9),' GPa']);
    disp(['vxylam = ', num2str(v_xy)]);
    disp(['End of results']);
    disp(' ');

elseif profile == 2 %If the failure mode has beeen selected then the following is ran: 

    failure=0;

    disp(['------------------------------']);
    disp(['Please select a Failure Criterion: ']);
    disp(['1: Maximum Stress Criterion: ']);
    disp(['2: Maximum Strain Criterion']);
    disp(['3: Tsai-Wu Criterion']);
    profile_2= input('Enter your selection: ');
    disp(['Increasing Load till failure. Please wait...']);
    disp(['-----------------------------'])

    stress = [];
    stress_12 = [];
    layer_Z = zeros(1,N_ply*2);
    strain=[];
    strain_12=[];
        
    %Material Properties
    %Stress limits
    Xt = 1040*10^6;
    Yt = 35*10^6;
    Xc = 570*10^6;
    Yc = 114*10^6;
    ss = 72*10^6;
    %Strain limits
    E1t = 0.021;
    E2t = 0.002;
    E1c = 0.011;
    E2c = 0.064;
    Gamma = 0.038;
    
    %Tsai-Wu, initial calculations
    F_1 = (1/Xt) - (1/Xc);
    F_11 = 1/(Xt * Xc);
    F_2 = (1/Yt) - (1/Yc);
    F_22 = 1/(Yt * Yc);
    F_6 = 0;
    F_66 = 1/(ss^2);
    F12s = -0.5;
    F_12 = F12s * (F_11 * F_22)^0.5;
  
    
    Nx=0; Ny=0; Nxy=0; %IMPORTANT: 14900 N/m is the initial guess of the failure load.
    Mx=0; My=0; Mxy=0; %The closer Nx is set to the failure load the faster the loop will run. BUT must be below the actual failure load, so if unsure set to 0.  
    
    
    while failure==0 %while no failure has occured:
    Nx = Nx+0.5;    %Keep on increasing the load by 0.5 N/m until a failure has occured. 
    N = [Nx ; Ny ; Nxy]; %Accuracy of result is increased by decreasing 0.5 to a lower value
    M = [Mx ; My ; Mxy];
       
    epsK = ABD_inv_star*[N ; M]; %Same as before
    
    for i = 0:N_ply-1 %Same loop as before but now with failure criteria
            %Zlayer
            layer_Z(2*i+2) = z_ply(i+2);
            layer_Z(2*i+1) = z_ply(i+1);

            %Strain (1,2) coordinates
            tempx_strain_12 = T_mat(theta_rad(i+1))*(epsK(1:3,1)+z_ply(i+1)*epsK(4:6, 1));
            tempy_strain_12 = T_mat(theta_rad(i+1))*(epsK(1:3,1)+z_ply(i+2)*epsK(4:6, 1));

            %stress (1,2) coordinates
            tempx_stress_12 = T_mat(theta_rad(i+1))*Q_bar(theta_rad(i+1))*(epsK(1:3,1)+z_ply(i+1)*epsK(4:6, 1));
            tempy_stress_12 = T_mat(theta_rad(i+1))*Q_bar(theta_rad(i+1))*(epsK(1:3,1)+z_ply(i+2)*epsK(4:6, 1));

            switch profile_2
                
            case 1  %Stress Failure test
                
                if tempx_stress_12(1)>= Xt || tempy_stress_12(1) >= Xt
                    disp(['Ply No: ',num2str(i+1), ', stress failure tension 1 '])
                    failure = 1;
                elseif -tempx_stress_12(1)>= Xc || -tempy_stress_12(1)>= Xc
                    disp(['Ply No: ',num2str(i+1), ', stress failure compression 1 '])
                    failure = 1;
                elseif tempx_stress_12(2)>= Yt || tempy_stress_12(2)>= Yt
                    disp(['Ply No: ',num2str(i+1), ', stress failure tension 2 '])
                    failure=1;
                elseif -tempx_stress_12(2)>= Yc || -tempy_stress_12(2)>= Yc
                    disp(['Ply No: ',num2str(i+1), ', stress failure compression 2 '])
                    failure=1;

               elseif abs(tempx_stress_12(3))>= ss || abs(tempy_stress_12(3))>= ss
                    disp(['Ply No: ',num2str(i+1), ', stress failure shear 2 '])
                    failure=1;
                end
                
            case 2 %Strain Failure test
                
                    if  tempx_strain_12(1)>= E1t || tempy_strain_12(1) >= E1t
                        disp(['Ply No: ', num2str(i+1), ', strain failure tension 1'])
                        failure = 1;
                    elseif -tempx_strain_12(1) >= E1c || -tempy_strain_12(1) >= E1c
                        disp(['Ply No: ',num2str(i+1), ', strain failure compression 1'])
                        failure = 1;
                    elseif tempx_strain_12(2)>= E2t || tempy_strain_12(2)>= E2t
                        disp(['Ply No: ',num2str(i+1), ', strain failure tension 2'])
                        failure = 1;
                    elseif -tempx_strain_12(2) >= E2c || -tempy_strain_12(2) >= E2c
                        disp(['Ply No: ',num2str(i+1), ', strain failure compression 2'])
                        failure = 1;
                    elseif abs(tempx_strain_12(3))>= Gamma || abs(tempy_strain_12(3))>= Gamma
                        disp(['Ply No: ',num2str(i+1), ', stress failure shear 2'])
                        failure = 1;
                    end
                    
            case 3 %Strain Failure test
          
                    if ((F_1 * tempx_stress_12(1)) + (F_2 * tempx_stress_12(2)) + (F_6 * tempx_stress_12(3))... 
                        + (F_11 * (tempx_stress_12(1))^2) + (F_22 * (tempx_stress_12(2))^2)... 
                        + (F_66 * (tempx_stress_12(3))^2) + (2 * F_12 * tempx_stress_12(1) * tempx_stress_12(2))) >= 1
                        disp(['Ply No: ',num2str(i+1), ', Tsai-Wu Failure'])
                        failure = 1;
                    elseif ((F_1 * tempy_stress_12(1)) + (F_2 * tempy_stress_12(2)) + (F_6 * tempy_stress_12(3))... 
                        + (F_11 * (tempy_stress_12(1))^2) + (F_22 * (tempy_stress_12(2))^2)... 
                        + (F_66 * (tempy_stress_12(3))^2) + (2 * F_12 * tempy_stress_12(1) * tempy_stress_12(2))) >= 1
                        disp(['Ply No: ',num2str(i+1), ', Tsai-Wu Failure'])
                        failure = 1;    
                    end
            end   
    end
    
    end
    
    disp(['Failure occured at Max Load: ',num2str(Nx*10^-3), ' kN/m'])
end


%------------------------
%Sandwich plate theory for symmetrical sandwich
%------------------------


N_ply_low = N_ply;
N_ply_up = N_ply;
A_sandwich = zeros(3,3);
B_sandwich = zeros(3,3);
C_sandwich = zeros(3, 3);
D_sandwich = zeros(3,3);

A_up = zeros(3,3);
C_up = zeros(3,3);
A_low = zeros(3,3);
C_low = zeros(3,3);

a =.75; b=.75; h_c=0.025; G_xz=38e6; G_yz=38e6; E_c =110e6;

for i = 1:N_ply_up;
  z_up= h_c/2:t_ply:h_c/2+h_f;
end

for i = 1:N_ply_up; 
    %formulaes from CLT to calculate ABD
    A_up = A_up + Q_bar(theta_rad(i)) * (z_up(i+1)-z_up(i));
    C_up = C_up + Q_bar(theta_rad(i)) * ((z_up(i+1)-z_up(i))*(z_up(i)+(z_up(i+1)-z_up(i))/2));
end;

for i = 1:N_ply_low;
  z_low = (-h_c/2-h_f):t_ply:-h_c/2;
end

for i = 1:N_ply_low; 
    %formulaes from CLT to calculate ABD
    A_low = A_low + Q_bar(theta_rad(i)) * (z_low(i+1)-z_low(i));
    C_low = C_low + Q_bar(theta_rad(i)) * ((z_low(i+1)-z_low(i))*(z_low(i)+(z_low(i+1)-z_low(i))/2));
end;

A_sandwich = A_low + A_up;
B_sandwich = (h_c+h_f)/2*A_up - (h_c+h_f)/2*A_low;
C_sandwich = C_low + C_up;
D_sandwich = (h_c + h_f)/2*C_up-(h_c + h_f)/2*C_low;

m=20;
n=20;
mn_matrix=zeros(m,n);

for m=1:20
    for n=1:20
       
        F11 = m^2*pi^2*D_sandwich(1,1)/a^2 + n^2*pi^2*D_sandwich(3,3)/b^2 + h_c*G_xz;
        F12 = m*n*pi^2*(D_sandwich(1,2)+D_sandwich(3,3))/(a*b);
        F13= m*pi*h_c*G_xz/a;
        F22 = n^2*pi^2*D_sandwich(2,2)/b^2 + m^2*pi^2*D_sandwich(3,3)/a^2+h_c*G_yz;
        F23 = n*pi*h_c*G_yz/b;
       % F33(n,m) = pi()^2*(m^2*h_c*G_xz/a^2 + n^2*h_c*G_yz/b^2 + m^2*N_0/a^2 + n^2*Ny/b^2);
        F33 = (F11*F23.^2+F22*F13.^2-2*F12*F13*F23)/(F11*F22-F12.^2);

        % Global Buckling, Load
        N_0 = h_c*(G_xz+(n/m)^2*(a/b)^2*G_yz) - a^2*F33/(pi^2*m^2);
        N_0_CLT = pi^2/a^2*(D(1,1)*m^2+2*(D(1,2)+2*D(3,3))*n^2+D(2,2)*n^4/m^2);

        mn_matrix(n,m)=N_0;


    end
end
     
[n,m]=find(mn_matrix==min(min(mn_matrix)));
N_0 = min(min(mn_matrix)); %Critical load

Q_constant =0.825;
Stress_wrinkle = Q_constant*((E_c*E_x*G_xz)/(1-(v_xy^2)))^(1/3);

disp(['Critical Buckling Load ',num2str(N_0),'N'])
disp(['Wrinkle Stress', num2str(Stress_wrinkle), 'Pa'])