set(0, 'DefaultFigureWindowStyle','docked')
clear all
close all

%default values for mesh grid, box reisstance, and neck size
mg = 20;
br = 1e-2;
bs = 0.2;

%Plot current vs mesh size
mg = linspace(5,100,100/5);
for i=1:length(mg)
    [I_max(i), I_avg(i)] = current_calculator(mg(i),br,bs);
end

figure(1)
plot(mg,I_max)
hold on
plot(mg,I_avg)
title('Current vs Mesh Grid sizes')
xlabel('X-axis Mesh Grid Size (#)')
ylabel('Current (Amps)')
legend('Maximum Current', 'Average Current')
    
%plot current vs bottleneck size
I_max = 0;
I_avg = 0;
mg = 20;
bs = linspace(0.1,0.8,8);
for i=1:length(bs)
    [I_max(i), I_avg(i)] = current_calculator(mg,br,bs(i));
end

figure(2)
plot(bs,I_max)
hold on
plot(bs,I_avg)
title('Current vs Bottleneck Size')
xlabel('Ratio of opening to total height')
ylabel('Current (Amps)')
legend('Maximum Current', 'Average Current')

%plot current vs box resistance
I_max = 0;
I_avg = 0;
mg = 20;
bs = 0.2;
br = linspace(1e-6,1,100);
for i=1:length(br)
    [I_max(i), I_avg(i)] = current_calculator(mg,br(i),bs);
end

figure(3)
plot(br,I_max)
hold on
plot(br,I_avg)
title('Current vs Box Conductivity')
xlabel('Conductivity')
ylabel('Current (Amps)')
legend('Maximum Current', 'Average Current')



function [J_max, J_avg] = current_calculator(meshgridsize, boxresistance, bottleneck_size)
v0 = 1;

%unequal array size (Defaults 20 and 30)
nx = meshgridsize;
ny = round(1.5*meshgridsize);
V = zeros(nx,ny);
G = sparse(nx*ny, nx*ny);
F = zeros(nx*ny,1);

%Easy adjustment for bottleneck size - useful for function
percent_open = bottleneck_size; %0.2 default
percent_boundary = (1 - percent_open)/2;
b_height_1 = round(percent_boundary*ny);
b_height_2 = round(ny-percent_boundary*ny);

%Make conductivity array 
cMap = zeros(nx,ny) + 1;
%Create bottleneck areas
BNx = [(2/5)*nx (3/5)*nx];
%BNx = [b_height_1 b_height_2];
BNy1 = [1 b_height_1];
BNy2 = [b_height_2 ny];
%Modify cMap values

for o = BNx(1):BNx(2)
    for p = BNy1(1):BNy1(2)
        cMap(o,p) = boxresistance; %default 1e-2
    end
    
    for u = BNy2(1):BNy2(2)
        cMap(o,u) = boxresistance;
    end
end

for i = 1:nx
    for j = 1:ny
        m = j+(i-1)*ny;
        
        %Edge cases (BC's nodes)
        if (i == 1)
            F(m) = v0;
            G(m,m) = 1;
        elseif (i == nx)
            F(m) = 0;
            G(m,m) = 1;
        elseif (j == 1)
            %Bottom side - only three resistors
            npy = j+1 + (i-1)*ny;
            npx = j + (i-1+1)*ny;
            nmx = j + (i-1-1)*ny;
            nmy = j-1 + (i-1)*ny;
            
            %Calculate resistances based on conductivity map
            rpx = (cMap(i,j) + cMap(i+1,j))/2;
            rmx = (cMap(i,j) + cMap(i-1,j))/2;
            rpy = (cMap(i,j) + cMap(i,j+1))/2;

            G(m,m) = -(rpx+rmx+rpy);
            G(m,npx) = rpx; 
            G(m,nmx) = rmx;
            G(m,npy) = rpy;
        elseif (j == ny)
            %top side - only three resistors
            nmy = j-1 + (i-1)*ny;
            npx = j + (i-1+1)*ny;
            nmx = j + (i-1-1)*ny;
            
            %Calculate resistances based on conductivity map
            rpx = (cMap(i,j) + cMap(i+1,j))/2;
            rmx = (cMap(i,j) + cMap(i-1,j))/2;
            rmy = (cMap(i,j) + cMap(i,j-1))/2;

            G(m,m) = -(rpx+rmx+rmy);
            G(m,npx) = rpx; 
            G(m,nmx) = rmx;
            G(m,nmy) = rmy;
        else
            %Inner region - four resistors
            %Define next/previous on the x and y axes
            npy = j+1 + (i-1)*ny;
            nmy = j-1 + (i-1)*ny;
            npx = j + (i-1+1)*ny;
            nmx = j + (i-1-1)*ny;
            
            %Calculate resistances based on conductivity map
            rpx = (cMap(i,j) + cMap(i+1,j))/2;
            rmx = (cMap(i,j) + cMap(i-1,j))/2;
            rpy = (cMap(i,j) + cMap(i,j+1))/2;
            rmy = (cMap(i,j) + cMap(i,j-1))/2;

            G(m,m) = -(rpx+rmx+rpy+rmy);
            G(m,npx) = rpx; 
            G(m,nmx) = rmx;
            G(m,npy) = rpy;
            G(m,nmy) = rmy;
        end
    end
end

%LU Decompsition
Varr = G\F;
Varr_temp = reshape(Varr,ny,nx);

%Output the current obtained
J = cMap' .*gradient(-Varr_temp);
%Peak & average current 
J_max = max(J, [], 'all');
J_avg = mean(J, 'all');

end