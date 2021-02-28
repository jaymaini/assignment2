% Jay Maini 101037537
%EIPA
set(0, 'DefaultFigureWindowStyle','docked')
clear all
close all

v0 = 1;

%unequal array size
nx = 20;
ny = 30;
V = zeros(nx,ny);
G = sparse(nx*ny, nx*ny);
F = zeros(nx*ny,1);
Inclusion = 0;

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
%         elseif (j == ny)
%             F(m) = 0;
%             G(m,m) = 1;
%         elseif(j == 1)
%             F(m) = 0;
%             G(m,m) = 1;
        else
            %Inner region (set by FD equation)
            %set back to -4 for default
            G(m,m) = -4;
            
            %Define next/previous on the x and y axes
            npx = j+1 + (i-1)*ny;
            nmx = j-1 + (i-1)*ny;
            npy = j + (i-1+1)*ny;
            nmy = j + (i-1-1)*ny;
            
            G(m,npx) = 1; 
            G(m,nmx) = 1;
            G(m,npy) = 1;
            G(m,nmy) = 1;
        end
    end
end

%LU Decompsition
Varr = G\F;


    Varr_temp = reshape(Varr,nx,ny);
    figure(1)
    surf(Varr_temp)


