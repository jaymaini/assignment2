% Jay Maini 101037537
%Analytical solution
set(0, 'DefaultFigureWindowStyle','docked')
clear all
close all

v0 = 1;

%unequal array size
nx = 20;
ny = 30;
V = zeros(nx,ny);

x = linspace(-nx/2,nx/2,nx);
y = linspace(0,ny,ny);
figure(1)

n = 1;
seriesSum = zeros(nx,ny);
a = ny;
b = nx/2;

while(n <= 90)
    for i = 1:nx
        for j = 1:ny
            seriesSum(i,j) = seriesSum(i,j) + ((cosh(n*pi*x(i)/a)/cosh(n*pi*b/a))*sin(n*pi*y(j)/a))/n;         
            V(i,j) = (4*v0/pi)*seriesSum(i,j);
        end
    end
    n = n + 2;
    surf(V)
    pause(0.05)
end


    


