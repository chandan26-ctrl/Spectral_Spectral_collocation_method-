% solving helmholz
clc; clear; close all;
x = (-1:2/80:1);
y = (-1:2/80:1);

sz = size(x);
N = sz(1,2);
xc = zeros(N,1);
for i = 1:N
    xc(i) = cos(pi*(i-1)/(N-1));
end
DX = zeros(N,N);
DX(1,1) = (2*(N-1)^2 +1)/6.0;
DX(N,N) = -DX(1,1);
c =zeros(N,1);

for i = 1:N
   if (i>1) && (i<N)
       c(i) = 1.0;
   else
       c(i) =2.0;
   end

end
for j =2:N-1
    DX(j,j) = -xc(j)/(2.0*(1.0-xc(j)^2));
end


for j =1:N
    for i = 1:N
        if i ~= j 
             DX(i,j)  = c(i)*(-1)^(i+j)/(c(j)*(xc(i) -xc(j)));
        end       
    end
end

% dd = DX*f;
yc = zeros(N,1);
for i = 1:N
    yc(i) = cos(pi*(i-1)/(N-1));
end
DY = zeros(N,N);
DY(1,1) = (2*(N-1)^2 +1)/6.0;
DY(N,N) = -DY(1,1);

for j =2:N-1
    DY(j,j) = -yc(j)/(2.0*(1.0-yc(j)^2));
end


for j =1:N
    for i = 1:N
        if i ~= j 
             DY(i,j)  = c(i)*(-1)^(i+j)/(c(j)*(yc(i) -yc(j)));
        end       
    end
end


dxx = DX*DX;
dyy = DY*DY;
I = eye(N);
   % L = kron(I,dxx) + kron(dyy,I);
   L = kron(dxx,I) + kron(I,dyy);

ES = N*N - (N-1);
EN = N*N;
B = -eye(N*N);
for i = 1:EN
    for j = 1:EN
        if  i == j
            if i <=N
                L(i,:) = 0.0;
                B(i,:) = 0.0;
                L(i,i) = 1.0;
            end
            if (ES<=i) && (i<=EN)
                 L(i,:) = 0.0;
                 B(i,:) = 0.0;
                 L(i,i) = 1.0;
            end
            if  (N<i)  && (i<ES)
                if rem(i,N) == 0
                    L(i,:) = 0;
                    B(i,:) = 0.0;
                    L(i,i) = 1.0;

                    L(i-N+1,:) = 0;
                    B(i-N+1,:) = 0.0;
                    L(i-N+1,i-N+1) = 1.0;
                end
            end
        end
    end
end


[V,D] = eig(L,B);

dd = diag(D)/(pi^2/4);

% DD = (abs(dd)).^2;
% ddr = real(dd);
% ddi = imag(dd)
% plot(ddr, ddi)
err = 1.0;
for i  = 1: (N-2)*(N-2)
      err1 = abs(dd(i) - 34);
    if err1 < err
        j = i;
        err = err1;
    end
   
end
xc = xc';
yc = yc';
[X,Y] = meshgrid(xc,yc);
[X,Y] = meshgrid(x,y);
 EIGF = V(:,j);
 plot(X, Y, '-k')
 hold on
  plot(Y, X, '-k')
  xlabel('X'); ylabel('Y');
NEIG = reshape(EIGF,N,N);
% neig1 = real(NEIG);
% NEIG = NEIG';
% contourf(X,Y,NEIG)
pcolor(X,Y,NEIG), shading interp
hold on
contour(X,Y,NEIG,15,'-.k','LineWidth',0.8)
colormap(flipud(jet))
colorbar
xlabel('X'); ylabel('Y');