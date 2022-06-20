% @ Muhammad Bin Sanaullah

% A Basic Implementation of Perlin Noise in 3D based on the algorithm and 
% procedures described in C. Basdogan, C. Ho, M. Srinivasan (1997) and
% Perlin (1985).

% parameters for gradient vector table, hashing numbers, spatial frequency
% and lattice cube dimensions
global G; global M; global freq;
freq = 1; dim = 25; % adjustable parameters

% Define a noise map array to create the Perlin Noise
noisemap = zeros(dim,dim);

% 256 pseudo-random unit vector table 
G = zeros(256,3);
for count = 1:256
    xv = rand(); yv = rand(); zv = rand();
    xv = interp1([0,1],[-1,1],xv);
    yv = interp1([0,1],[-1,1],yv);
    zv = interp1([0,1],[-1,1],zv);
    v = [xv, yv, zv]; v = v/norm(v);
    G(count,:) = v;
end

% Pseudo-random numbers between 0 to 255 inclusive 
M = randperm(255); M(end+1) = 0; 

% Testing for various points in lattice cube
cnt1 = 1;
for i = 1:(dim/255):dim
    cnt2 = 1;
    for j = 1:(dim/255):dim
        noisemap(cnt1, cnt2) = pnoise(i,j,1);
        cnt2 = cnt2 + 1;
    end
    cnt1 = cnt1 + 1;
end

%figure; imagesc(noisemap);
figure; mesh(noisemap); zlim([-2 2]);
str = sprintf('Perlin Noise with Spatial Frequency = %d & Lattice Cube Dimension = %d', freq, dim);
title(str);

% Perlin Noise Function
function outputn = pnoise(x,y,z)

    % Adjusting Spatial Frequency 
    global freq; %instead of using wavelength 
    x = x*freq; y = y*freq; z = z*freq;

    % Mapping Collision Point Coordinates with Lattice Cells 
    % to compute collision coordinates 
    i = floor(x); j = floor(y); k = floor(z);
    u = x - i; v = y - j; w = z - k;

    % Weight Function for smooth interpolation
    function out = weight(a,b,c)
        out = (1 - 3*(a.^2) + 2*(a.^3))*(1 - 3*(b.^2) + 2*(b.^3))*(1 - 3*(c.^2) + 2*(c.^3));
    end

    % using M to index G as a Permutation Table
    function val = grad(i,j,k)
        global G; global M;
        arg = mod(mod(mod(M(i),256)+M(j),256)+M(k),256);
        if (arg == 0)
            arg = 1;
        end
        val = G(arg);
    end

    % Basic Linear Interpolation
    function out = lerp(t,a,b)
        out = a + t*(b-a);
    end

    % Trilinear Interpolation for 8 Lattice Points takes the 3D input
    % to provide a scalar output that denotes the Perlin Noise at that
    % point inside the lattice cube in 3D space
    outputn = lerp(w, lerp(v, lerp(u, weight(u,v,w)*grad(i,j,k), ...
                                      weight(1-u,v,w)*grad(i+1,j,k)), ...
                              lerp(u, weight(1-u,1-v,w)*grad(i+1,j+1,k), ...
                                      weight(1-u,1-v,1-w)*grad(i+1,j+1,k+1))), ...
                      lerp(v, lerp(u, weight(u,1-v,1-w)*grad(i,j+1,k+1), ...
                                      weight(u,v,1-w)*grad(i,j,k+1)), ...
                              lerp(u, weight(1-u,v,1-w)*grad(i+1,j,k+1), ...
                                      weight(u,1-v,w)*grad(i,j+1,k))));   
end









