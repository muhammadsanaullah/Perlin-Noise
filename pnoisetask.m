
dim = 300;
% Define a noise map array to create the Perlin Noise image
noisemap = zeros(dim,dim);

global G; global M;

%256 pseudo-random unit vector table
G = zeros(256,3);
for count = 1:256
    xv = rand(); yv = rand(); zv = rand();
    xv = interp1([0,1],[-1,1],xv);
    yv = interp1([0,1],[-1,1],yv);
    zv = interp1([0,1],[-1,1],zv);
    v = [xv, yv, zv]; v = v/norm(v);
    G(count,:) = v;
end

%Pseudo-random numbers between 0 to 255 inclusive
M = randperm(255); M(end+1) = 0;


%n = noise(250,250,250);

for i = 1:dim
    for j = 1:dim
        for k = 1:5
            noisemap(i,j) = noise(i,j,k);
        end
    end
end

% %Plotting Perlin Noise
% for i = 1:10
%     for j = 1:10
%         for k = 1:3
%             noisemap(i,j,k) = noise(i,j,k);
%         end
%     end
% end

figure; imagesc(noisemap);
figure; mesh(noisemap);

% xlin = linspace(min(noisemap), max(noisemap), 10);
% ylin = linspace(min(noisemap), max(noisemap), 10);
% [X,Y] = meshgrid(xlin, ylin);
% Z = griddata(x,y,z,X,Y,'v4');
% mesh(X,Y,Z);

%Perlin Noise Function
function outputn = noise(x,y,z)

    %Mapping Collision Point Coordinates with Lattice Cells
    i = floor(x); j = floor(y); k = floor(z);
    u = x - i; v = y - j; w = z - k;

    %Weight Function
    function out = weight(a,b,c)
        out = (1 - 3*(a.^2) + 2*(a.^3))*(1 - 3*(b.^2) + 2*(b.^3))*(1 - 3*(c.^2) + 2*(c.^3));
    end

    function val = grad(i,j,k)

        global G; global M;
        arg = mod(mod(mod(M(i),256)+M(j),256)+M(k),256);
        if (arg == 0)
            arg = 1;
        end
        val = G(arg);

    end

    % Interpolating for 8 lattice points
    outputn = weight(u,v,w)*grad(i,j,k) + weight(1-u,v,w)*grad(i+1,j,k) + ... 
            weight(1-u,1-v,w)*grad(i+1,j+1,k) + weight(1-u,1-v,1-w)*grad(i+1,j+1,k+1)...
            + weight(u,1-v,1-w)*grad(i,j+1,k+1) + weight(u,v,1-w)*grad(i,j,k+1)...
            + weight(1-u,v,1-w)*grad(i+1,j,k+1) + weight(u,1-v,w)*grad(i,j+1,k);


% 
%     xc = noise_vector(1,1); yc = noise_vector(1,2); zc = noise_vector(1,3);
% 
%         %lattice points
%         % unit random vectors generate a unit Sphere with radius = 1
%         % a Bounding Box Cube of length = 2 with central point 0,0,0
%         % gives these lattice points
%         lp1 = [-1,-1,-1];
%         lp2 = [-1,-1,1];
%         lp3 = [-1,1,-1];
%         lp4 = [-1,1,1];
%         lp5 = [1,-1,-1];
%         lp6 = [1,-1,1];
%         lp7 = [1,1,-1];
%         lp8 = [1,1,1];
% 
%         %distance vectors
%         % distance vectors from the collision point to the lattice corners
%         dv1 = [(xc - lp1(1)), (yc - lp1(2)), (zc - lp1(3))];
%         dv2 = [(xc - lp2(1)), (yc - lp2(2)), (zc - lp2(3))];
%         dv3 = [(xc - lp3(1)), (yc - lp3(2)), (zc - lp3(3))];
%         dv4 = [(xc - lp4(1)), (yc - lp4(2)), (zc - lp4(3))];
%         dv5 = [(xc - lp5(1)), (yc - lp5(2)), (zc - lp5(3))];
%         dv6 = [(xc - lp6(1)), (yc - lp6(2)), (zc - lp6(3))];
%         dv7 = [(xc - lp7(1)), (yc - lp7(2)), (zc - lp7(3))];
%         dv8 = [(xc - lp8(1)), (yc - lp8(2)), (zc - lp8(3))];
% 
%         %dot products with noise gradient vectors
%         a = dot(dv1, noise_vector);
%         b = dot(dv2, noise_vector);
%         c = dot(dv3, noise_vector);
%         d = dot(dv4, noise_vector);
%         e = dot(dv5, noise_vector);
%         f = dot(dv6, noise_vector);
%         g = dot(dv7, noise_vector);
%         h = dot(dv8, noise_vector);
% 
%         %interpolation
%         function out = lerp(t,a,b)
%             out = a + t*(b-a);
%         end
% 
%         outputn = lerp(w, lerp(v, lerp(u, a, b), lerp(u, c, d)),...
%                           lerp(v, lerp(u, e, f), lerp(u, g, h)));
        
end









