clear all;
clc;
close all;

%[x,y,z] = sphere(50);

N = 10;
% thetavec = logspace(0,log10(pi+1),50)-1;
% thetavec = [0,thetavec,pi];
phivec = linspace(0,2*pi,2*N);
dphi = phivec(2) - phivec(1);
thetavec = [0];
dthetaPole = pi/200;
dtheta = dthetaPole
thetavec = [thetavec(end),thetavec(end) + dtheta];
while thetavec(end) < pi
    if thetavec(end) > (pi - dthetaPole);
        break;
    end
    if (thetavec(end) < pi/6) || (thetavec(end) > 12*pi/13)
        dtheta =  dphi * sin(thetavec(end))
        if length(thetavec) == 2
            thetavec(2) = (thetavec(2) + dtheta)/2;
            dtheta = (thetavec(2) + dtheta)/2;
            dthetaPole = dtheta;
        end
    else
        dtheta = pi/(2*N);
    end
     thetavec = [thetavec,thetavec(end) + dtheta];
end
thetavec = [thetavec,pi];
%thetavec = [0,thetavec,2*pi];
[th, ph] = meshgrid(thetavec,phivec);
R = ones(size(th)); % should be your R(theta,phi) surface in general

x = R.*sin(th).*cos(ph);
y = R.*sin(th).*sin(ph);
z = R.*cos(th);

figure;
surf(x,y,z);
axis vis3d
surf(x,y,z);
fvc = surf2patch(x,y,z,'triangle');

%[vertices, m, n] = unique(fvc.vertices,"row");
% removed = setdiff(1:size(n),m);
% 
% [top, top_idx] = ismember(vertices, [0,0,1], 'rows');
% top_idx = find(top_idx);
% [bottom, bottom_idx] = ismember(vertices, [0,0,-1], 'rows');
% bottom_idx = find(bottom_idx);
faces = fvc.faces;
vertices = [0,0,1;0,0,-1];
for old_vInd = 1:size(fvc.vertices,1)
    if norm(fvc.vertices(old_vInd,:) - [0,0,1]) < 1e-3
        faces(faces == old_vInd) = 1;
    elseif norm(fvc.vertices(old_vInd,:) - [0,0,-1]) < 1e-2
        faces(faces == old_vInd) = 2;
    elseif (abs(fvc.vertices(old_vInd,2)) < 5e-10) && (fvc.vertices(old_vInd,2) < 0)
        fvc.vertices(old_vInd,2) = inf;
        [right_vInd,~] = ismember(vertices, [fvc.vertices(old_vInd,1),0,fvc.vertices(old_vInd,3)], 'rows');
        vertices = [vertices; fvc.vertices(right_vInd,:)];
        faces(faces == old_vInd) = find(right_vInd);
    else
        vertices = [vertices; fvc.vertices(old_vInd,:)];
        new_vInd = size(vertices,1);
        faces(faces == old_vInd) = new_vInd; 
    end
end


[newfaces, m, n] = unique(faces,"row");

removeList = [];
for j = 1:size(newfaces,1)
    if length(newfaces(j,:)) ~= length(unique(newfaces(j,:)))
        removeList = [removeList,j];
    end
end
newfaces(removeList,:) = [];

%patch('Faces',newfaces,'Vertices',vertices);
plywrite("C:\Users\Kieran\Desktop\DDG_membrane\build\bin\Release\input-file\UVsphere.ply",newfaces,vertices);









