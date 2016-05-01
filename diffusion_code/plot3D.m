function pts = plot3D(bdryPts,ctrlPts,dt,moviename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plot3D(bdryPts,ctrlPts,dt,moviename)
%
% INPUT:
%   bdryPts
%   ctrlPts
%   dt
%   moviename - if moviename is not provided no movie is saved
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




N = size(bdryPts,3);
nofBdryPts = size(bdryPts,1);
nofCtrlPts = size(ctrlPts,1);
c = jet(N)

figure(1)
hold off


% variable to determine whether to create a movie
create_movie = false;
if (nargin>3)
    create_movie = true;
end


if create_movie
  mov(1:N) = struct('cdata', [],'colormap', []);
end



T = [0:dt:N/dt];
pts = zeros(N*nofBdryPts,3);
for i=1:N
    pts((i-1)*nofBdryPts+1:i*nofBdryPts,:) = [bdryPts(:,:,i) dt*i*ones(nofBdryPts,1)];
    temp = [bdryPts(:,:,i) dt*i*ones(nofBdryPts,1);bdryPts(1,:,i) dt*i];
    temp1 = [ctrlPts(:,:,i) dt*i*ones(nofCtrlPts,1);ctrlPts(1,:,i) dt*i];


    view(45,30)
    


    h1 = plot3(temp(:,3),temp(:,1),temp(:,2),'Color','k','Linewidth',3);
    view(45,30)
    a = axis;
    a(2) = 40;
    % axis([0    40    20    70    35    55]);
    axis([0 40 10 50 10 50])
    axis([0 40 10 70 10 70])
    % axis equal
    % axis([0 40 min(min(bdryPts(:,1,:)))-5 max(max(bdryPts(:,1,:)))+5 min(min(bdryPts(:,2,:)))-5 max(max(bdryPts(:,2,:)))+5])
    % axis(a)

    hold on
    h2 = plot3(temp1(:,3),temp1(:,1),temp1(:,2),'o','Linewidth',2,'MarkerFaceColor','k','MarkerEdgeColor','k');
 
    set(gca,'ytick',[])
    set(gca,'ztick',[])



    pause(0.01)
    if create_movie
        mov(i) = getframe(gcf);
    end 
  
    set(h1,'Color',c(i,:), 'Linewidth',2);
    set(h2,'Color','k');

 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% optionally plot the ellipses for multiple_OU
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if false
set(h1,'Color','k', 'Linewidth',2);
a = [0:0.1:2*pi+0.1];
ellipse_lr = [10*sin(a)' 5*cos(a)'] + 30;
ellipse_ud = [5*sin(a)' 10*cos(a)'] + 30;
plot3(30.1*ones(64,1),1.5*ellipse_ud(:,1),1.5*ellipse_ud(:,2),'Linewidth',2)
plot3(30.1*ones(64,1),1.5*ellipse_lr(:,1),1.5*ellipse_lr(:,2),'Linewidth',2)
end

mov(i) = getframe(gcf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% optionally plot the dumbbell for OU
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set(h1,'Color','k', 'Linewidth',2);
%load('matlab.mat','dumbbell')
%plot3(30.1*ones(82,1),1.5*dumbbell(:,1),1.5*dumbbell(:,2),'Linewidth',2)

%mov(i) = getframe(gcf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if create_movie
    movie2avi(mov, moviename, 'compression', 'None','FPS',120,'quality',100);
end
