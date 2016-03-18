function drawAircraft(uu,P)
    % process inputs to function
    pn       = uu(1);       % inertial North position     
    pe       = uu(2);       % inertial East position
    pd       = uu(3);           
    u        = uu(4);       
    v        = uu(5);       
    w        = uu(6);       
    phi      = uu(7);       % roll angle         
    theta    = uu(8);       % pitch angle     
    psi      = uu(9);       % yaw angle     
    p        = uu(10);       % roll rate
    q        = uu(11);       % pitch rate     
    r        = uu(12);       % yaw rate    
    t        = uu(13);       % time
    
    WL = 1000;

    % define persistent variables 
    persistent drone_handle;
    
    % first time function is called, initialize plot and persistent vars
    if t==0,
        figure(1), clf
        drone_handle = drawDroneModel(P.drawV,P.drawF,P.drawpatchcolors,...
                                               pn,pe,pd,phi,theta,psi,...
                                               [],'normal');
        title('Unmanned Aerial Vehicle')
        xlabel('East')
        ylabel('North')
        zlabel('-Down')
        view(32,47)  % set the vieew angle for figure
        %{
        set(gcf,'Renderer','opengl')
        hL1 = light('Position',[-WL,-WL,-WL],...
                    'Color', [1,1,1]);
        hL2 = light('Position',[WL,WL,WL],...
                    'Color', [1,1,1]);
        material shiny
        %}
        axis([-WL,WL,-WL,WL,-WL,WL])
        view(60,15)
        hold on
        
    % at every other time step, redraw base and rod
    else 
        drawDroneModel(P.drawV,P.drawF,P.drawpatchcolors,...
                           pn,pe,pd,phi,theta,psi,...
                           drone_handle);
    end
end

  
%=======================================================================
% drawDroneModel
% return handle if 3rd argument is empty, otherwise use 3rd arg as handle
%=======================================================================
%
function handle = drawDroneModel(V,F,patchcolors,...
                                     pn,pe,pd,phi,theta,psi,...
                                     handle,mode)
  V = rotate(V', phi, theta, psi)';  % rotate drone
  V = translate(V', pn, pe, pd)';  % translate drone
  % transform vertices from NED to XYZ (for matlab rendering)
  R = [...
      0, 1, 0;...
      1, 0, 0;...
      0, 0, -1;...
      ];
  V = V*R;
  
  if isempty(handle),
  handle = patch('Vertices', V, 'Faces', F,...
                 'FaceVertexCData',patchcolors,...
                 'FaceColor','flat',...
                 'EraseMode', mode,...
                 'EdgeColor', 'flat');%,...
                 %'FaceLighting', 'phong');
  else
    set(handle,'Vertices',V,'Faces',F);
    %shading interp
    drawnow
  end
end

%%%%%%%%%%%%%%%%%%%%%%%
function XYZ=rotate(XYZ,phi,theta,psi);
  % define rotation matrix
  R_roll = [...
          1, 0, 0;...
          0, cos(phi), -sin(phi);...
          0, sin(phi), cos(phi)];
  R_pitch = [...
          cos(theta), 0, sin(theta);...
          0, 1, 0;...
          -sin(theta), 0, cos(theta)];
  R_yaw = [...
          cos(psi), -sin(psi), 0;...
          sin(psi), cos(psi), 0;...
          0, 0, 1];
  R = R_yaw*R_pitch*R_roll;
  % rotate vertices
  XYZ = R*XYZ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% translate vertices by pn, pe, pd
function XYZ = translate(XYZ,pn,pe,pd)
  XYZ = XYZ + repmat([pn;pe;pd],1,size(XYZ,2));
end

  