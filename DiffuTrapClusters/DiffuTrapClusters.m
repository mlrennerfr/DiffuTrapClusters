function varargout = DiffuTrapClusters(varargin)
% DIFFUTRAPCLUSTERS MATLAB code for DiffuTrapClusters.fig
%
% GUI to launch simulations of lateral diffusion and trappnig on 2D.
%
% Marianne Renner 2020
% Last Modified by GUIDE v2.5 21-Jun-2020 15:25:02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DiffuTrapClusters_OpeningFcn, ...
                   'gui_OutputFcn',  @DiffuTrapClusters_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DiffuTrapClusters_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

guidata(hObject, handles);

%------------------------------------------------------------------------
function varargout = DiffuTrapClusters_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gopushbutton_Callback(hObject, eventdata, handles)

DoSimulate(handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gridpushbutton_Callback(hObject, eventdata, handles)
% create datahexagrid file with the possible positions of sites, following
% an hexagonal grid

disp('Generation of a circular hexagonal grid of binding sites')

% dialog box
prompt = {'Diameter of the cercle (nm):','Desired distance between sites (nm):'};
num_lines= 1;
dlg_title = 'Circular hexagonal grid';
def = {'300','15'}; % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
exit=size(answer);
if exit(1) == 0
    return; 
end
sizesyn=str2num(answer{1});
h_goal=str2num(answer{2}); 

disp(' ')
disp(['Calculating grid of ',num2str(sizesyn),' nm in diameter; sites separated ',num2str(h_goal),' nm.'])

pcercle=calculhexgridcerclesim(sizesyn,h_goal);

disp(['Total number of possible sites: ',num2str(size(pcercle,1))])
disp(' ')
    
save(['datahexgrid-',num2str(sizesyn),'-',num2str(h_goal),'.txt'],'pcercle','-ascii')

disp(['File datahexgrid-',num2str(sizesyn),'-',num2str(h_goal),'.txt saved'])
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GUI 
%------------------------------------------------------------------------

function nromol1_Callback(hObject, eventdata, handles)
function nromol1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function nromol2_Callback(hObject, eventdata, handles)
function nromol2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Diff1_Callback(hObject, eventdata, handles)
function Diff1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Diff2_Callback(hObject, eventdata, handles)
function Diff2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Pbind1_Callback(hObject, eventdata, handles)
function Pbind1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Pbind2_Callback(hObject, eventdata, handles)
function Pbind2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Pfree1_Callback(hObject, eventdata, handles)
function Pfree1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Pfree2_Callback(hObject, eventdata, handles)
function Pfree2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bound1_Callback(hObject, eventdata, handles)
function bound1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bound2_Callback(hObject, eventdata, handles)
function bound2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sizesimspace_Callback(hObject, eventdata, handles)
function sizesimspace_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bsites_Callback(hObject, eventdata, handles)
function bsites_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function trajlength_Callback(hObject, eventdata, handles)
function trajlength_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sizemol1_Callback(hObject, eventdata, handles)
function sizemol1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sizemol2_Callback(hObject, eventdata, handles)
function sizemol2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function newPbind_Callback(hObject, eventdata, handles)
function newPbind_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function framePbind_Callback(hObject, eventdata, handles)
function framePbind_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function simurounds_Callback(hObject, eventdata, handles)
function simurounds_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function distsites_Callback(hObject, eventdata, handles)
function distsites_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Dslots_Callback(hObject, eventdata, handles)
function Dslots_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function syndiam_Callback(hObject, eventdata, handles)
function syndiam_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function changePradiobutton_Callback(hObject, eventdata, handles)

function hexaradiobutton_Callback(hObject, eventdata, handles)
hexaval=get(hObject,'Value') ;
if hexaval==0
    if isfield(handles,'disthexatest')
        set(handles.disthexatest,'String',' ');  
    end
    set(handles.syndiam,'Enable','on');  
else
    set(handles.syndiam,'Enable','off');  
end
guidata(hObject, handles);

function SPTradiobutton_Callback(hObject, eventdata, handles)

function nroclupopupmenu_Callback(hObject, eventdata, handles)
function nroclupopupmenu_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pcercle=calculhexgridcerclesim(sizesyn,h_goal)
% function function p_cercle=calculhexgridcerclesim(sizesyn,h_goal)
%
% sizesyn: size of cercle, in nm.
% h_goal: distance between slots
%
% pcercle: positions (in nm) of the slots taking into account hexagonal
% lattice, on a cercle of diameter sizesyn
%
% uses scripts of John Burkardt (2005)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

centro=sizesyn/2;
sizesyn=sizesyn+round(sizesyn/5);

box(:,1)=[1,1];
box(:,2)=[sizesyn,sizesyn]; % starts with a square

%nodes=round(sizesyn/h_goal);
nodes_per_layer=round(sizesyn/h_goal); % initial value
%[nodes_per_layer,n]=adjustnodes(n,nodes,box); % optimize nodes_per_layer to have n slots

if (nodes_per_layer < 1)
    layers = 0;
    hx = 0.0;
    hy = 0.0;
    n = 0;
elseif (nodes_per_layer == 1)
    layers = 1;
     hx = box(1,2) - box(1,1);
    hy = box(2,2) - box(2,1);
   n = 1;
else
    hx = (box(1,2) - box(1,1)) / (nodes_per_layer - 1);
    hy = sqrt (3) * hx/2;
    layers = 1 + floor (( box(2,2) - box(2,1)) / hy);
    n = nodes_per_layer* floor ((layers + 1) / 2) +   (nodes_per_layer - 1) * floor ((layers) / 2);
end

% optimize h
nodes_per_layer = 1 + round ( ( box(1,2) - box(1,1) ) / h_goal );

%  Check whether roundoff means we could use one less node per layer.
if ( 2 < nodes_per_layer )
    nodes_per_layer2 = nodes_per_layer - 1;
    h = ( box(1,2) - box(1,1) ) / ( nodes_per_layer2 - 1 );
    if ( h <= h_goal )
      nodes_per_layer = nodes_per_layer2;
      return
    end
end
h = ( box(1,2) - box(1,1) ) / ( nodes_per_layer - 1 );
n = nodes_per_layer* floor ((layers + 1) / 2) +   (nodes_per_layer - 1) * floor ((layers) / 2);

p = hex_grid_points ( nodes_per_layer, layers, n, box )';

% clean forbidden slots and write data
pnew=[];
countlayer=1;
paux=p(1,:);
pivot=1;
for i=2:size(p,1)
        if p(i,2)==p(i-1,2)
            paux=[paux;p(i,:)]; %even line
        else  % odd line
            if rem(countlayer,2)>0 
                pnew=[pnew; round(paux(:,:)) ];
            else
                pivot=-pivot;
                for j=1:size(paux,1)
                    if pivot>0
                        if rem(j,2)>0 
                            pnew=[pnew; round(paux(j,:)) ];
                        end
                    else
                        if rem(j,2)>0 
                        else
                            pnew=[pnew; round(paux(j,:)) ];
                        end
                    end
                end
            end
            paux=p(i,:);
            countlayer=countlayer+1;
        end
end

% cercle
pcercle=[];
for i=2:size(pnew,1)
    distcentro=sqrt((centro-pnew(i,1))^2+(centro-pnew(i,2))^2);
    if distcentro<centro
        pcercle=[pcercle;pnew(i,:)];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ nodes_per_layer, n ] = adjustnodes ( n_goal,nodes, box)
%function [ nodes_per_layer, n ] = adjustnodes ( n_goal,nodes, box)
% This routine experiments with various  values until it is convinced 
% it has the value of NODES_PER_LAYER   that comes as close as possible 
% to producing N nodes.
%
%  Modified from   John Burkardt
%    BOX(2,2), the lower and upper corners of the rectangular region.
%    N_GOAL, the desired number of nodes.
%    NODES_PER_LAYER, the number of nodes per layer
%    N, the number of nodes in the mesh.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  m = 2;
  nodes_per_layer_low = 0;
  n_low = 0;
  nodes_per_layer = nodes;
  nodes_per_layer_high = n_goal;
  n_high = n_goal * n_goal;

  while ( 1 )
      
    n = hex_grid_n ( nodes_per_layer, box );
    if ( n == n_goal )
      break
    end
    if ( n < n_goal )
      nodes_per_layer_low = nodes_per_layer;
      n_low = n;
    else
      nodes_per_layer_high = nodes_per_layer;
      n_high = n;
    end
    if ( nodes_per_layer_low + 1 == nodes_per_layer_high )
      if ( n - n_low <= n_high - n )
        nodes_per_layer = nodes_per_layer_high;
        n = n_high;
      else
        nodes_per_layer = nodes_per_layer_low;
        n = n_low;
      end
      break
    end
    nodes_per_layer = ...
      round ( ( nodes_per_layer_low + nodes_per_layer_high ) / 2 );
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = hex_grid_points ( nodes_per_layer, layers, n, box )
%*****************************************************************************80
%
%% HEX_GRID_POINTS returns coordinate box hex grid points.
%    This routine determines the coordinates of the elements of
%    a hexagonal grid in the unit square.
%    A hexagonal grid is defined in the coordinate box [A,B] x [C,D].
%    All nodes of the grid lie on one of LAYERS horizontal lines.
%    The first of these lines is from (A,C) to (B,C), and each
%    successive line is HY units higher.
%    On all the odd numbered lines, there are NODES_PER_LAYER points,
%    equally spaced from A to B, with a spacing of HX.
%    On the even numbered lines, there are NODES_PER_LAYER-1 points,
%    whose values are the midpoints of successive intervals on
%    an odd numbered line.  (The grid is staggered).
%    HY = HX * sqrt ( 3 ) / 2.
%
%  Licensing: This code is distributed under the GNU LGPL license.
%  Modified: 08 March 2005
%  Author:  John Burkardt
%  Parameters:
%    Input, integer NODES_PER_LAYER, the number of grid points on the first
%    horizontal layer of points.
%    Input, integer LAYERS, the number of horizontal layers.
%    Input, integer N, the total number of hex grid points.
%    Input, real BOX(2,2), the values of A, B, C and D
%    that define the coordinate box.
%    Output, real P(2,N), the coordinates of the
%    mesh points, listed one horizontal layer at a time.


  ndim = 2;
  if ( nodes_per_layer < 1 )
    p = [];
    return
  end
  if ( nodes_per_layer == 1 )
    p(1:ndim,1) = ( box(1:ndim,1) + box(1:ndim,2) ) / 2.0;
    return
  end

  [ hx, hy ] = hex_grid_h ( nodes_per_layer, box );
  k = 0;

  for j = 1 : layers
    y = box(2,1) + hy * ( j - 1 );
    jmod = mod ( j, 2 );

    if ( jmod == 1 )
      for i = 1 : nodes_per_layer
        x = box(1,1) + ( box(1,2) - box(1,1) ) * ( i - 1 ) / (nodes_per_layer - 1 );
        k = k + 1;
        if ( k <= n )
          p(1,k) = x;
          p(2,k) = y;
        end
      end
    else
      for i = 1 : nodes_per_layer-1
        x = box(1,1) + (box(1,2) - box(1,1)) * (2 * i - 1) / ( 2 * nodes_per_layer - 2 );
        k = k + 1;
        if ( k <= n )
          p(1,k) = x;
          p(2,k) = y;
        end
      end
    end
  end
  return
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ hx, hy ] = hex_grid_h ( nodes_per_layer, box )

%*****************************************************************************80
%
%% HEX_GRID_H computes the coordinate box hex grid spacings.
%    This routine determines the values of HX and HY from
%    the fundamental hexagonal grid parameter NODES_PER_LAYER.
%    A hexagonal grid is defined in the coordinate box [A,B] x [C,D].
%    All nodes of the grid lie on one of LAYERS horizontal lines.
%    The first of these lines is from (A,C) to (B,C), and each
%    successive line is HY units higher.
%    On all the odd numbered lines, there are NODES_PER_LAYER points,
%    equally spaced from A to B, with a spacing of HX.
%    On the even numbered lines, there are NODES_PER_LAYER-1 points,
%    whose values are the midpoints of successive intervals on
%    an odd numbered line.  (The grid is staggered).
%    HY = HX * sqrt ( 3 ) / 2.
%
%  Licensing:
%    This code is distributed under the GNU LGPL license.
%  Modified:
%    08 March 2005
%  Author:
%    John Burkardt
%  Parameters:
%    Input, integer NODES_PER_LAYER, the number of grid points on the first
%    horizontal layer of points.
%    Input, real BOX(2,2), the values of A, B, C and D
%    that define the coordinate box.
%    Output, real HX, the spacing between grid points
%    on a horizontal line.
%    Output, real HY, the spacing between horizontal lines.
%
  if ( nodes_per_layer < 1 )

    hx = 0.0;
    hy = 0.0;

  elseif ( nodes_per_layer == 1 )

    hx = box(1,2) - box(1,1);
    hy = box(2,2) - box(2,1);

  else

    hx = ( box(1,2) - box(1,1) ) / ( nodes_per_layer - 1 );
    hy = hx * sqrt ( 3.0 ) / 2.0;

  end 

  return
%end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
