function DoSimulate(handles)
% 
% Called by GUI DiffuTrapClusters 
%
% it creates un space with a given number of binding sites (slots) in 1,2,4 or 7 clusters
% if slots are distributed followin an hexagonal grid, it needs a datahexgrid file to generate the simulation space.
%
%
% it can save SPT-like trajectories
% it contabilizes the number of slots occupied at each time point
% it contabilizes the number of bindings in real time (simulation time)
%
% 
% Marianne Renner 02/19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------
% for future use
%-------------------------------------------------------------------------

% this can be used to introduce confinement (i.e. picket and fences)
Pextra=0; %desactivated
Psyn=0;  %desactivated
distconfextra=50;
distconfsyn=10;
numbersyn(1)=0;
numbersyn(2)=0;

%third group possible (fixed)
nroslotoccupied=0; 

%-------------------------------------------------------------------------

current_dir=cd;
tours=str2num(get(handles.simurounds,'String'));

disp('Simulating diffusion and trapping in 2D')
disp(' ')

%check parameters

% SPTlike trajectories
SPTlike=get(handles.SPTradiobutton,'Value');
if SPTlike==1
    prompt = {'Pixel size (nm):','Time between frames (ms):','Localization accuracy (nm):'};
    num_lines= 1;
    dlg_title = 'SPT-like trajectories';
    def = {'167','75','20'}; % default values
    answer  = inputdlg(prompt,dlg_title,num_lines,def);
    exit=size(answer);
    if exit(1) == 0
        return; 
    end
    szpx=str2num(answer{1}); %pixel size
    realtill=str2num(answer{2}); %pixel size
    paccu=str2num(answer{3});
else
   realtill=75;           
   szpx=167;           
end

%simulation space & trajectories

%------------------------------------
till=1; %1 ms  simulation time period
%------------------------------------

sizesimpx=str2num(get(handles.sizesimspace,'String')); %size simulation space in µm
sizesimspace=sizesimpx*szpx; %
distsites=str2num(get(handles.distsites,'String'));
nroslots=str2num(get(handles.bsites,'String')); %number of slots

%number of clusters    
clusteroption=get(handles.nroclupopupmenu,'Value'); %only some values accepted
switch clusteroption
    case 1 
        nronanoclu=1;
    case 2 
        nronanoclu=2;
    case 3 
        nronanoclu=4;
    case 4 
        nronanoclu=7;
end

lengthtrc=str2num(get(handles.trajlength,'String')); % length trajectories in ms
lengthtrcframes=lengthtrc/realtill; % length trajectories in SPT frames

nrotrajectories(1)=str2num(get(handles.nromol1,'String')); % number of trajectories
nrotrajectories(2)=str2num(get(handles.nromol2,'String'));

diammol(1)=str2num(get(handles.sizemol1,'String')); %size molecules
diammol(2)=str2num(get(handles.sizemol2,'String'));

% diffusion and trapping 
diffcoef(1)=str2num(get(handles.Diff1,'String')); % D diffusing molecules 
diffcoef(2)=str2num(get(handles.Diff2,'String'));
Dslot=str2num(get(handles.Dslots,'String')); % D bound molecules

pfree(1)=str2num(get(handles.Pfree1,'String')); % probfree 
pfree(2)=str2num(get(handles.Pfree2,'String'));

probabilidadinterac(1)=str2num(get(handles.Pbind1,'String')); % Pbind 
probabilidadinterac(2)=str2num(get(handles.Pbind2,'String'));

numeroocup(1)=str2num(get(handles.bound1,'String')); % number bound molecules at t=0 
numeroocup(2)=str2num(get(handles.bound2,'String'));
slotsocupinicio=sum(numeroocup(:));       

% changes in binding
dochange=get(handles.changePradiobutton,'Value');
defaultprob=probabilidadinterac(1);

timeprob=str2num(get(handles.framePbind,'String'));
valuechangeprob=str2num(get(handles.newPbind,'String'));

timechangeprob(1)=timeprob;
timechangeprob(2)=valuechangeprob;
timechangeprob(3)=pfree(1);

% changes of molecule number
if size(nrotrajectories,2)>1
    if numeroocup(1)>nrotrajectories(1) 
        disp('Error: verify the number of trajectories')
    else
        if size(diffcoef,2)<2 || size(pfree,2)<2 || size(probabilidadinterac,2)<2
            disp('Error: verify the number of trajectories')
        end
    end
end

%data sites positions
hexadistribution=get(handles.hexaradiobutton,'Value');
if hexadistribution==1
    [filedatagrid,path] = uigetfile('*datahexgrid*','Load datahexgrid file');
    if isequal(filedatagrid,0)
        return
    end
    cd(path)
    [namefile, ~] = strtok(filedatagrid,'.');
    [~, remain] = strtok(namefile,'-');
    [stringsyndiam, remain2] = strtok(remain,'-');
    rem=remain2(2:size(remain2,2));
    dist=str2num(rem);
    syndiam=str2num(stringsyndiam);
    set(handles.syndiam,'String',syndiam);  %  # actualize window
    set(handles.distsites,'String',distsites);  
    set(handles.disthexatext,'String',filedatagrid);  

else
  filedatagrid=[];
  path=cd;
  syndiam=str2num(get(handles.syndiam,'String')); 
end

% initialize %%%%%%%%%%%%%change this, not needed any more
if size(nrotrajectories,2)>1 % twogroups
    groupos=1;
    nrotraj=nrotrajectories(1)+nrotrajectories(2);
else
    groupos=0;
    nrotraj=nrotrajectories;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

for seriescounter=1:tours

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng('shuffle') %re-initialize rnd machine

waitbarhandle=waitbar( 0,'Please wait...','Name',['Simulating synaptic trajectories (serie # ',num2str(seriescounter),' of ',num2str(tours),')']) ;
synapse=[];

% create syn --------------------------------------------------

if hexadistribution==1
    % binding sites distributed on an hexagonal grid (datahexgrid file)
    [synapse,~,matrixslots,~,mask,~]=synnanohexamask2(filedatagrid,nroslots,sizesimspace,nronanoclu);
else
    %binding sites distributed randomly
    [synapse,~,matrixslots,~,mask,~]=synnanorandmask2(distsites,nroslots,syndiam,sizesimspace,nronanoclu);
end

if isempty(synapse)==1
    return
end

if max(max(synapse))==0
    return
end
auxslots=matrixslots;

% occupied at the beginning (third group of molecules, immobile)
if nroslotoccupied>0
    if nroslotoccupied>nroslots
        nro=(nroslots)-nroslotoccupied;
        matrixslots(:,2)=1;
        matrixslots(:,3)=10000; % 
        matrixslots(:,4)=1; % frame=1
        matrixslots(:,5)=10000; % nrotraj
        matrixslots(:,6)=3; % group 3
        matrixslots(:,10)=0; % contabilize bindings to this slot
        indexno=zeros(nro,1);
        for cont=1:nro
            indexno(cont)=floor(rand*(nroslots))+1;
        end
        if isempty(indexno)==0
            matrixslots(indexno,2:6)=0;
        end
    else
        for cont=1:nroslotoccupied+1
            indexoui(cont)=floor(rand*(nroslots))+1; 
        end
        matrixslots(indexoui,2)=1;
        matrixslots(indexoui,3)=10000; % 
        matrixslots(indexoui,4)=1; % frame=1
        matrixslots(indexoui,5)=10000; % nrotraj
        matrixslots(indexoui,6)=3; % group 3
        matrixslots(:,10)=0; % contabilize bindings to this slot

    end
end
%-------------------------------------------------------------------------

pos{nrotraj,2}=[]; % positions all trajectories, all simulation frames 
compileduration=[]; % duration of each interaction + number slot
occupation=[]; % % occupied slots
newtrcdata=[]; %preallocate!!!! % trc files, in wanted frame rate
allnewtrcdata=[]; %preallocate!!!! % real trc files

matrixposmol=zeros(nrotraj,3);  % positions mol at each frame

%-------------------------------------------------------------------------

% first point 
nro=1;
firstpos=1;
prob=zeros(nro,1);
group=zeros(nrotraj,1);
countslots=zeros(nro,1);
D=zeros(nro,1);
sizeconf=zeros(nro,1);
first=zeros(nro,1);
probconf=zeros(nro,1); %for future use

while nro<nrotraj+1
    if exist('waitbarhandle') %#ok<EXIST>
        waitbar(1/lengthtrcframes,waitbarhandle,'Initializing...');
    end
    
    % chooses an initial position free of molecules -----------------------
    forbidden=1;
    while forbidden
           pos{nro}(1,1)=((rand-0.5)*sizesimspace)+sizesimspace/2; %initial pos in x at random
           pos{nro}(1,2)=((rand-0.5)*sizesimspace)+sizesimspace/2; %initial pos in y at random 
           matrixposmol(nro,1)=nro;
           matrixposmol(nro,2)=pos{nro}(1,1);
           matrixposmol(nro,3)=pos{nro}(1,2);
           if nro>1
               forbidden=checkconflict2(matrixposmol, nro, firstpos,diammol,group);
           else
               forbidden=0;
           end
    end
    
    %counting binding
    matrixbinding{nro}=[];
    
   % slots at the beginning?----------------------------------------
   
    if groupos==1  % groups?
        if nro<nrotrajectories(1)+1 %first group
            Dvalue=diffcoef(1);
            probfree=pfree(1);
            probinterac=probabilidadinterac(1); %#ok<NASGU>
            % slots or syn at the beginning?
            if numbersyn(1)>0
                if nro<numbersyn(1)+1 % first trajs 1st group on syn
                    pos=findsyn(synapse,indexpossyn,pos,nro);
                else
                    if numeroocup(1)>0
                        if nro<numeroocup(1)+numbersyn(1)+1 % first trajs 1st group on slots
                            % position of slot
                            [pos,auxslots]=findslot(nroslots,auxslots,pos,nro);
                            prob(nro)=1;
                        end
                    end
                end
            else
                if numeroocup(1)>0
                    if nro<numeroocup(1) % first trajs 1st group on slots
                        % position of slot
                        [pos,auxslots]=findslot(nroslots,auxslots,pos,nro);
                        prob(nro)=1;
                    end
                end
            end
            group(nro)=1;
        else
            Dvalue=diffcoef(2);
            probfree=pfree(2);
            probinterac=probabilidadinterac(2);
            % slots or syn at the beginning?
            if numbersyn(2)>0
                if nro<numbersyn(2)+nrotrajectories(1)+1 % first trajs 2nd group on syn
                    pos=findsyn(synapse,indexpossyn,pos,nro);
                else
                    if numeroocup(2)>0
                        if nro<numeroocup(2)+nrotrajectories(1)+numbersyn(2)+1 % first trajs 1st group on slots
                            % position of slot
                            [pos,auxslots]=findslot(nroslots,auxslots,pos,nro);
                            prob(nro)=1;
                        end
                    end
                end
            else
                if numeroocup(2)>0
                    if nro<numeroocup(2)+nrotrajectories(1) % first trajs 1st group on slots
                        % position of slot
                        [pos,auxslots]=findslot(nroslots,auxslots,pos,nro);
                        % duration!
                        prob(nro)=1;
                    end
                end
            end
            group(nro)=2;
        end
    else % no groups
        %group(nro)=1;
        Dvalue=diffcoef(1);
        probfree=pfree(1);
        probinterac=probabilidadinterac(1);
        if numeroocup(1)>0
            if nro<numeroocup(1) % first trajs on slots
                % position of slot
               [pos,auxslots]=findslot(nroslots,auxslots,pos,nro);
                % duration!
                prob(nro)=1;
            end
        end
        group(nro)=1;
    end 
    countslots(nro)=0;
    control=1  ;
    
    % syn at the beginning? -----------------------------------------------------------
    
    if synapse(ceil(pos{nro}(1,2)),ceil(pos{nro}(1,1)))==0 %transposed!!!!!!!!
        %extra
        D(nro)=Dvalue;        
        prob(nro)=Pextra; sizeconf(nro)=distconfextra;
    else
        if synapse(round(pos{nro}(1,2)),round(pos{nro}(1,1)))==1000 %syn not slot
            D(nro)=Dvalue;
            prob(nro)=Psyn; sizeconf(nro)=distconfsyn;
        elseif synapse(round(pos{nro}(1,2)),round(pos{nro}(1,1)))==2 %slot area
            D(nro)=Dvalue;
            prob(nro)=0.99; sizeconf(nro)=5; % for future use: differences inwithin clusters
        else %slot or around slot
            nroslot=abs(synapse(round(pos{nro}(1,2)),round(pos{nro}(1,1))));
            indexslot=find(matrixslots(:,1)==nroslot);
            % slot occupied?
            if matrixslots(indexslot,2)==0 %free slot : it occupies it
                D(nro)=Dslot;
                center{nro}=pos{nro}(1,:);
                prob(nro)=1; sizeconf(nro)=5; %even if it diffuses, remains confined around site
                countslots(nro)=countslots(nro)+1;
                matrixslots(indexslot,1)=nroslot; % nroslot
                matrixslots(indexslot,2)=1; % interacting
                matrixslots(indexslot,3)=probfree; % 
                matrixslots(indexslot,4)=1; % frame=1
                matrixslots(indexslot,5)=nro; % nrotraj
                matrixslots(indexslot,6)=group(nro); % groups
                matrixslots(indexslot,10)=0; % first one not taken into account
                matrixbinding{nro}=[1 1 group(nro)];
                
            elseif matrixslots(indexslot,1)>0
                % forbidden
                control=0;
                D(nro)=Dvalue;
                prob(nro)=Psyn; 
                sizeconf(nro)=distconfsyn;
            end
        end %elseif
    end % synapse  
    
    if control==1    
        loc{nro}(1)=synapse(ceil(pos{nro}(1,2)),ceil(pos{nro}(1,2))); 
        first(nro)=1;
        probconf(nro)=0;
        newtrcdata=[newtrcdata; nro 1 pos{nro}(1,1) pos{nro}(1,2) D(nro) synapse(ceil(pos{nro}(1,2)),ceil(pos{nro}(1,1))) mask(ceil(pos{nro}(1,2)),ceil(pos{nro}(1,1))) group(nro)];
        nro=nro+1;
    end
    
end % loop trajectories first point
%-----------------------------------------------------------------------
% loop frames
%-----------------------------------------------------------------------
firstpos=0;

if exist('waitbarhandle') %#ok<EXIST>
    waitbar(2/lengthtrcframes,waitbarhandle,'Frame # 1');
end

nexttill=realtill/till+1; % first point
counttill=2;

for i=2: lengthtrc % loop all frames from frame=2
    
     %change of probability
     if dochange>0
         if i>timechangeprob(1)
             probabilidadinterac(1)=timechangeprob(2);
         else
             probabilidadinterac(1)=defaultprob; % when several series!
         end
     end

    %  all trajectories at frame=i
    for nro=1:nrotraj             
    
        % groups?     
        if group(nro)==1
                Dvalue=diffcoef(1);
                probfree=pfree(1);
                probinterac=probabilidadinterac(1);
        elseif group(nro)==2
                Dvalue=diffcoef(2);
                probfree=pfree(2);
                probinterac=probabilidadinterac(2);
        else
            Dvalue=diffcoef(1);        
            probfree=pfree(1);
            probinterac=probabilidadinterac(1);
        end
        
        %binding status
        databinding=matrixbinding{nro};
        
        % for future use
        if probconf(nro)==prob(nro) % same conf than before
        else
            if prob(nro)==1
                first(nro)=1; %start of strong confinement
            else
                first(nro)=0; % not strong confinement
            end
        end
        probconf(nro)=prob(nro) ;
        
        % displacements - pixel simulation=1 nm
        % %alternative way: 
        %ang=rand*360/57.29; %radians
        %Displ=abs(normrnd(0,sqrt(2*D(nro)*(till/1000)/((1/1000)^2)))); % pixel/frame   
        
         % displacements - pixel simulation = 1 nm
        Displx=normrnd(0,sqrt(2*D(nro)*(till/1000)/((1/1000)^2))); %from µ2/m to pixel/time point
        Disply=normrnd(0,sqrt(2*D(nro)*(till/1000)/((1/1000)^2))); 
        Displ=sqrt(Displx^2+Disply^2); %step distance
        ang=atan(Disply/Displx); %angle (to check for sites in the trajectory)

        % for future use: confinement
        conf=rand; %random probability  
        if conf<probconf(nro)  
            if probconf(nro)==1 % immo or ++ conf: stop periods
                if first(nro)==1
                    center{nro}=pos{nro}(1,:); %#ok<*AGROW> % center of the confinement area
                    first(nro)=0;
                end
            else % confinement soft: obstacles
                center{nro}=[];
            end
            if Displ>sizeconf(nro) % displacement in x outside of the confinement area/obstacle
                veces=floor(Displ/sizeconf(nro));
                Displ=sizeconf(nro)-(Displ-veces*sizeconf(nro));
                controldepas=1;
            else
                controldepas=0;
            end
        else
            first(nro)=1;        
            center{nro}=[];
        end %conf   
        
        % new positions
        if isempty(center{nro})==0 && controldepas==1  % confined
          %  newposx=center{nro}(1)+Displ*sin(ang); %alternative way
          %  newposy=center{nro}(2)+Displ*cos(ang);
            newposx=center{nro}(1)+Displx;
            newposy=center{nro}(2)+Disply;
        else  % not confined
          %  newposx=pos{nro}(1,1)+Displ*sin(ang); %alternative way
          %  newposy=pos{nro}(1,2)+Displ*cos(ang);
            newposx=pos{nro}(1,1)+Displx;
            newposy=pos{nro}(1,2)+Disply;
        end
        matrixposmol(nro,2)=newposx;
        matrixposmol(nro,3)=newposy;
        
        %check next point
        [possynx,possyny,newposx, newposy]=defineposmol2(newposx, newposy, sizesimspace);
        
        %----------------------------------------------------------------
        % check forbidden: positions around    
        forbidden=checkconflict2(matrixposmol, nro, firstpos,diammol,group);
        
        if forbidden==1
                % new position
                if D(nro)==Dslot % does not move
                  if size(center,2)>nro
                      newposx=center{nro}(1);
                      newposy=center{nro}(2);
                  else
                      newdispl=Displ-(Displ/10); 
                      if isempty(center{nro})==0 && controldepas==1  % confined
                          newposx=center{nro}(1)+newdispl*sin(ang);
                          newposy=center{nro}(2)+newdispl*cos(ang);
                      else  % not confined
                          newposx=pos{nro}(1,1)+newdispl*sin(ang);
                          newposy=pos{nro}(1,2)+newdispl*cos(ang);
                      end                      
                  end
                      
                else
                    newdispl=Displ-(Displ/2); 
                    if isempty(center{nro})==0 && controldepas==1  % confined
                        newposx=center{nro}(1)+newdispl*sin(ang);
                        newposy=center{nro}(2)+newdispl*cos(ang);
                    else  % not confined
                        newposx=pos{nro}(1,1)+newdispl*sin(ang);
                        newposy=pos{nro}(1,2)+newdispl*cos(ang);
                    end
                end
                [possynx,possyny,newposx, newposy]=defineposmol2(newposx, newposy, sizesimspace);
                
                matrixposmol(nro,2)=newposx;
                matrixposmol(nro,3)=newposy;
                
                 % control again
                forbidden=checkconflict2(matrixposmol, nro, firstpos,diammol,group);
                if forbidden==1                    
                        % cannot move
                        newposx=floor(pos{nro}(1,1));
                        newposy=floor(pos{nro}(1,2)); 
                        [possynx,possyny,newposx, newposy]=defineposmol2(newposx, newposy, sizesimspace);
                        matrixposmol(nro,2)=newposx;
                        matrixposmol(nro,3)=newposy;
                end % forbidden 2
         end % forbidden 1
         
         % error control: impossible to move
         if isnan(possynx)
             possynx=pos{nro}(1,1);
         end
         if isnan(possyny)
             possynx=pos{nro}(1,2);
         end
        %--------------------------------------------------------------
        
        % check syn         
       
        indextraj=find(matrixslots(:,5)==nro);
        
        if isempty(indextraj)==0  % already on slot
            free=rand; 
            
            if probfree>free %gets free
                compileduration=[compileduration; nro i matrixslots(indextraj,4) matrixslots(indextraj,1) matrixslots(indextraj,6)];
                matrixslots(indextraj,2)=0;
                matrixslots(indextraj,3)=0; 
                matrixslots(indextraj,4)=0; 
                matrixslots(indextraj,5)=0;
                matrixslots(indextraj,6)=0; % groups
                D(nro)=Dvalue;
                prob(nro)=Psyn; sizeconf(nro)=distconfsyn;
                
            else %still interacting 
                matrixslots(indextraj,4)=matrixslots(indextraj,4)+1; % one more frame
                D(nro)=Dslot;
                prob(nro)=1;
                
                if isempty(databinding)==0
                    databinding(size(databinding,1),2)= databinding(size(databinding,1),2)+1;
                end

            end
            
        else % indextraj; not stabilized 
            
            if synapse(possynx,possyny)==0
                %extra
                D(nro)=Dvalue;
                prob(nro)=Pextra; sizeconf(nro)=distconfextra;
            else
                % it was outside slot
                if synapse(possynx,possyny)==1000 %syn not slot
                    D(nro)=Dvalue;
                    prob(nro)=Psyn; sizeconf(nro)=distconfsyn;

                    %verify during the path (1nm)
                    pas=1; %1 pixel... 
                    while pas<Displ
                        xpas=round((pas*cos(ang))+pos{nro}(1,1));
                        ypas=round((pas*sin(ang))+pos{nro}(1,2));   
                        if synapse(ypas,xpas)<0 % slot found on the way
                            % probability of interaction
                            free=rand; 
                            if probinterac>free
                                nroslot=abs(synapse(possynx,possyny));
                                indexslot=find(matrixslots(:,1)==nroslot);
                                if nroslot<1000 && matrixslots(indexslot,2)==0 % slot free: it binds
                                    D(nro)=Dslot;
                                    prob(nro)=1; 
                                    sizeconf(nro)=5;
                                    newposx=xpas;
                                    newposy=ypas;
                                    possynx=round(ypas);
                                    possyny=round(xpas);
                                    center{nro}=[xpas ypas];   
                                    countslots(nro)=countslots(nro)+1;
                                    matrixslots(indexslot,1)=nroslot; % nroslot
                                    matrixslots(indexslot,2)=1; % interacting
                                    matrixslots(indexslot,3)=probfree; 
                                    matrixslots(indexslot,4)=1; % frame=1
                                    matrixslots(indexslot,5)=nro; % nrotraj    
                                    matrixslots(indexslot,6)=group(nro); % groups  
                                    matrixslots(indexslot,10)=matrixslots(indexslot,10)+1;
                                    databinding=[databinding; i i group(nro)]; %new binding
                                    
                                elseif matrixslots(indexslot,2)>0 
                                    if nroslot<1000 && matrixslots(indexslot,5)~=nro % other occupant
                                        %obstacle                    
                                        newposx=pos{nro}(1,1); 
                                        newposy=pos{nro}(1,2);
                                        pas=Displ; % obstacle=stop
                                    end
                                end
                            else % it does not interact
                                D(nro)=Dvalue;
                                prob(nro)=0.99; sizeconf(nro)=5;
                            end % probinteract
                            pas=Displ;
                        else % synapse
                            pas=pas+sqrt(2);
                        end %synapse
                    end % while
                
                elseif synapse(possynx,possyny)<0  % it reached one slot
                    
                    % probability of interaction
                    free=rand; 
                    if probinterac>free %has to interact
                        % slot occupied?     
                        nroslot=abs(synapse(possynx,possyny))  ;
                        indexslot=find(matrixslots(:,1)==nroslot);
                        
                        if matrixslots(indexslot,2)==0 %free slot: it binds
                            
                            center{nro}=[newposx newposy];
                            D(nro)=Dslot;
                            prob(nro)=1;                
                            countslots(nro)=countslots(nro)+1;
                            matrixslots(indexslot,1)=nroslot; % nroslot
                            matrixslots(indexslot,2)=1; % interacting
                            matrixslots(indexslot,3)=probfree;
                            matrixslots(indexslot,4)=1; % frame=1
                            matrixslots(indexslot,5)=nro; % nrotraj
                            matrixslots(indexslot,6)=group(nro); % groups 
                            matrixslots(indexslot,10)=matrixslots(indexslot,10)+1;
                            databinding=[databinding; i i group(nro)]; %new binding
                            
                        else %indexslot
                        end %index
                    else %probinterac
                        D(nro)=Dvalue;
                        prob(nro)=Psyn; sizeconf(nro)=distconfsyn;
                    end %probinterac
                end %synapse
            end %on slot?
        end %indextraj
        
        pos{nro}(1,1)=newposx;
        pos{nro}(1,2)=newposy;
        
        matrixposmol(nro,2)=newposx;
        matrixposmol(nro,3)=newposy;

        if i==nexttill % SPT trajectory
           newtrcdata=[newtrcdata; nro counttill newposx newposy D(nro) synapse(possynx,possyny) mask(possynx,possyny) group(nro)];
        end
        
        matrixbinding{nro}=databinding;
        maxnro=nro;
    end % trajectories
    
    if exist('waitbarhandle') %
        waitbar(i/lengthtrc,waitbarhandle,['Point # ',num2str(i)]);
    end

    if i==nexttill % SPT trajectory
        counttill=counttill+1;
        nexttill=counttill*realtill/till; % VER
    end
    %---------------------------------------------------------------------
    
end% points trajectories
if exist('waitbarhandle') %#ok<EXIST>
   close(waitbarhandle);
end

%--------------------------------------------------------------------------

%count binding
resbinding=[];
for i=1:maxnro
    data=matrixbinding{i};
    if isempty(data)==0
        resbinding=[resbinding; i * ones(size(data,1),1) data ];
    end
end

%--------------------------------------------------------------------------
%count occupation slots
newtrcdata=sortrows(newtrcdata,1);
syntrc=[];
for gg=1:max(newtrcdata(:,1))
    index=find(newtrcdata(:,1)==gg);
    if isempty(index)==0
        data=newtrcdata(index,:);
        for hh=1:size(data,1)
            if data(hh,5)==Dslot
               data(hh,6)=10;
            elseif data(hh,6)<0
               data(hh,6)=1000;
           end
       end
       syntrc=[syntrc;data];
    end
end
loctrc=newtrcdata;
loctrc(:,6)=loctrc(:,7);
for i=1:max(syntrc(:,2)) %all frames
    index=find(syntrc(:,2)==i);
    trcframe=syntrc(index,:);  
    indexsyn=find(trcframe(:,6)==10);
    if isempty(indexsyn)==0
        index1=find(trcframe(indexsyn,8)==1);
        index2=find(trcframe(indexsyn,8)==2);
        if isempty(index1)==1
            count1=0;
        else
            count1=size(index1,1);
        end
        if isempty(index2)==1
            count2=0;
        else
            count2=size(index2,1);
        end
    else
        count1=0;
        count2=0;
    end
    occupation=[occupation; i*(realtill/1000) count1 count2];            
end % all frames

%save folder and name
if isdir('simsyn\trc'); else mkdir('simsyn\trc'); end
cd('simsyn')
d=dir('reportsim*');
st = {d.name};
nroorder=size(st,2)+1;
savename = ['sim-',num2str(nroorder)] ;
save([savename,'-slots.txt'],'occupation','-ascii');

%--------------------------------------------------------------------------
% SPT-like trajectories
%--------------------------------------------------------------------------

if SPTlike==1
    % noise
    for t=1:size(newtrcdata,1)
        noisex=normrnd(0,paccu/2);
        noisey=normrnd(0,paccu/2);   
        newtrcdata(t,3)=newtrcdata(t,3)+noisex ;
        newtrcdata(t,4)=newtrcdata(t,4)+noisey ;
    end
    % convertion to desired pixel size
    for t=1:size(newtrcdata,1)
        newtrcdata(t,3)=newtrcdata(t,3)/szpx ;
        newtrcdata(t,4)=newtrcdata(t,4)/szpx ;
    end
    trcdata=newtrcdata;
    newtrcdata(:,6)=newtrcdata(:,7);
    
    % create and save trajectories files
    cd('trc')
    save([savename,'.syn.trc'],'syntrc','-ascii'); % with data slots
    cd(path)
    cd('simsyn')
end %if SPTlike==1

%-----------------------------------------------------------------------
% save image simulation space
synapsegray = imadd(synapse,1); % white
synapsegray(round(sizesimspace/2),round(sizesimspace/2))=0; % black point for displaying reasons
imwrite(synapsegray,[savename,'-sites.tif'],'compression','none');
clear synapsegray

%-------------------------------------------------------------------------
%report simulation 
name=['report',savename,'.txt'];
fi = fopen(name,'wt');
if fi<3
    error('File not found or readerror.');
end

report{1}= ['Simulation number ',num2str(nroorder)];
report{2}= '  ';
report{3}= ['Simulation space : ',num2str(sizesimpx),' µm. Synapse diameter: ',num2str(syndiam),' nm.'];
report{4}= ['Trajectory length (points): ',num2str(lengthtrc)];

if hexadistribution==1
    report{5}= ['Number of binding sites: ',num2str(nroslots),' in ',num2str(nronanoclu),' cluster(s). File grid: ',filedatagrid];
else
    report{5}= ['Number of binding sites: ',num2str(nroslots),' in ',num2str(nronanoclu),' cluster(s).'];
end

report{6}= '  ';
report{7}= 'Group 1:';
report{8}= ['Size of molecules: ',num2str(diammol(1)),'  nm. D: ',num2str(diffcoef(1)),'(µm2/s). ',num2str(numeroocup(1)),' molecules bound at t=0.'];
report{9}= ['Probability of binding :', num2str(defaultprob),'. Probability to get free: ',num2str(pfree(1)),'.' ];
if dochange>0
    report{10}= ['Probability of binding changes at ', num2str(timechangeprob(1)),' frames, to be ', num2str(valuechangeprob)];
end

if nrotrajectories(2)>0  % two groups
    report{11}= 'Group 2:';
    report{12}= ['Size of molecules: ',num2str(diammol(2)),'  nm. D: ',num2str(diffcoef(2)),'(µm2/s). ',num2str(numeroocup(2)),' molecules bound at t=0.'];
    report{13}= ['Probability of binding :', num2str(probabilidadinterac(2)),'. Probability to get free: ',num2str(pfree(2)),'.'];

     report{14}= '  ';
    if SPTlike==1
        report{15}= ['SPT-like trajectories saved with integration time: ',num2str(realtill),' ms, pixel size: ',num2str(szpx),' nm   Loc accuracy : ', num2str(paccu),'  nm'];
    else
        report{15}=' ';
    end

else
    if SPTlike==1
        report{11}= ['SPT-like trajectories saved with integration time: ',num2str(realtill),' ms, pixel size: ',num2str(szpx),' nm   Loc accuracy : ', num2str(paccu),'  nm'];
    else
        report{11}=' ';
    end
end
    
if nrotrajectories(2)>0 % two groups 
    for celda=1:15
        fseek(fi,10,0);
        fprintf(fi,'%10s\n',report{celda});
    end
else
    for celda=1:11
        fseek(fi,10,0);
        fprintf(fi,'%10s\n',report{celda});
    end
end
    
fclose(fi);
disp(['Simulation #',num2str(nroorder),' saved'])

clear trcdata newtrc newtrcdata allnewtrcdata synapse newsynapse newsynapse2 report occupation %VERIF

cd(path);

end % tours

disp('Done')
disp(' ')

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pos,auxslots]=findslot(nroslots,matrixslots,pos,nro);

aux=matrixslots;
indexnotzero=find(aux(:,1)>0);
matrixslotsclean=matrixslots(indexnotzero,:);

positionclean=floor(rand*size(matrixslotsclean,1));
if positionclean<1
    positionclean=1;
end
nroslot=matrixslotsclean(positionclean,1);
position=find(matrixslots(:,1)==nroslot);
   
pos{nro}(1,1)=matrixslots(position,8); %initial pos in slot
pos{nro}(1,2)=matrixslots(position,7); %

auxslots=matrixslots(find(matrixslots(:,1)~=nroslot),:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% eof
