function [newsynapse,indexsyn,matrixslots,listeobs,finalmask,synapse]=synnanorandmask2(distsites,nroslots,syndiam,sizesimspace,n)
% function [newsynapse,indexsyn,matrixslots,listeobs,finalmask,synapse]=synnanorandmask2(distsites,nroslots,syndiam,sizesimspace,n)
%
% creates a distribution of binding sites in 1,2,4 or 7 clusters
% sites are randomly distributed in clusters
%
% called by DoSimulate.m
% input:
%   file: datahexgrid
%   nroslots: number of binding istes (total)
%   sizesimspace: simulation area
%   n: number of clusters
% 
% Marianne Renner 01/20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

finalmask=[];
indexsyn=[];
listeobs=[];
newsynapse=zeros(sizesimspace,sizesimspace);
diam=syndiam+1 ;
centro=round(diam/2);
radio=floor(diam/2);
matrixslots=zeros(nroslots-1,12); 
synapse=zeros(diam+1,diam+1); 

resto=zeros(diam+1,ceil((sizesimspace-diam)/2));

for i=1:diam+1
    for j=1:diam+1
        distcentro=sqrt((centro-i)^2+(centro-j)^2);
        if distcentro<diam/2
            synapse(i,j)=1000;
        end
    end
end

synnanoclu=synapse; % mask to determine obstacle sites
poscentro=round(sizesimspace/2);
centersyn=[poscentro poscentro; poscentro poscentro]; 

if n==1 %one cluster
   % sidey=350;
    sidey=syndiam;
    sidex=sidey;
    distbetweenx=0; %nm entre nanoclu (1 seul)
    distbetweeny=0;
    poscentx(1)=centro; 
    poscenty(1)=centro;
    
elseif n==2
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %      x x
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if distsites==10
        sidey=150; %max distance in x or y
        sidex=sidey;
        poscentx(1)=-95+centro; 
        poscenty(1)=centro;
        poscentx(2)=95+centro;
        poscenty(2)=centro;

    elseif distsites==15
        
        sidey=230; %max distance in x or y
        sidex=sidey;
        poscentx(1)=-150+centro; 
        poscenty(1)=centro;
        poscentx(2)=150+centro;
        poscenty(2)=centro;
        
    else
        disp('distance between sites can only be 10 or 15 nm') 
    end

elseif n==4
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %      x x
        %      x x 
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if distsites==10 % ATT!!! only to organize nanolcusters in the same way than before 
        sidey=110; %max distance in x or y
        sidex=sidey;
        poscentx(1)=-80+centro; 
        poscenty(1)=-80+centro;
        poscentx(2)=80+centro;
        poscenty(2)=-80+centro;
        poscentx(3)=-80+centro;
        poscenty(3)=80+centro;
        poscentx(4)=80+centro;
        poscenty(4)=80+centro;
        
    elseif distsites==15 % ATT!!! only to organize nanolcusters in the same way than before 
        
        sidey=160; %max distance in x or y
        sidex=sidey;
        poscentx(1)=-120+centro; 
        poscenty(1)=-120+centro;
        poscentx(2)=120+centro;
        poscenty(2)=-120+centro;
        poscentx(3)=-120+centro;
        poscenty(3)=120+centro;
        poscentx(4)=120+centro;
        poscenty(4)=120+centro;
        
    else
        disp('distance between sites can only be 10 or 15 nm') 
    end

elseif n==7
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %      x x
    %     x x x 
    %      x x
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if distsites==10
        sidex=80;
        sidey=sidex;
        
        poscentx(1)=-50+centro; 
        poscenty(1)=-100+centro;
        
        poscentx(2)=50+centro;
        poscenty(2)=-100+centro;
        
        poscentx(3)=-100+centro;
        poscenty(3)=centro;
        
        poscentx(4)=centro;
        poscenty(4)=centro;
        
        poscentx(5)=100+centro; 
        poscenty(5)=centro;
        
        poscentx(6)=-50+centro;
        poscenty(6)=100+centro;
        
        poscentx(7)=50+centro;
        poscenty(7)=100+centro;
        
    elseif distsites==15
        
        sidex=120;
        sidey=sidex;
        
        poscentx(1)=-75+centro; 
        poscenty(1)=-150+centro;
        
        poscentx(2)=75+centro;
        poscenty(2)=-150+centro;
        
        poscentx(3)=-150+centro;
        poscenty(3)=centro;
        
        poscentx(4)=centro;
        poscenty(4)=centro;
        
        poscentx(5)=150+centro; 
        poscenty(5)=centro;
        
        poscentx(6)=-75+centro;
        poscenty(6)=150+centro;
        
        poscentx(7)=75+centro;
        poscenty(7)=150+centro;

    else
        disp('distance between sites can only be 10 or 15 nm') 
    end

end

totnumbernano=ceil(nroslots/n);
incrementnro=ceil(totnumbernano/100)*100;

numbersites=0;
maxsites=[]; %max number of sites (for area calculation)
        
for i=1:n % all nanoclusters

    countslots=incrementnro*(i-1)+1;

    % position of nanocluster's center fixed
    centx=poscentx(i);
    centy=poscenty(i);
    lastslot=nroslots*(i-1);
    if numbersites+totnumbernano>nroslots % too many...
        totnumbernano=nroslots-numbersites;
    end
    count=1;
    iter=1;
    
    while count<totnumbernano+1
        
        % random distribution, sites separated 10 or 15 nm
        % determine possible area (cercle around center)
        % utiliser rand pos generation as for initial positions, check
        % conflict (less than 10 or 15 nm)
        
         % chooses an initial position free of molecules -----------------------
         forbidden=1;
         iter=0;
         control=1;
         while forbidden
             % random angle and distance from center
              angles=rand*360/57.29;
              distx=(rand*sidex/2)*sin(angles);
              disty=(rand*sidey/2)*cos(angles);
              
              posi=round(distx+centx);
              posj=round(disty+centy);
              
             if count>1
                 forbidden=checkconflictsites(matrixslots, count, distsites,posi,posj);
             else
                 forbidden=0;
             end
             
             iter=iter+1;
             if iter>500000 %!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 disp('not succeded')
                 control=0;
                 break
             end
         end

         if control==1
             synapse(posi-1:posi+1,posj-1:posj+1)=-countslots; %forbidden area
             matrixslots(count+lastslot,1)=countslots;
             matrixslots(count+lastslot,11)=posi;
             matrixslots(count+lastslot,12)=posj;
             count=count+1;
             countslots=countslots+1;
             numbersites=numbersites+1;
         end
      
       maxsites=[maxsites; posi posj]; %%%%%%%%%%

       iter=iter+1;
       if iter>500000
          disp('not succeded')
           return
       end
    end
    
    t=1;
    m=1;
    startiii=centersyn(m,1)-radio;
    stopiii=centersyn(m,1)+radio;
    startj=centersyn(m,2)-radio;
    stopj=centersyn(m,2)+radio;
    
    for iii=startiii : stopiii
        k=1;
        for j=startj:stopj
            
            if synapse(t,k)<0 && synapse(t,k)>-5000 %obstacle sites not count in matrixslots!
                indexslot=find(matrixslots(:,1)==abs(synapse(t,k)));
                matrixslots(indexslot,7)=iii;
                matrixslots(indexslot,8)=j;
                matrixslots(indexslot,9)=i; %nro cluster
            end
            k=k+1;
        end
        t=t+1;
    end

end % loop nanoclusters

% figure
% imshow(synnanoclu,'InitialMagnification','fit')
   
aux=[resto synapse resto];
aux2=zeros(size(resto,2),size(aux,2));
newsynapse=[aux2;aux;aux2];

matrixslots=matrixslots(find(matrixslots(:,1)~=0),:);
  
%mask
[k,maxarea] = convhull(maxsites(:,1),maxsites(:,2)); %area calculated with ALL posible sites (not only those occupied)

mask=zeros(size(synapse,1),size(synapse,2));

mask = roipoly(mask, maxsites(k,1),maxsites(k,2)) ;
aux=[resto mask resto];
aux2=zeros(size(resto,2),size(aux,2));
finalmask=[aux2;aux;aux2];
finalmask=newsynapse;

indexsyn=[];

for i=1:size(newsynapse,1)
    index=find(newsynapse(i,:)>999);
    if isempty(index)==0
       indexsyn=[indexsyn; i*ones(size(index,2),1) index'];
    end
end

clear  synapse2 numberpergroup

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%