function [newsynapse,indexsyn,matrixslots,listeobs,finalmask,synapse]=synnanohexamask2(file,nroslots,sizesimspace,n)
%function [newsynapse,indexsyn,matrixslots,listeobs,finalmask,synapse]=synnanohexamask2(file,nroslots,sizesimspace,n)
%
% creates a distribution of binding sites in 1,2,4 or 7 clusters, based on
% hexagonal grid
% called by DoSimulate.m
% input:
%   file: datahexgrid
%   nroslots: number of binding istes (total)
%   sizesimspace: simulation area
%   n: number of clusters
% 
% Marianne Renner 01/20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize variables
newsynapse=[];
indexsyn=[];
matrixslots=[];
listeobs=[];
finalmask=[];
synapse=[];

% data grid
[namefile, ~] = strtok(file,'.');
[~, remain] = strtok(namefile,'-');
[stringsyndiam, remain2] = strtok(remain,'-');
rem=remain2(2:size(remain2,2));
dist=str2num(rem);
diamarea=str2num(stringsyndiam); % ok for one cluster

if n==1
    syndiam=ceil(diamarea*1.5);
elseif n==2
    syndiam=ceil(diamarea * 2.5);
elseif n==4
    syndiam=ceil(diamarea * 2.5);
elseif n==7
    syndiam=ceil(diamarea * 4);
end

datasqr=load(file);

newsynapse=zeros(sizesimspace,sizesimspace);
diam=syndiam+1 ;
centro=round(diam/2);
radio=floor(diam/2);
matrixslots=zeros(nroslots-1,12); 
synapse=zeros(diam+1,diam+1); %total image

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

if nroslots>size(datasqr,1)*n
    disp('Error : not enough available slots with this synaptic size')
    disp(' ')
    return
end
datasqr(:,3)=zeros(size(datasqr,1),1);

sidey=diamarea; %size cluster
sidex=sidey;

if n==1 %one cluster
   % sidey=300;
    %sidey=diamarea;
   % sidex=sidey;
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

    if dist==10
       % sidey=120; %max distance in x or y 
       % sidex=sidey;
        poscentx(1)=-85+centro; 
        poscenty(1)=centro;
        poscentx(2)=85+centro;
        poscenty(2)=centro;

    elseif dist==15
        
      %  sidey=210; %max distance in x or y 
       % sidex=sidey;
        poscentx(1)=-125+centro; 
        poscenty(1)=centro;
        poscentx(2)=125+centro;
        poscenty(2)=centro;
        
    else
        disp('Error : distance between sites can only be 10 or 15 nm')
        disp(' ')
        return
    end

elseif n==4
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %      x x
        %      x x 
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if dist==10
      %  sidey=90; %max distance in x or y
      %  sidex=sidey;
        poscentx(1)=-67+centro; 
        poscenty(1)=-67+centro;
        poscentx(2)=67+centro;
        poscenty(2)=-67+centro;
        poscentx(3)=-67+centro;
        poscenty(3)=67+centro;
        poscentx(4)=67+centro;
        poscenty(4)=67+centro;
        
    elseif dist==15
        
      %  sidey=150; %max distance in x or y
      %  sidex=sidey;
        poscentx(1)=-110+centro; 
        poscenty(1)=-110+centro;
        poscentx(2)=110+centro;
        poscenty(2)=-110+centro;
        poscentx(3)=-110+centro;
        poscenty(3)=110+centro;
        poscentx(4)=110+centro;
        poscenty(4)=110+centro;
        
    else
        disp('Error : distance between sites can only be 10 or 15 nm') 
        disp(' ')
        return
    end

elseif n==7
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %      x x
    %     x x x 
    %      x x
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if dist==10
        sidex=75;
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
        
    elseif dist==15
        
        sidex=100;
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
        disp('Error : distance between sites can only be 10 or 15 nm') 
        disp(' ')
        return
    end

end

totnumbernano=ceil(nroslots/n);
incrementnro=ceil(totnumbernano/100)*100;
numbersites=0;
maxsites=[]; %max number of sites (for area calculation)

        
for i=1:n % all nanoclusters

    countslots=incrementnro*(i-1)+1;
    datasqr(:,3)=0;

    % position of nanocluster's center fixed
    centx=poscentx(i);
    centy=poscenty(i);
    posx=round(centx-sidex/2);
    posy=round(centy-sidey/2);
    
    lastslot=nroslots*(i-1);
    if numbersites+totnumbernano>nroslots % too many...
        totnumbernano=nroslots-numbersites;
    end
    count=1;
    iter=1;

    while count<totnumbernano+1
      index=round(rand*size(datasqr,1));
       if index==0; index=1; end
       if index>size(datasqr,1); index=size(datasqr,1); end
       if datasqr(index,3)==0
            posi=datasqr(index,1)+posx;
            posj=datasqr(index,2)+posy;
            
            if posi <0 || posj <0 || posi>size(synapse,1) || posj>size(synapse,2) 
                disp('Error : verify the number of sites of the hexagonal grid, the size and the number of clusters')
                disp(' ')
                return
            end
            
            synapse(posi-1:posi+1,posj-1:posj+1)=-countslots; %forbidden area
            
            datasqr(index,3)=1; % already chosen
            matrixslots(count+lastslot,1)=countslots;
            matrixslots(count+lastslot,11)=posi;
            matrixslots(count+lastslot,12)=posj;
            count=count+1;
            countslots=countslots+1;
            numbersites=numbersites+1;
       end
       maxsites=[maxsites; datasqr(:,1)+posx datasqr(:,2)+posy];

       iter=iter+1;
       if iter>500000
           disp('Not succeded to find enough positions for sites. Verify data')
           disp(' ')
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

aux=[resto synapse resto];
aux2=zeros(size(resto,2),size(aux,2));
newsynapse=[aux2;aux;aux2];
%   figure
%   imshow(synapse,'InitialMagnification','fit')

matrixslots=matrixslots(find(matrixslots(:,1)~=0),:);
 
%mask
[k,maxarea] = convhull(maxsites(:,1),maxsites(:,2)); %area calculated with ALL posible sites (not only those occupied)
mask=zeros(size(synapse,1),size(synapse,2));
mask = roipoly(mask, maxsites(k,2),maxsites(k,1)) ;
aux=[resto mask resto];
aux2=zeros(size(resto,2),size(aux,2));
finalmask=[aux2;aux;aux2];

indexsyn=[];
for i=1:size(newsynapse,1)
    index=find(newsynapse(i,:)>999);
    if isempty(index)==0
       indexsyn=[indexsyn; i*ones(size(index,2),1) index'];
    end
end

clear  synapse2 numberpergroup

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%