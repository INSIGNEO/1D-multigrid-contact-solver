% solveOneBody(5,3)
% solveTwoBody(3,4,1)
% coarsTwoBody(3,4,1,2)

M=5;
N=5;
nlev=3;
x20=1;

mlTwoBodySolver(M,N,nlev,x20)

% [x1,u1,x2,u2]=solveTwoBody(M,N,x20);
% u1
% u2

%Two Body Functions
function [x1,u1,x2,u2]=solveTwoBody(M,N,x20)
Nel1=2^M; Nno1=Nel1+1; Nel2=2^N; Nno2=Nel2+1; %Defines number of elemements and nodes could be better noted with M

x1=(0:(2^M))/(2^M);%lengh vector 1
ki=0.5;kf=5; %defining stiffness bounds
k1=(2^M)*(ki+(kf-ki)*0.5*(x1(1:end-1)+x1(2:end))); %defining stiffness vector

x2=(2^(N-M))*(0:(2^N))/(2^N);%lengh vector 2
kfix=2.5;% stiffness second spring
k2=(2^N)*(kfix+0.0*0.5*(x2(1:end-1)+x2(2:end)));%stiffness vector feel like it could be simplar
x2=x20+x2;%shifted spring
g0=x2(1)-x1(end); %gap

ii=[1;reshape(repmat(2:(Nno1-1),3,1),(Nno1-2)*3,1);Nno1;Nno1];%element colum
jj=(2:(Nno1-1)); jj=[jj-1;jj;jj+1]; %element ordering row
jj=[1;reshape(jj,(Nno1-2)*3,1);Nno1-1;Nno1]; %change to vector
vv=[k1(1:end-1);k1(2:end)]; 
vv=[-vv(1,:); vv(1,:)+vv(2,:); -vv(2,:)];
vv=[reshape(vv,(Nno1-2)*3,1);-k1(end);k1(end)];%values for stiffness matrix
if(length(k1)==1)
    vv=[k1(1);vv];
else
    vv=[k1(1)+k1(2);vv];%applying BC?
end
Kmat=sparse(ii,jj,vv,Nno1+Nno2,Nno1+Nno2,(Nno1+Nno2)*3);%forming stiffness matrix

ii=[1;1;reshape(repmat(2:(Nno2-1),3,1),(Nno2-2)*3,1);Nno2]; %element colum K2
jj=(2:(Nno2-1)); jj=[jj-1;jj;jj+1];%element ordering row
jj=[1;2;reshape(jj,(Nno2-2)*3,1);Nno2];%change to vector
vv=[k2(1:end-1);k2(2:end)];
vv=[-vv(1,:); vv(1,:)+vv(2,:); -vv(2,:)];
vv=[k2(1);-k2(1);reshape(vv,(Nno2-2)*3,1)];%values for stiffness matrix
if(length(k2)==1)
    vv=[vv;k2(end)];
else
    vv=[vv;k2(end-1)+k2(end)];%applying BC?
end
Kmat=Kmat+sparse(Nno1+ii,Nno1+jj,vv,Nno1+Nno2,Nno1+Nno2,(Nno1+Nno2)*3); %adjoining matrices

U0=-g0-0.1;%-gap-spring length (penalty lenth?)
inc=10; %incremtns
kc=1e3*0.5*(k1(end)+k2(1)); %penalty spring
tol=1e-6; %tolerance
Kmatc=Kmat+sparse(Nno1+[0;0;1;1],Nno1+[0;1;0;1], ...
    kc*[-1;1;1;-1],Nno1+Nno2,Nno1+Nno2,4); %stiffness matrix inculuding penalty springs
rc=kc*sparse(Nno1+(0:1)',[1;1],[-1;1],Nno1+Nno2,1,2); %isolated just penalty springs (oposite sign) pulling canceler
for ii=1:inc %iterate till solution final
    if(length(k2)==1) % penalty value
        vv=(ii*U0/inc)*k2(end);
    else
        vv=(ii*U0/inc)*(k2(end-1)+k2(end));
    end
    
    b=sparse(Nno1+Nno2,1,vv,Nno1+Nno2,1,1); %BC and location
    u=Kmat\b; %displacment finded
    overc=u(Nno1+1)-u(Nno1)+g0; %over corrention is has been pushed in too far? (over closure)
    du=inf;
    
    if (overc<0)%if penertarting
        while (max(abs(du))>tol*abs(U0/inc))%moving back spring
            du=Kmatc\(b+rc*overc-Kmat*u); %small push back factor
            u=u+du;%push back
            overc=u(Nno1+1)-u(Nno1)+g0; %redfine over correction
        end
    end
end
u1=u(1:Nno1);u2=u(Nno1+(1:Nno2)); %final answer
end

%One Body Functions

function [x,u]=solveOneBody(N,nlev)
    Nel=2^N; %Defines Number of Elements
    Nno=Nel+1; %Defines Numnber of Nodes
    
    x=(0:(2^N))/(2^N); %Displacments Vector (0-1)
    k=(2^N)*(0.5+4.5*0.5*(x(1:end-1)+x(2:end))); %Stiffness vector for nodes defined by function
    
    
    ii=[1;reshape(repmat(2:(Nno-1),3,1),1,(Nno-2)*3,1)';Nno]; %ii Repmat - Replicate and tile an array ReShape-Reshapes repmat 1XNel+1X1 ' was missing so dimentions not work added it seesm to work now
    jj=(2:(Nno-1));jj=[jj-1;jj;jj+1]; %node numbering file
    jj=[1;reshape(jj,(Nno-2)*3,1);Nno]; %forming triplet vector
    vv=-[k(1:end-1);k(2:end)]; %concainated two spring rows (K) together of 1,2(combine with next line can likely be done)
    vv=[vv(1,:); -vv(1,:)-vv(2,:); vv(2,:)]; %vv combined with (k1+k2) in central row
    vv=[k(1)+k(2);reshape(vv,(Nno-2)*3,1);k(end-1)+k(end)]; %k(1)+k(2) BC and vv combines into vector ending with k(n-1)+k(n) BC
    Kmat=sparse(ii,jj,vv,Nno,Nno,2+(Nno-2)*3);% forms tridigonal stiffness matrix 
  
    b=[zeros(Nno-1,1);-0.01*(k(end-1)+k(end))];%no.element long zero matrix and concat -(k(n-1)+k(n)) on end BC
    
    u=[zeros(Nno-1,1);-0.01]; %no.element long zero matrix and concat -0.01 first guess
    
    % u=-0.01*(0:(Nno-1))'/(Nno-1);
    
    u=mlSolve(0.61,Kmat,b,u,2,nlev); %Solve functions
    
    % m=diag((diag(Kmat)).^0.5);
    % u=pcg(Kmat,b,[],[],m,m,u);
    % u=Kmat\b;
end

function u=mlSolve(w,K,b,u,n,lev)

    if lev>1 %if level is greater that 1
        
        u=jacSmooth(w,K,b,u,n); %jacobi
        [inter,rest,Kc]=mlOper(K); %coarsen
        rc=rest*(b-K*u); %residuals
        u=u+inter*mlSolve(w,Kc,rc,0*rc,n,lev-1);
        u=jacSmooth(w,K,b,u,n);
    else %at largest level use PCG to solve 
        m=diag(diag(K));
        u=pcg(K,b,[],[],m);
    end
end


%Two MG Code

function [KmatM,KmatN,uM,uN,bM,bN,kc,g0,kN]=matrixCreation(M,N,x20)%Creates M and N matrix and u

    Mel=2^M; Mnod=Mel+1; Nel=2^N; Nnod=Nel+1; %Defines number of elemements and nodes could be better noted with M

    xM=(0:(Mel))/(Mel);%lengh vector 1
    ki=0.5;kf=5; %defining stiffness bounds
    kM=(Mel)*(ki+(kf-ki)*0.5*(xM(1:end-1)+xM(2:end))); %defining stiffness vector

    xN=(2^(N-M))*(0:(2^N))/(2^N);%lengh vector 2
    kfix=2.5;% stiffness second spring
    kN=(Nel)*(kfix+0.0*0.5*(xN(1:end-1)+xN(2:end)));%stiffness vector feel like it could be simplar
    
    xN=x20+xN;%shifted spring
    g0=xN(1)-xM(end); %gap
    
    ii=[1;reshape(repmat(2:(Mnod-1),3,1),(Mnod-2)*3,1);Mnod;Mnod];%element colum

    jj=(2:(Mnod-1)); jj=[jj-1;jj;jj+1]; %element ordering row
    jj=[1;reshape(jj,(Mnod-2)*3,1);Mnod-1;Mnod]; %change to vector
    vv=[kM(1:end-1);kM(2:end)]; 
    vv=[-vv(1,:); vv(1,:)+vv(2,:); -vv(2,:)];
    vv=[reshape(vv,(Mnod-2)*3,1);-kN(end);kN(end)];%values for stiffness matrix
    
    if(length(kM)==1)
        vv=[kM(1);vv];
    else
        vv=[kM(1)+kM(2);vv];%applying BC?
    end
 
    KmatM = sparse(ii,jj,vv,Mnod,Mnod,(Mnod)*3); %forming stiffness matrix


    ii=[1;1;reshape(repmat(2:(Nnod-1),3,1),(Nnod-2)*3,1);Nnod]; %element colum K2
    jj=(2:(Nnod-1)); jj=[jj-1;jj;jj+1];%element ordering row
    jj=[1;2;reshape(jj,(Nnod-2)*3,1);Nnod];%change to vector
    vv=[kN(1:end-1);kN(2:end)];
    vv=[-vv(1,:); vv(1,:)+vv(2,:); -vv(2,:)];
    vv=[kN(1);-kN(1);reshape(vv,(Nnod-2)*3,1)];%values for stiffness matrix
    
    if(length(kN)==1)
        vv=[vv;kN(end)];
    else
        vv=[vv;kN(end-1)+kN(end)];%applying BC?
    end

    KmatN = sparse(ii,jj,vv,Nnod,Nnod,(Nnod)*3);
%     KmatN = sparse(Mnod+ii,Mnod+jj,vv,Mnod+Nnod,Mnod+Nnod,(Mnod+Nnod)*3); %forming N matrix
    
  
    %Defining U guess
    uM=-0.01*(0:(Mnod-1))'/(Mnod-1);
    uN=-0.01*(0:(Nnod-1))'/(Nnod-1);

    %Defining B
    bM=[zeros(Mnod-1,1);-0.01*(kM(end-1)+kM(end))];%no.element long zero matrix and concat -(k(n-1)+k(n)) on end BC
    bN=[zeros(Nnod-1,1);-0.01*(kN(end-1)+kN(end))];

    kc=1e3*0.5*(kM(end)+kN(1)); %penalty spring


end

function [KcM,KcN,ucM,ucN,bcM,bcN,interM,interN]=coarsTwoBody(nlev,KcM,KcN,uM,uN,bM,bN)

    ucM = uM;
    bcM = bM;
    ucN = uN;
    bcN = bN;

    
    for n = 1:nlev
        
    [interM,restM,KcM]=mlOper(KcM);
    ucM = restM*ucM;
    bcM = restM*bcM;

    end
    
    for n = 1:nlev
        
    [interN,restN,KcN]=mlOper(KcN);
    ucN = restN*ucN;
    bcN = restN*bcN;


    end
    
%     [interN,restN,KcN]=mlOper(KmatN)
    
    
% KmatN =Kmat+sparse(Nno1+ii,Nno1+jj,vv,Nno1+Nno2,Nno1+Nno2,(Nno1+Nno2)*3); %adjoining matrices

    
end

function [inter,rest,Kc]=mlOper(K)
    
    sz1=size(K,1);
    szc1=floor((sz1+1)/2);%removed plus one from here seems to work but could creat bug down line? %added floor instead for rounding to stop errors

    ii=[1;reshape(repmat((1:2:sz1-4),3,1)+(1:3)',[],1);sz1];
    jj=[1;reshape(repmat(2:(szc1-1),3,1),[],1);szc1];
    ss=[1;repmat([0.5;1;0.5],szc1-2,1);1]; 
    
    
    inter=sparse(ii,jj,ss,sz1,szc1,3*szc1+2);
    [ii,jj,ss]=find(inter(2:end-1,2:end-1));
    %inter grid to refine mesh
    rest=sparse([1;jj+1;szc1],[1;ii+1;sz1],...
        [1;0.5*ss;1],szc1,sz1,3*szc1+2);
    %restriction grid 
    Kc=rest*K*inter;

    %coauresend C
end

function [uSolve]=FirstContactSolver(KcMat,uSolve,bSolve,kc,g0,McNod,NcNod,kN)



    U0=-g0-0.1;%-gap-spring length (penalty lenth?)
    inc=10; %incremtent
    tol=1e-6; %tolerance

    kMatCon=KcMat+sparse(McNod+[0;0;1;1],McNod+[0;1;0;1], ...
        kc*[-1;1;1;-1],McNod+NcNod,McNod+NcNod,4); %stiffness matrix inculuding penalty springs

    rc=kc*sparse(McNod+(0:1)',[1;1],[-1;1],McNod+NcNod,1,2); %isolated just penalty springs (oposite sign) pulling canceler

    for ii=1:inc %iterate till solution final
        if(length(kN)==1) % penalty value
            vv=(ii*U0/inc)*kN(end);
        else
            vv=(ii*U0/inc)*(kN(end-1)+kN(end));
        end
        
        b=sparse(McNod+NcNod,1,vv,McNod+NcNod,1,1); %BC and locatiom
%         uSolve=Kmat\b; %displacment finded
        uSolve = pcg(KcMat,b);
        
        overc=uSolve(McNod+1)-uSolve(McNod)+g0; %over corrention is has been pushed in too far? (over closure)
        du=inf;
        
        if (overc<0)%if penertarting
            while (max(abs(du))>tol*abs(U0/inc))%moving back spring
                du=kMatCon\(b+rc*overc-KcMat*uSolve); %small push back factor
%                 du = pcg(kMatCon,(b+rc*overc-KcMat*uSolve)); %chech how to use pcg here
                uSolve=uSolve+du;%push back
                overc=uSolve(McNod+1)-uSolve(McNod)+g0; %redfine over correction
            end 
        end
    end
%     u1=u(1:Nno1);u2=u(Nno1+(1:NcNod)); %final answer
end

function uNext=jacSmooth(w,K,b,uPrev,n)
    for ii=1:n %for all ii
        uNext=uPrev+w*diag(1./diag(K))*(b-K*uPrev); %jacobien solver formula on wikipedia
        uPrev=uNext;
    end
end


function [x1,u1,x2,u2] = mlTwoBodySolver(M,N,nlev,x20)

    [KmatM,KmatN,uM,uN,bM,bN,kc,g0,kN] = matrixCreation(M,N,x20); %Matrix Creation
    KcM = KmatM;
    KcN = KmatN; %Redfining matrix to be coursened matrix
   
    [KcM,KcN,ucM,ucN,bcM,bcN,interM,interN]=coarsTwoBody(nlev,KcM,KcN,uM,uN,bM,bN); %Coarsening bodies together

    McNod = size(ucM,1); 
    NcNod = size(ucN,1); %No. of coarsened nodes
    
    [iiM, jjM, vvM] = find(KcM);
    [iiN, jjN, vvN] = find(KcN);
    
    KcMat = sparse([iiM; iiN + size(KcM,1)], [jjM; jjN+size(KcM,1)], [vvM; vvN]); %Concat K matrix together (check logic)
    uSolve= vertcat(ucM,ucN);%combinging first guess
    bSolve= vertcat(bcM,bcN);%this line may not be needed

    [uSolve]=FirstContactSolver(KcMat,uSolve,bSolve,kc,g0,McNod,NcNod,kN) %Coarsest level of contact

    %Refined levels of contact
    for currentLev = 1:(nlev-1)

        currentLev = nlev-currentLev;%So it counts down

        uGM=uSolve(1:McNod); %Spliting U for M and N
        uGN=uSolve(McNod+(1:NcNod));

        uGM=interM*uGM; %interpolatioing guess using previous matrix
        uGN=interN*uGN;

        [KcM,KcN,ucM,ucN,bcM,bcN,interM,interN]=coarsTwoBody(currentLev,KmatM,KmatN,uM,uN,bM,bN);
        
        ucM = uGM;%override gentered u with guessed u
        ucN = uGN;
        
        McNod = size(ucM,1); 
        NcNod = size(ucN,1); %No. of coarsened nodes
    
        [iiM, jjM, vvM] = find(KcM);%break down coarsended matrix to be reco
        [iiN, jjN, vvN] = find(KcN);
        
        KcMat = sparse([iiM; iiN + size(KcM,1)], [jjM; jjN+size(KcM,1)], [vvM; vvN]); %Concat K matrix together (check logic)
        uSolve= vertcat(ucM,ucN);%combinging first guess
        bSolve= vertcat(bcM,bcN);%this line may not be needed
        
        

    end



    
    u1=uSolve(1:McNod);
    u2=uSolve(McNod+(1:NcNod)); %final answer

end
