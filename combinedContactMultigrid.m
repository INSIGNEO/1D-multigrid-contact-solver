% solveOneBody(5,3)
% solveTwoBody(3,4,1)
coarsTwoBody(3,4,1,2)


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
        %finding it hard to follow this part of the code
        u=jacSmooth(w,K,b,u,n);
        [inter,rest,Kc]=mlOper(K);
        rc=rest*(b-K*u);
        u=u+inter*mlSolve(w,Kc,rc,0*rc,n,lev-1);
        u=jacSmooth(w,K,b,u,n);
    else %at largest level use PCG to solve 
        m=diag(diag(K));
        u=pcg(K,b,[],[],m);
    end
end


function uNext=jacSmooth(w,K,b,uPrev,n)
    for ii=1:n %for all ii
        uNext=uPrev+w*diag(1./diag(K))*(b-K*uPrev); %jacobien solver formula on wikipedia
        uPrev=uNext;
    end
end

%Two MG Code

function []=coarsTwoBody(M,N,x20,nlev)

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
 
    KmatM = sparse(ii,jj,vv,Mnod+Nnod,Mnod+Nnod,(Mnod+Nnod)*3); %forming stiffness matrix


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
    KmatN = sparse(Mnod+ii,Mnod+jj,vv,Mnod+Nnod,Mnod+Nnod,(Mnod+Nnod)*3); %adjoining matrices
    
    KcM = KmatM;
    KcN = KmatN;
    
    for n = 1:nlev
        
    [interM,restM,KcM]=mlOper(KcM)
    
    end
    
    for n = 1:nlev
        
    [interN,restN,KcN]=mlOper(KcN)
    
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
