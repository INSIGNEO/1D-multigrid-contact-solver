

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

function [inter,rest,Kc]=mlOper(K)
    sz1=size(K,1);
    szc1=(sz1+1)/2;
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

function uNext=jacSmooth(w,K,b,uPrev,n)
    for ii=1:n %for all ii
        uNext=uPrev+w*diag(1./diag(K))*(b-K*uPrev); %jacobien solver formula on wikipedia
        uPrev=uNext;
    end
end

