% solveOneBody(5,3)
% solveTwoBody(3,4,1)
M=6;
N=6;
nlev=1;
x20=1;
nMax=120;  

vecJacNo=zeros(nMax,1);
vecError=zeros(nMax,1);

% legendstrings=cell(10,1);

hold on 

[x,u]=solveTwoBody(M,N,x20);



% for nlev =1:4
% 
%     for n = 1:1:nMax
%     
%         [xml,uml]=solveTwoBodyML(M,N,x20,nlev,n);
% 
%         error= max(abs(u-uml));
% 
%         jacNo=n;
%         vecJacNo(n) = jacNo;
%         vecError(n)= error;
%     end
% 
%     errorChart(nlev)=plot(vecJacNo,vecError);
%     legendstrings{nlev}=sprintf('%s %d','Multi Level',nlev);
%     vecJacNo=[];
%     vecError=[];
% 
% end
% 
% legend(errorChart,legendstrings);
% title('Error and Jacobien Smoother Iterations')
% xlabel('Jacobien Smoother Iterations')
% ylabel('Max error')


% [x,u]=solveTwoBody(10,10,x20);
% anagraph=plot(x,u);
% [x,u]=solveTwoBody(3,3,x20);
% anagraph=plot(x,u);
% [x,u]=solveTwoBody(6,6,x20);
% anagraph=plot(x,u);
% 
% legend('Analitical Soultion','16 Node Solution','128 Node Solution');
% title('Displacment of a spring')
% 
% ylabel('U Displacement (units)')


% % legend(legendstrings);
% % title('Displacement Caused by Increments')
% % xlabel('Starting Poition of Spring')
% % ylabel('U Displacement (units)')



%Two Body Functions
function [x,u]=solveTwoBody(M,N,x20)

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
x = vertcat (x1',x2');
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
%     incgraph(ii)=plot(x,u);
%     legendstrings{ii}=sprintf('%s %d','Increment ',ii);
end

% plot (x,u)

end


function uNext=jacSmooth(w,K,b,uPrev,n)
    for ii=1:n %for all ii
        uNext=uPrev+w*diag(1./diag(K))*(b-K*uPrev); %jacobien solver formula on wikipedia
        uPrev=uNext;
    end
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


function u=mlSolve(w,K,b,u,n,lev)
    u=jacSmooth(w,K,b,u,5); %jacobi twice

    if lev>1 %if level is greater that 1

        [inter,rest,Kc]=mlOper(K);

        rc=rest*(b-K*u);  %compute residual and restrict it to coasrer grid
        
        u=u+inter*mlSolve(w,Kc,rc,0*rc,n,lev-1); %use mlsolve to compute update %interpolate update and add vecotor to solution

    else %if coarsest 
        m=diag(diag(K));
        u=pcg(K,b,[],[],m);%use pcg
        u=jacSmooth(w,K,b,u,n);
%         u=K\b;
    end

    u=jacSmooth(w,K,b,u,5); %jacobi twice
end


function [x,u]=solveTwoBodyML(M,N,x20,nlev,n)
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

u=-0.01*(0:(Nno1+Nno2-1))'/(Nno1+Nno2-1);
% x = 1:(Nno1+Nno2);
%     x = x.*(2/(Nno1+Nno2+1));
%     plot (x,u)

for ii=1:inc %iterate till solution final
    if(length(k2)==1) % penalty value
        vv=(ii*U0/inc)*k2(end);
    else
        vv=(ii*U0/inc)*(k2(end-1)+k2(end));
    end
    
    b=sparse(Nno1+Nno2,1,vv,Nno1+Nno2,1,1); %BC and location
%     u=Kmat\b; %displacment finded

    u = mlSolve(0.61,Kmat,b,u,n,nlev);

    overc=u(Nno1+1)-u(Nno1)+g0; %over corrention is has been pushed in too far? (over closure)
    du=u*tol*10;
    count = 0;

    if (overc<0)%if penertarting
        while (max(abs(du))>tol*abs(U0/inc)) && count < 10 %moving back spring
%             du=Kmatc\(b+rc*overc-Kmat*u); %small push back factor
            du = mlSolve(0.61,Kmatc,(b+rc*overc-Kmat*u),du,n,nlev);
            u=u+du;%push back
            overc=u(Nno1+1)-u(Nno1)+g0; %redfine over correction
            count = count + 1;
        end
    end
%     x = 1:(Nno1+Nno2);
%     x = x.*(2/(Nno1+Nno2+1));
%     plot (x,u)
    
end

x = vertcat (x1',x2');
% plot (x,u)


end
