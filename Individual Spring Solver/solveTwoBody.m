
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
