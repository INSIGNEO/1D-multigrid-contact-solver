function [x1,u1,x2,u2]=solveTwoBody(M,N,x20)
Nel1=2^M; Nno1=Nel1+1; Nel2=2^N; Nno2=Nel2+1;

x1=(0:(2^M))/(2^M);
ki=0.5;kf=5;
k1=(2^M)*(ki+(kf-ki)*0.5*(x1(1:end-1)+x1(2:end)));

x2=(2^(N-M))*(0:(2^N))/(2^N);
kfix=2.5;
k2=(2^N)*(kfix+0.0*0.5*(x2(1:end-1)+x2(2:end)));
x2=x20+x2;
g0=x2(1)-x1(end);

ii=[1;reshape(repmat(2:(Nno1-1),3,1),(Nno1-2)*3,1);Nno1;Nno1];
jj=(2:(Nno1-1)); jj=[jj-1;jj;jj+1];
jj=[1;reshape(jj,(Nno1-2)*3,1);Nno1-1;Nno1];
vv=[k1(1:end-1);k1(2:end)];
vv=[-vv(1,:); vv(1,:)+vv(2,:); -vv(2,:)];
vv=[reshape(vv,(Nno1-2)*3,1);-k1(end);k1(end)];
if(length(k1)==1)
    vv=[k1(1);vv];
else
    vv=[k1(1)+k1(2);vv];
end
Kmat=sparse(ii,jj,vv,Nno1+Nno2,Nno1+Nno2,(Nno1+Nno2)*3);

ii=[1;1;reshape(repmat(2:(Nno2-1),3,1),(Nno2-2)*3,1);Nno2];
jj=(2:(Nno2-1)); jj=[jj-1;jj;jj+1];
jj=[1;2;reshape(jj,(Nno2-2)*3,1);Nno2];
vv=[k2(1:end-1);k2(2:end)];
vv=[-vv(1,:); vv(1,:)+vv(2,:); -vv(2,:)];
vv=[k2(1);-k2(1);reshape(vv,(Nno2-2)*3,1)];
if(length(k2)==1)
    vv=[vv;k2(end)];
else
    vv=[vv;k2(end-1)+k2(end)];
end
Kmat=Kmat+sparse(Nno1+ii,Nno1+jj,vv,Nno1+Nno2,Nno1+Nno2,(Nno1+Nno2)*3);

U0=-g0-0.1;
inc=10;
kc=1e3*0.5*(k1(end)+k2(1));
tol=1e-6;
Kmatc=Kmat+sparse(Nno1+[0;0;1;1],Nno1+[0;1;0;1], ...
    kc*[-1;1;1;-1],Nno1+Nno2,Nno1+Nno2,4);
rc=kc*sparse(Nno1+(0:1)',[1;1],[-1;1],Nno1+Nno2,1,2);
for ii=1:inc
    if(length(k2)==1)
        vv=(ii*U0/inc)*k2(end);
    else
        vv=(ii*U0/inc)*(k2(end-1)+k2(end));
    end
    b=sparse(Nno1+Nno2,1,vv,Nno1+Nno2,1,1);
    u=Kmat\b;
    overc=u(Nno1+1)-u(Nno1)+g0;
    du=inf;
    if (overc<0)
        while (max(abs(du))>tol*abs(U0/inc))
            du=Kmatc\(b+rc*overc-Kmat*u);
            u=u+du;
            overc=u(Nno1+1)-u(Nno1)+g0;
        end
    end
end
u1=u(1:Nno1);u2=u(Nno1+(1:Nno2));
end
