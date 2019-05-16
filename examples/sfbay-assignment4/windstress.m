function [U10,taux,tauy]=windstress(n,dt)

getvariables

t=n*dt;
if(omega_0==0)
    U10 = Uw0;
else    
    U10 = max(0,Uw0*cos(omega_0*(n+.5)*dt));
end
tau_w = rho_air/rho0*Cdw*U10^2;
taux = tau_w*cos(thetaW);
tauy = tau_w*sin(thetaW);
