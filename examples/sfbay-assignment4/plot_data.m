global C_vals Q_vals A_vals U_vals

n_avg = fix(24*3600/dt);
num_days = fix(nsteps/n_avg);

Flux = zeros(1,num_days);
for m=1:num_days
    f = sum(Q_vals{1}.*C_vals{1},2);
    Flux(m) = sum(f(1+(m-1)*n_avg:m*n_avg))*dt/1e3;
end

fprintf('Flux on 8th day is %.2f ton/day\n',Flux(end));

figure(1)
plot([1:num_days],Flux);
xlabel('Day');
ylabel('Flux (ton/day)');

t0 = [1:length(C_out)]*dt*nout/86400;

figure(2)
subplot(3,1,1)
plot(t0,C_out*1000);
ylabel('C (mg/L)');
subplot(3,1,2)
plot(t0,(h_out-hoffset));
ylabel('h (m)');
subplot(3,1,3)
plot(t0,zb_out*1000);
ylabel('zb (mm)');
xlabel('t (day)');






