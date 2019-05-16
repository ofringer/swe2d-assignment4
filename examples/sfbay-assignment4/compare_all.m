fnames = {'theta_25_wind_1_tides_1',...
          'theta_25_wind_0_tides_1',...
          'theta_25_wind_1_tides_0',...
          'theta_90_wind_1_tides_1',...
          'theta_90_wind_0_tides_1',...
          'theta_90_wind_1_tides_0'};
legend_text = {'-25, wind+tides',...
               '-25, no wind',...
               '-25, no tides',...
               '90, wind+tides',...
               '90, no wind',...
               '90, no tides'};

cs = {'k-','k--','k:','r-','r--','r:'};

nout=3;
n_avg = fix(24*3600/dt);
num_days = fix(nsteps/n_avg);
t0 = [1:nsteps/nout]*dt*nout;

Flux = zeros(length(fnames),num_days);
C_all = zeros(length(fnames),length(t0));
zb_all = zeros(length(fnames),length(t0));

for fn=1:length(fnames)
    load(sprintf('solution_data/%s',fnames{fn}));

    C_all(fn,:) = C_out;
    zb_all(fn,:) = zb_out;
    for m=1:num_days
        f = sum(Q_vals{1}.*C_vals{1},2);
        Flux(fn,m) = sum(f(1+(m-1)*n_avg:m*n_avg))*dt/1e6;
    end
end

figure(1)
clf;
hold on;
set(gca,'box','on');
for fn=1:length(fnames)
    plot([1:num_days],Flux(fn,:),cs{fn});
end
legend(legend_text);
xlabel('t (days)');
ylabel('Daily flux (kton)');
    
figure(2)
clf;

subplot(2,1,1)
hold on;
set(gca,'box','on');
for fn=1:length(fnames)
    plot(t0,C_all(fn,:)*1000,cs{fn});
end
ylabel('C (mg/L');

subplot(2,1,2)
hold on;
set(gca,'box','on');
for fn=1:length(fnames)
    plot(t0,zb_all(fn,:)*1000,cs{fn});
end
ylabel('z_b (mm)');
xlabel('t (days)');

for fn=1:length(fnames)
    fprintf('%s: %.2e kton/day\n',legend_text{fn},Flux(fn,end));
end



