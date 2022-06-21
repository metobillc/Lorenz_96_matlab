function nature_plots(Zt, tsteps, tF, spinup)
% Sample plots to evaluate nature run

[Nsteps, Nx] = size(Zt);
fprintf('Size of Zt is %d x %d\n',size(Zt))
figure;
set(gcf,'Position',[100 50 500 700])
% suptitle was apparently deprecated, so
sgtitle('Sample Lorenz ''96 Nature Run Plots')

% Top panel is time series for a particular location
subplot(311)
initial=1:spinup;
spunup=(spinup+1):Nsteps;
node = 1;
plot(tsteps(initial), Zt(initial, node), 'r--')
hold on
plot(tsteps(spunup), Zt(spunup, node), 'b-')
ylabel(sprintf('Time series, Node %d',node))
xlabel('Time [s]')
ylim([-10 20])
grid on

% Middle panel is a state plot at an early time, which may be
% during spinup (if so, plot style will change)
subplot(312)
final_time = size(Zt,1);
time_fraction = 0.05;
plttime = round(time_fraction * final_time);
if plttime > spinup
    cspec = 'b-';
else
    cspec = 'r--';
end
plot(1:Nx, Zt(plttime, :), cspec)
ylabel(sprintf('Ring Map, T= %6.2f',plttime*tF))
xlabel('Location')
ylim([-10 20])
xlim([1 Nx])
grid on

% Bottom panel is a state plot at the final time
subplot(313)
plttime = round(final_time);
if plttime > spinup
    cspec = 'b-';
else
    cspec = 'r--';
end
plot(1:Nx, Zt(plttime, :), cspec)
ylabel(sprintf('Ring Map, T= %6.2f',plttime*tF))
xlabel('Location')
ylim([-10 20])
xlim([1 Nx])
grid on
end