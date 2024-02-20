clear all
close all

%% Basic Execution
%%% Each column displays the values for one dimension.
% Settings IDSFID Signal
N = 300;
dim = 2;
T =  [5, 5];
Ts = 1;
lambda = [0.5, 0.5];
%
%%% Signal Generation
[u, uProxy, yProxy, levels, visitsLevels] = IDSFIDARX(N, dim, T, Ts, lambda);

%% Execution with Advanced Settings
%
% N = 100;
% dim = 2;
% T =  [5, 5];
% Ts = 1;
% lambda = [1, 1];
% M = [30, 30];
% amplitudeConstraints = [0 -5; 5 0];                                      % Minimum (1st row) and maximum (2nd row) value for each input.
% velocityConstraints = [-1 -2; 2 3];                                      % Minimum (1st row) and maximum (2nd row) rate  of change per time step for each input.
% additionalConstraints = 'yes';
% existingData = [ones(5,1)*5, ones(5,1)*-5];
% Lmax = [30, 40];
% n = [4, 4];
% 
% %%% Signal Generation
% [u, uProxy, yProxy, levels, visitsLevels] = IDSFIDARX(N, dim, T, Ts, lambda, M, amplitudeConstraints, velocityConstraints, additionalConstraints);



%% Signal Properties
LineWidth = 2;
sg = 20;

for jj = 1:dim
    % Signal trajectory of the output.
    if jj == 1
        figure
        k = 1:N;
        kTicks = [1, N];
        yTicks = [0, 1];
        plot(k, yProxy, 'LineWidth', LineWidth, 'Color', [0 0.4470 0.7410], "LineStyle","-");
        set(gca, 'XTickMode', 'manual', 'XTick', kTicks, 'xlim', [1, N],'fontsize',sg);
        set(gca, 'YTickMode', 'manual', 'YTick', yTicks, 'ylim', [-0.05, 1.05],'fontsize',sg);
        set(gca,'TickLabelInterpreter','latex')
        xlabel('$k$','interpreter','latex','fontsize',sg);
        ylabel('$\tilde{y}(k)$','interpreter','latex','fontsize',sg);
        pos1=1120; pos2=0; width=560; height=210;
        set(gcf,'position',[pos1,pos2,width,height])
    end

    % Signal trajectory of the inputs.
    urange = max(u(:,jj)) - min(u(:,jj));
    figure
    k = 1:N;
    kTicks = [1, N];
    uTicks = [min(u(:,jj)), max(u(:,jj))];
    plot(k, u(:,jj), 'LineWidth', LineWidth, 'Color', [0 0.4470 0.7410], "LineStyle","-");
    set(gca, 'XTickMode', 'manual', 'XTick', kTicks, 'xlim', [1, N],'fontsize',sg);
    set(gca, 'YTickMode', 'manual', 'YTick', uTicks, 'ylim', [min(u(:,jj))-0.05*urange, max(u(:,jj))+0.05*urange],'fontsize',sg);
    set(gca,'TickLabelInterpreter','latex')
    xlabel('$k$','interpreter','latex','fontsize',sg);
    ylabel('$u(k)$','interpreter','latex','fontsize',sg);
    pos1=0; pos2=0; width=560; height=210;
    set(gcf,'position',[pos1,pos2,width,height])
    titles = ['Input', num2str(jj)];
    title(titles);

    %%% Proxy regressor space of the inputs.
    figure
    yTicks = [0, 1]; %ceil(max(yProxy*10))/10];
    uTicks = [0, 1]; %max(uProxy(:,jj))];
    xLabel = ['$\tilde{u}_{', num2str(jj), '}(k-1)$'];
    scatter(uProxy(:,jj), yProxy, 'LineWidth', LineWidth, 'Color', [0 0.4470 0.7410])
    set(gca, 'XTickMode', 'manual', 'XTick', uTicks, 'xlim', [-0.025, 1+0.025],'fontsize',sg);
    set(gca, 'YTickMode', 'manual', 'YTick', yTicks, 'ylim', [-0.025, 1+0.025],'fontsize',sg);
    set(gca,'TickLabelInterpreter','latex')
    xlabel(xLabel,'interpreter','latex','fontsize',sg);
    ylabel('$\tilde{y}(k-1)$','interpreter','latex','fontsize',sg);
    pos1=0; pos2=500; width=550; height=550;
    set(gcf,'position',[pos1,pos2,width,height])
    titles = ['Input', num2str(jj)];
    title(titles);

    %%% Spectrum
    fs = 1/Ts;
    df = fs/(N);
    fk = 0 : df : fs/2;
    udetr = detrend(u(:,jj));
    U = fft(udetr);
    UamplSp = 2*df*abs(U)/fs;

    fTicks = linspace(fk(1), fk(end), 6);
    Umin = 0; Umax = max(UamplSp)+0.025*max(UamplSp);
    UTicks = [Umin, Umax];

    figure
    plot(fk, UamplSp(1:length(fk)),'LineWidth',LineWidth-0.5,'Color',[0 0.4470 0.7410]);
    set(gca, 'XTickMode', 'manual', 'XTick', fTicks, 'xlim', [fk(1), fk(end)],'fontsize',sg);
    set(gca, 'YTickMode', 'manual', 'YTick', UTicks, 'ylim', [0, Umax],'fontsize',sg);
    set(gca,'TickLabelInterpreter','latex')
    xlabel('$f \ [\mathrm{Hz}]$','interpreter','latex','fontsize',sg);
    ylabel('$|U(f)| $','interpreter','latex','fontsize',sg);
    pos1=562; pos2=0; width=560; height=210;
    set(gcf,'position',[pos1,pos2,width,height])
    titles = ['Input', num2str(jj)];
    title(titles);

    %%% Visited Levels

    % Recommendations for determining the number of levels, M:
    % Ideally, each level should be visited exactly once, but achieving
    % this can be challenging. Particularly when dealing with multiple
    % inputs, it is recommended to generate a Latin Hypercube that is
    % (slightly) overpopulated (ensuring that each level is visited more
    % than once) rather than underpopulated.

    XTicks = [0, median(levels{jj}), max(levels{jj})];
    YTicks = 0 : 1 : max(visitsLevels{jj})+1;
    xLabel = ['$\tilde{u}_{', num2str(jj), '}-$value'];

    figure
    bar(levels{jj}, visitsLevels{jj});
    set(gca, 'XTickMode', 'manual', 'XTick', XTicks, 'xlim', [0-1/(2*length(visitsLevels{jj})), max(levels{jj}+1/(2*length(visitsLevels{jj})))],'fontsize',sg);
    set(gca, 'YTickMode', 'manual', 'YTick', YTicks, 'ylim', [0, max(visitsLevels{jj})+1],'fontsize',sg);
    set(gca,'TickLabelInterpreter','latex')
    xlabel(xLabel,'interpreter','latex','fontsize',sg);
    ylabel('Number of Visits','interpreter','latex','fontsize',sg);
    pos1=562; pos2=500; width=550; height=420;
    set(gcf,'position',[pos1,pos2,width,height])
    titles = ['Input', num2str(jj)];
    title(titles);
end
