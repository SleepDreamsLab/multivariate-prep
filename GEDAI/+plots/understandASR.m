function [] = understandASR(EEG, colors, opts)
arguments (Input)
    EEG
    colors
    opts.Stage = 'N2'
    opts.BoutDuration = [] % in min
    opts.ScoringDuration = 30 % in s
    opts.K = [45]
    opts.iwin = [1]
    opts.A = [];
    opts.B = []; 
end

    % Extract EEG of sleep stage
    EEG_SleepStage = pop_epoch(EEG, opts.Stage, [0 opts.ScoringDuration]);
    EEG_SleepStage = eeg_epoch2continuous(EEG_SleepStage); 

    % Bout Length
    if isempty(opts.BoutDuration)
        BoutLength = EEG_SleepStage.pnts - 1;
    else
        BoutLength = opts.BoutDuration * 60 * EEG.srate;    
    end
    BoutStarts = 1 : BoutLength : (EEG_SleepStage.pnts);
    
    % --- The Loop ---
    %%
    for iwin = opts.iwin
        
        % Extract the current window
        idxwindow = BoutStarts(iwin) : BoutStarts(iwin) + BoutLength - 1;
        
        % Extract respective data
        X = EEG_SleepStage;
        X.data = EEG_SleepStage.data(:, idxwindow);
        X = eeg_checkset(X);

        % Find calibration data
        [XC, mask, Meta] = plots.asr.clean_windows(X); %[ XC, mask] = clean_windows(X);


                %% ----- Figure 1 -----
                figure; set(gcf, 'Color', 'w')
        
                % Layout
                nrow=4; ncol=1;
                layout = tiledlayout(nrow, ncol, ...
                    'TileSpacing', 'compact', ...
                    'Padding', 'compact');
        
                %%% Data of all epochs of respective sleep stage
                nexttile
                Times_min = linspace(0, EEG_SleepStage.pnts/EEG_SleepStage.srate/60, EEG_SleepStage.pnts);         
                plot(Times_min, EEG_SleepStage.data')
                hold on; xline([opts.BoutDuration*(iwin-1), opts.BoutDuration*iwin], 'k-', 'LineWidth', 2)
                xtickformat('%gmin')
                ylabel('\muV')
                xlim([Times_min(1) Times_min(end)])
                title(sprintf('All %s epochs of all %d channels', opts.Stage, EEG_SleepStage.nbchan))
        
                %%% Data of one bout
                nexttile
                Times_min = linspace(0, X.pnts/X.srate/60, X.pnts);    
                plot(Times_min, X.data')
                xtickformat('%gmin')
                ylabel('\muV')
                xlim([Times_min(1) Times_min(end)])
                title(sprintf('Selected %dmin bout', opts.BoutDuration))
                
                %%% RMS values
                nexttile
                Times_min = linspace(0, X.pnts/X.srate/60, Meta.nWindows);
                plot(Times_min, Meta.zRMS', '.', 'MarkerSize', 5, 'HandleVisibility','off')
                hold on; L1=yline([Meta.zThresholds], 'k-', 'DisplayName', sprintf('Calibration thresholds = [%.1f, %.1f]', Meta.zThresholds(1), Meta.zThresholds(2)))
                xtickformat('%gmin')
                ylabel('RMS (z-standardized)')
                xlim([Times_min(1) Times_min(end)])
                title(sprintf('RMS vaues (%ds windows, %.2fs overlap)', Meta.Windowlen, Meta.Overlap))
                legend(L1(1))
        
                %%% Calibration data
                nexttile
                Times_min = linspace(0, X.pnts/X.srate/60, X.pnts);    
                Times_min = Times_min(mask);
                plot(Times_min, X.data(:, mask)', '.', 'MarkerSize', 2)
                xtickformat('%gmin')
                ylabel('\muV')
                xlim([Times_min(1) Times_min(end)])
                title(sprintf('Calibration data (keeping %.1f%%)', sum(mask)/numel(mask)*100))
        
                
        % Compute covariance matrix
        DefaultState = plots.asr.asr_calibrate(XC.data,XC.srate,30,[],opts.B,opts.A);        

        % Eigenvector decomposition
        [VC,DC]             = eig(DefaultState.M);
        [DCsorted, idx]     = sort(diag(DC), 'descend');

        % Wrong way
        M = cov(XC.data');
        [VC_,DC_]           = eig(M);
        [DCsorted_, idx_]   = sort(diag(DC_), 'descend');    
                

                %% Figure 2
                figure; set(gcf, 'Color', 'w')        
        
                % Layout
                nrow=2; ncol=2;
                layout = tiledlayout(nrow, ncol, ...
                    'TileSpacing', 'compact', ...
                    'Padding', 'compact');       
        
                % Filter applied
                nexttile
                [h, f] = freqz(State45.B, State45.A, 1024, EEG.srate);
                plot(f, abs(h)*100, 'LineWidth', 2);
                grid on;
                ylabel('Magnitude (%)');
                xtickformat('%gHz')
                xlim([0 EEG.srate/2]);
                xticks([4, 12, 20:20:100])   
                title('Filter response')
        
                % Covariance matrix
                I3=nexttile
                imagesc(State45.M)
                set(gca, 'YDir', 'normal')
                xlabel('Channels')
                ylabel('Channels')
                colorbar()
                title(sprintf('Covariance matrix\n(as computed by ASR)'))  
        
                % Eigenvector
                I3=nexttile
                imagesc(VC(:, idx));
                set(gca, 'YDir', 'normal')
                xlabel('Principal component')
                ylabel('Channels')
                colorbar()
                title(sprintf('Eigenvector \n(Channel weights on principal components)'))
                colormap(I3, colors.maps.div)
        
                % Eigenvalue
                I3=nexttile
                imagesc(DC(:, idx));
                set(gca, 'YDir', 'normal')
                xlabel('Principal component')
                ylabel('Channels')
                colorbar()
                title(sprintf('Eigenvector \n(Explained variance by principal components)'))
                colormap(I3, colors.maps.lin) 

                

                %% Figure 3
                figure; set(gcf, 'Color', 'w')
        
                % Layout
                nrow=2; ncol=2;
                layout = tiledlayout(nrow, ncol, ...
                    'TileSpacing', 'compact', ...
                    'Padding', 'compact');
        
                % Actual covariance matric
                nexttile
                imagesc(M)
                set(gca, 'YDir', 'normal')
                xlabel('Channels')
                ylabel('Channels')
                colorbar()      
                title(sprintf('Covariance matrix\n(as computed by cov() )'))
                
                % Covariance matrix how ASR computes it
                I3=nexttile
                imagesc(State45.M)
                set(gca, 'YDir', 'normal')
                xlabel('Channels')
                ylabel('Channels')
                colorbar()
                title(sprintf('Covariance matrix\n(as computed by ASR)'))
                % colormap(I3, colors.maps.lin)
                % caxis([0 max(State.M, [], 'all')])
        
                % Eigenvalue
                I3=nexttile
                imagesc(VC_(:, idx_));
                set(gca, 'YDir', 'normal')
                xlabel('Principal component')
                ylabel('Channels')
                colorbar()
                title(sprintf('Eigenvector \n(as computed by cov())'))
                colormap(I3, colors.maps.div)          
                caxis([-.7 .7])
        
                % Eigenvector
                I3=nexttile
                imagesc(VC(:, idx));
                set(gca, 'YDir', 'normal')
                xlabel('Principal component')
                ylabel('Channels')
                colorbar()
                title(sprintf('Eigenvector \n(as computed by ASR)'))
                colormap(I3, colors.maps.div)
                caxis([-.7 .7])













        % Project data into PC space
        XP = XC;
        XP.data = VC(:, idx)' * XP.data;

        % Recounstruct
        
        %%
        windowlength = 1;         
        % round(EEG.srate/3)
        [Y30, outstate, MetaR] = plots.asr.asr_process(X.data, EEG.srate, State30, windowlength, [], []);
        [Y20, ~, MetaR] = plots.asr.asr_process(X.data, EEG.srate, State20, windowlength, [], []);
        [Y45, ~, MetaR] = plots.asr.asr_process(X.data, EEG.srate, State45, windowlength, [], []);
        
        % shift signal content back (to compensate for processing delay)
        Y30(:, 1:size(outstate.carry,2)) = [];
        Y30(:, end+1:end+size(outstate.carry,2)) = 0;
        Y20(:, 1:size(outstate.carry,2)) = [];
        Y20(:, end+1:end+size(outstate.carry,2)) = 0;
        Y45(:, 1:size(outstate.carry,2)) = [];
        Y45(:, end+1:end+size(outstate.carry,2)) = 0;


        istep = 951;
        % do a PCA to find potential artifact components
        [V,D] = eig(MetaR.Xcov(:,:,istep));
        % [D, order] = sort(diag(D), 'descend'); V = V(:,order);
        [D,order] = sort(reshape(diag(D),1,EEG.nbchan), 'descend'); 
        V = V(:,order);



    


        %%% Figure 4
        %%
        figure; set(gcf, 'Color', 'w')

        % Layout
        nrow=4; ncol=1;
        layout = tiledlayout(nrow, ncol, ...
            'TileSpacing', 'compact', ...
            'Padding', 'compact');

        %%% Project data into PC space
        nexttile
        Times_min = linspace(0, X.pnts/X.srate/60, X.pnts);    
        Times_min = Times_min(mask);
        plot(Times_min, XP.data', '.', 'MarkerSize', 2)
        xtickformat('%gmin')
        ylabel('\muV')
        xlim([Times_min(1) Times_min(end)])
        title(sprintf('Calibration data in PC space'))
  
        %%% RMS values
        nexttile
        Times_min = interp1(1:numel(Times_min), Times_min, linspace(1, numel(Times_min), size(State45.RMS, 2)), 'linear');
        
        % Times_min = resample(Times_min, numel(Times_min), size(State.RMS, 2));
        % Times_min = linspace(0, X.pnts/X.srate/60, size(State.RMS, 2));        
        plot(Times_min, State45.RMS', '.', 'MarkerSize',3 )
        xlim([Times_min(1) Times_min(end)])
        xtickformat('%gmin')
        ylabel('\muV')        
        title(sprintf('RMS values in PC space (%.1fs windows, %.2f%% overlap)', State45.Windowlen, State45.Overlap))
 
        %%% Means
        nexttile
        % Times_min = linspace(0, X.pnts/X.srate/60, size(State.RMS, 2));        
        % plot(State.Mu', '.', 'MarkerSize',3 )
        errorbar(numel(State45.Mu):-1:1, State45.Mu, State45.SD, '.', ...
            'MarkerSize', 12, ...
            'LineWidth', 1, ...
            'LineStyle', 'none', ...
            'CapSize', 10);        
        xlim([0 EEG.nbchan+1])
        xlabel('Channel (in PC space)')
        ylabel('\muV')     
        title(sprintf('Mean%cSTD RMS value per channel in PC space', 177))

        %%% T
        nexttile
        hold on
        plot(flip(diag(State30.TRaw)), '-', 'LineWidth', 2, 'DisplayName', 'k=30')
        plot(flip(diag(State45.TRaw)), '-', 'LineWidth', 2, 'DisplayName', 'k=45')
        xlabel('Channel (in PC space)')
        ylabel('\muV')
        title(sprintf('Threshold value per channel in PC space (Mean+k*STD)'))  
        xlim([0 EEG.nbchan+1])
        legend()     

 %%
        %%% Figure 5
        figure; set(gcf, 'Color', 'w')

        % Layout
        nrow=5; ncol=1;
        layout = tiledlayout(nrow, ncol, ...
            'TileSpacing', 'compact', ...
            'Padding', 'compact');

        %%% Data of one bout
        nexttile    
        Times_min = linspace(0, X.pnts/X.srate/60, X.pnts);    
        plot(Times_min, X.data')
        xtickformat('%gmin')
        ylabel('\muV')
        xlim([Times_min(1) Times_min(end)])
        title(sprintf('Selected %dmin bout', opts.BoutDuration))  

        Start = round((istep*32-31)/EEG.srate);
        End   = Start + MetaR.Windowlen;
        % xline([Start/60 End/60], 'k-', 'LineWidth', 2)

        % Draw a very thin patch
        y_lims = ylim; % Get current axis limits

        p = patch([Start/59 End/61 End/61 Start/59], ...
                  [y_lims(1) y_lims(1) y_lims(2) y_lims(2)], 'r');    
        set(p, 'FaceAlpha', 0.2, 'EdgeColor', 'none');        
       
        %%% Data of window
        nexttile    
        Times_s = linspace(Start, End, MetaR.Windowlen*EEG.srate);    
        plot(Times_s, X.data(:, Start*EEG.srate:End*EEG.srate-1)')
        xtickformat('%gs')
        ylabel('\muV')
        xlim([Times_s(1) Times_s(end)])
        title(sprintf('Selected window of %.1fs (updated every %d samples)', MetaR.Windowlen, MetaR.Stepsize))        
                  
        %%%
        I3=nexttile;
        % Times_s = interp1(1:numel(Times_s), Times_s, linspace(1, numel(Times_s), size(MetaR.Xcov)), 'linear');        
        imagesc(V)
        xlabel('Principal component')
        ylabel('Channel')
        colormap(I3, colors.maps.div)
        cbar=colorbar
        ylabel(cbar, 'Channel weights')
        title('Eigenvectors in selected window')

        %%% Tresholds and actually explained variance
        nexttile
        hold on
        plot(sum((State45.T*V).^2), 'DisplayName', 'k=45', 'LineWidth',2)
        plot(sum((State30.T*V).^2), 'DisplayName', 'k=30', 'LineWidth',2)
        plot(sum((State20.T*V).^2), 'DisplayName', 'k=20', 'LineWidth',2)
        plot(D, 'k-.', 'LineWidth',2, 'DisplayName', 'Actually explain variance')
        xlabel('Principal component')
        ylabel('Explained variance')
        legend('Location', 'Best')
        title('Remove artefactual components')
        xlim([1 20])

        %%% Cleaned window
        nexttile    
        Times_s = linspace(Start, End, MetaR.Windowlen*EEG.srate);    
        plot(Times_s, Y30(:, Start*EEG.srate:End*EEG.srate-1)')
        % Times_min = linspace(0, X.pnts/X.srate/60, X.pnts);    
        % plot(Times_min, Y')
        xtickformat('%gs')
        ylabel('\muV')
        xlim([Times_s(1) Times_s(end)])
        title(sprintf('Selected window cleaned (k=30)'))               


        %%
        %%% Figure 6
        figure; set(gcf, 'Color', 'w')

        % Layout
        nrow=4; ncol=1;
        layout = tiledlayout(nrow, ncol, ...
            'TileSpacing', 'compact', ...
            'Padding', 'compact');

        %%% Data of all epochs of respective sleep stage
        T1=nexttile;     
        Times_min = linspace(0, X.pnts/X.srate/60, X.pnts);    
        plot(Times_min, X.data')
        xtickformat('%gmin')
        ylabel('\muV')
        xlim([Times_min(1) Times_min(end)])
        title(sprintf('Selected %dmin bout', opts.BoutDuration))
        
        title(sprintf('Original data'))
        
        %%% K=45
        T2=nexttile;     
        plot(Times_min, Y45')
        xtickformat('%gmin')
        ylabel('\muV')
        xlim([Times_min(1) Times_min(end)])
        title(sprintf('k=45'))
  
        %%% K=30
        T3=nexttile;     
        plot(Times_min, Y30')
        xtickformat('%gmin')
        ylabel('\muV')
        xlim([Times_min(1) Times_min(end)])
        title(sprintf('k=30'))   

        %%% K=30
        T4=nexttile;     
        plot(Times_min, Y20')
        xtickformat('%gmin')
        ylabel('\muV')
        xlim([Times_min(1) Times_min(end)])
        title(sprintf('k=20'))   
        linkaxes([T1, T2, T3, T4])        


        %%
        %%% Figure 7
        figure; set(gcf, 'Color', 'w')

        % Layout
        nrow=4; ncol=1;
        layout = tiledlayout(nrow, ncol, ...
            'TileSpacing', 'compact', ...
            'Padding', 'compact');

        %%% Data of all epochs of respective sleep stage
        T1=nexttile;     
        Times_min = linspace(0, X.pnts/X.srate/60, X.pnts);    
        plot(Times_min, X.data')
        xtickformat('%gmin')
        ylabel('\muV')
        xlim([Times_min(1) Times_min(end)])
        title(sprintf('Selected %dmin bout', opts.BoutDuration))
        
        title(sprintf('Original data'))
        
        %%% K=45
        T2=nexttile;     
        plot(Times_min, X.data'-Y45')
        xtickformat('%gmin')
        ylabel('\muV')
        xlim([Times_min(1) Times_min(end)])
        title(sprintf('k=45'))
  
        %%% K=30
        T3=nexttile;     
        plot(Times_min, X.data'-Y30')
        xtickformat('%gmin')
        ylabel('\muV')
        xlim([Times_min(1) Times_min(end)])
        title(sprintf('k=30'))   

        %%% K=30
        T4=nexttile;     
        plot(Times_min, X.data'-Y20')
        xtickformat('%gmin')
        ylabel('\muV')
        xlim([Times_min(1) Times_min(end)])
        title(sprintf('k=20'))   
        linkaxes([T1, T2, T3, T4])                
        
    end
end