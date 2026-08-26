% EEG_INTERP - interpolate data channels
%
% Usage: EEGOUT = eeg_interp(EEG, badchans, method, t_range);
%
% Inputs: 
%     EEG      - EEGLAB dataset
%     badchans - [integer array] indices of channels to interpolate.
%                For instance, these channels might be bad.
%                [chanlocs structure] channel location structure containing
%                either locations of channels to interpolate or a full
%                channel structure (missing channels in the current 
%                dataset are interpolated).
%     method   - [string] method used for interpolation (default is 'spherical').
%                'invdist'/'v4' uses inverse distance on the scalp
%                'spherical' uses superfast spherical interpolation. 
%                'sphericalKang' uses Kang et al. 2015 parameters 
%                'spacetime' uses griddata3 to interpolate both in space 
%                and time (very slow and cannot be interrupted).
%     t_range  - [integer array with just two elements] time interval of the
%                badchans which should be interpolated. First element is
%                the start time and the second element is the end time.
%     params   - [1x3 array] Parameters [lambda m n]. Defaults are
%                'spherical'     -> lambda=0, m = 4, n =7 (Perin et al., 1989)
%                'sphericalKang' -> lambda=1e-8, m = 3 , n = 50 (Kang et al., 2015)
%                if the 'params' input is used, then a spherical
%                interpolation is used with these parameters (method is
%                ignored). See also https://github.com/sccn/eeglab/issues/500
% Output: 
%     EEGOUT   - data set with bad electrode data replaced by
%                interpolated data
%
% Note: See also the sphericalSplineInterpolate function of the clean_rawdata
%       EEGLAB plugin (in the private folder) which approximate Legendre 
%       polynomials with up to 500 iteration.                  .
%
% References:
% Perrin, F., Pernier, J., Bertrand, O., & Echallier, J. F. (1989). Spherical 
% splines for scalp potential and current density mapping. Electroencephalography 
% and clinical neurophysiology, 72(2), 184–187. https://doi.org/10.1016/0013-4694(89)90180-6
%
% Kang, S.S., Lano, T.J. & Sponheim, S.R. Distortions in EEG interregional phase 
% synchrony by spherical spline interpolation: causes and remedies. Neuropsychiatr 
% Electrophysiol 1, 9 (2015). https://doi.org/10.1186/s40810-015-0009-5
%
% Author: Arnaud Delorme, CERCO, CNRS, Mai 2006-
%
% -------------------------------------------------------------------------
% LOCAL PATCH (multivariate-prep/patches) - memory footprint only, same output.
%
% Stock eeg_interp interpolates a full-length continuous recording in one go, and
% the spherical branch is what makes that expensive. For a 256-channel night at
% 250 Hz (~7e6 samples, ~7 GB in single) the original code allocates, on top of the
% two copies of the recording it already holds:
%
%   * the mean-free copy of the data, and again the [values;zeros] copy   (~2x data)
%   * the spline coefficients C, one row per good channel                 (~1x data)
%   * inside the per-bad-channel loop, repmat(Gsph(j,:)', [1 nsamples])
%     and the C.*repmat product - both DOUBLE and full length             (~4x data,
%     re-allocated for every interpolated channel)
%   * tmpdata = zeros(nchan, pnts, trials), i.e. DOUBLE regardless of the class of
%     the data                                                          (~2x data),
%     and with EEGLAB's option_single on, eeg_checkset then casts that whole double
%     array back down to single - one more full-size copy. With option_single off
%     the recording simply stays double for the rest of the caller.
%
% That peaks at ~50 GB for a night's recording, which is where the 97 % RAM came from.
%
% Two changes here:
%
%   1. The spherical spline is applied as the linear operator it is. The whole
%      chain (subtract the channel mean, solve for C, weight by Gsph, add the mean
%      back) is linear in the data, so it collapses into one nBad x nGood matrix W
%      that is built once from the electrode coordinates - see spheric_spline_matrix.
%      Applying it is a single matrix product, run in time chunks of ~256 MB, written
%      straight into the output array. No full-length temporary survives a chunk.
%   2. The output array is allocated in the class of EEG.data ('like'), so single
%      data stays single instead of being silently promoted to double.
%
% Result is the same interpolation (same Gelec/Gsph/pinv, differences at single
% precision rounding level) at ~2x the data size instead of ~7x. The non-spherical
% methods ('invdist'/'v4'/'spacetime') are untouched and still go through
% badchansdata. The original spheric_spline subfunction is gone - nothing called it
% any more.
% -------------------------------------------------------------------------

% Copyright (C) Arnaud Delorme, CERCO, 2006, arno@salk.edu
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function EEG = eeg_interp(ORIEEG, bad_elec, method, t_range, params)

    if nargin < 2
        help eeg_interp;
        return;
    end
    EEG = ORIEEG;
    
    if nargin < 3
        disp('Using spherical interpolation');
        method = 'spherical';
    end
    
    if nargin < 4
        t_range = [ORIEEG.xmin ORIEEG.xmax];
    end
    if nargin < 5
        if strcmpi(method, 'spherical')
            params = [0 4 7];
        elseif strcmpi(method, 'sphericalKang')
            params = [1e-8 3 50];
        elseif strcmpi(method, 'sphericalCRD')
            params = [1e-5 4 500];
        end
    else
        if length(params) ~=3 
            error('parameters array must have 3 elements')
        end
        method = 'spherical';
    end

    % check channel structure
    tmplocs = ORIEEG.chanlocs;
    if isempty(tmplocs) || isempty([tmplocs.X])
        error('Interpolation require channel location');
    end
    
    if isstruct(bad_elec)
       
        % remove duplicate channel labels
        % -------------------------------
        allLabels = { bad_elec.labels };
        if length(unique(allLabels)) ~= length(allLabels)
            [~,uniqueInd] = unique(allLabels);
            bad_elec = bad_elec(uniqueInd);
        end

        % add missing channels in interpolation structure
        % -----------------------------------------------
        lab1 = { bad_elec.labels };
        tmpchanlocs = EEG.chanlocs;
        lab2 = { tmpchanlocs.labels };
        [~, tmpchan] = setdiff_bc( lab2, lab1);
        tmpchan = sort(tmpchan);
        
        % From 'bad_elec' using only fields present on EEG.chanlocs
        fields = fieldnames(bad_elec);
        [~, indx1] = setxor(fields,fieldnames(EEG.chanlocs)); clear tmp;
        if ~isempty(indx1)
            bad_elec = rmfield(bad_elec,fields(indx1));
            fields = fieldnames(bad_elec);
        end
        
        if ~isempty(tmpchan)
            newchanlocs = [];
            for index = 1:length(fields)
                if isfield(bad_elec, fields{index})
                    for cind = 1:length(tmpchan)
                        fieldval = getfield( EEG.chanlocs, { tmpchan(cind) },  fields{index});
                        newchanlocs = setfield(newchanlocs, { cind }, fields{index}, fieldval);
                    end
                end
            end
            newchanlocs(end+1:end+length(bad_elec)) = bad_elec;
            bad_elec = newchanlocs;
        end
        if length(EEG.chanlocs) == length(bad_elec), return; end
        
        lab1 = { bad_elec.labels };
        tmpchanlocs = EEG.chanlocs;
        lab2 = { tmpchanlocs.labels };
        [~, badchans] = setdiff_bc( lab1, lab2);
        fprintf('Interpolating %d channels...\n', length(badchans));
        if isempty(badchans), return; end
        goodchans      = sort(setdiff(1:length(bad_elec), badchans));
       
        % re-order good channels
        % ----------------------
        [~, tmp2, neworder] = intersect_bc( lab1, lab2 );
        [~, ordertmp2] = sort(tmp2);
        neworder = neworder(ordertmp2);
        EEG.data = EEG.data(neworder, :, :);

        % looking at channels for ICA
        % ---------------------------
        %[tmp sorti] = sort(neworder);
        %{ EEG.chanlocs(EEG.icachansind).labels; bad_elec(goodchans(sorti(EEG.icachansind))).labels }
        
        % update EEG dataset (add blank channels)
        % ---------------------------------------
        if ~isempty(EEG.icasphere)
            
            [~, sorti] = sort(neworder);
            EEG.icachansind = sorti(EEG.icachansind);
            EEG.icachansind = goodchans(EEG.icachansind);
            EEG.chaninfo.icachansind = EEG.icachansind;
            
            % TESTING SORTING
            %icachansind = [ 3 4 5 7 8]
            %data = round(rand(8,10)*10)
            %neworder = shuffle(1:8)
            %data2 = data(neworder,:)
            %icachansind2 = sorti(icachansind)
            %data(icachansind,:)
            %data2(icachansind2,:)
        end
        % { EEG.chanlocs(neworder).labels; bad_elec(sort(goodchans)).labels }
        %tmpdata                  = zeros(length(bad_elec), size(EEG.data,2), size(EEG.data,3));
        %tmpdata(goodchans, :, :) = EEG.data;
        
        % looking at the data
        % -------------------
        %tmp1 = mattocell(EEG.data(sorti,1));
        %tmp2 = mattocell(tmpdata(goodchans,1));
        %{ EEG.chanlocs.labels; bad_elec(goodchans).labels; tmp1{:}; tmp2{:} }
        %EEG.data      = tmpdata;
        
        EEG.chanlocs  = bad_elec;

    else
        badchans  = bad_elec;
        if iscell(badchans)
            badchans = eeg_chaninds(EEG, badchans);
        end
        goodchans = setdiff_bc(1:EEG.nbchan, badchans);
        if strcmpi(method, 'sphericalfast')
            EEG.data(badchans,:) = [];
            EEG.nbchan = length(goodchans);
        else
            oldelocs  = EEG.chanlocs;
            EEG       = pop_select(EEG, 'nochannel', badchans);
            EEG.chanlocs = oldelocs;
            disp('Interpolating missing channels...');
        end
    end

    % find non-empty good channels
    % ----------------------------
    origoodchans = goodchans;
    chanlocs     = EEG.chanlocs;
    nonemptychans = find(~cellfun('isempty', { chanlocs.theta }));
    [tmp indgood ] = intersect_bc(goodchans, nonemptychans);
    goodchans = goodchans( sort(indgood) );
    datachans = getdatachans(goodchans,badchans);
    badchans  = intersect_bc(badchans, nonemptychans);
    if isempty(badchans), return; end
    
    % scan data points
    % ----------------
    splineW = [];   % local patch: set by the spherical branch, see the header
    if strcmpi(method, 'spherical') || strcmpi(method, 'sphericalfast') || strcmpi(method, 'sphericalKang') || strcmpi(method, 'sphericalCRD')
        % get theta, rad of electrodes
        % ----------------------------
        tmpgoodlocs = EEG.chanlocs(goodchans);
        xelec = [ tmpgoodlocs.X ];
        yelec = [ tmpgoodlocs.Y ];
        zelec = [ tmpgoodlocs.Z ];
        rad = sqrt(xelec.^2+yelec.^2+zelec.^2);
        xelec = xelec./rad;
        yelec = yelec./rad;
        zelec = zelec./rad;
        tmpbadlocs = EEG.chanlocs(badchans);
        xbad = [ tmpbadlocs.X ];
        ybad = [ tmpbadlocs.Y ];
        zbad = [ tmpbadlocs.Z ];
        rad = sqrt(xbad.^2+ybad.^2+zbad.^2);
        xbad = xbad./rad;
        ybad = ybad./rad;
        zbad = zbad./rad;
        if isempty(xbad)
            error('Trying to interpolate a channel without coordinates assigned to it');
        end

        % Local patch: build the spline operator now, apply it further down straight
        % into the output array. splineW being non-empty is what marks this branch.
        splineW      = spheric_spline_matrix( xelec, yelec, zelec, xbad, ybad, zbad, params);
        badchansdata = [];
    elseif strcmpi(method, 'spacetime') % 3D interpolation, works but x10 times slower
        disp('Warning: if processing epoch data, epoch boundary are ignored...');
        disp('3-D interpolation, this can take a long (long) time...');
        tmpgoodlocs = EEG.chanlocs(goodchans);
        tmpbadlocs = EEG.chanlocs(badchans);
        [xbad ,ybad]  = pol2cart([tmpbadlocs.theta],[tmpbadlocs.radius]);
        [xgood,ygood] = pol2cart([tmpgoodlocs.theta],[tmpgoodlocs.radius]);
        pnts = size(EEG.data,2)*size(EEG.data,3);
        zgood = 1:pnts;
        zgood = repmat(zgood, [length(xgood) 1]);    
        zgood = reshape(zgood,numel(zgood),1);
        xgood = repmat(xgood, [1 pnts]); xgood = reshape(xgood,numel(xgood),1);
        ygood = repmat(ygood, [1 pnts]); ygood = reshape(ygood,numel(ygood),1);
        tmpdata = reshape(EEG.data(datachans,:,:), numel(EEG.data(datachans,:,:)),1);
        zbad = 1:pnts;
        zbad = repmat(zbad, [length(xbad) 1]);     
        zbad = reshape(zbad,numel(zbad),1);
        xbad = repmat(xbad, [1 pnts]); xbad = reshape(xbad,numel(xbad),1);
        ybad = repmat(ybad, [1 pnts]); ybad = reshape(ybad,numel(ybad),1);
        badchansdata = griddatan([ygood, xgood, zgood], tmpdata,[ybad, xbad, zbad], 'nearest'); % interpolate data                                            
    else 
        % get theta, rad of electrodes
        % ----------------------------
        tmpchanlocs = EEG.chanlocs;
        [xbad ,ybad]  = pol2cart([tmpchanlocs( badchans).theta],[tmpchanlocs( badchans).radius]);
        [xgood,ygood] = pol2cart([tmpchanlocs(goodchans).theta],[tmpchanlocs(goodchans).radius]);

        fprintf('Points (/%d):', size(EEG.data,2)*size(EEG.data,3));
        badchansdata = zeros(length(badchans), size(EEG.data,2)*size(EEG.data,3));
        
        for t=1:(size(EEG.data,2)*size(EEG.data,3)) % scan data points
            if mod(t,100) == 0, fprintf('%d ', t); end
            if mod(t,1000) == 0, fprintf('\n'); end;          
        
            %for c = 1:length(badchans)
            %   [h EEG.data(badchans(c),t)]= topoplot(EEG.data(goodchans,t),EEG.chanlocs(goodchans),'noplot', ...
            %        [EEG.chanlocs( badchans(c)).radius EEG.chanlocs( badchans(c)).theta]);
            %end
            tmpdata = reshape(EEG.data, size(EEG.data,1), size(EEG.data,2)*size(EEG.data,3) );
            if strcmpi(method, 'invdist'), method = 'v4'; end
            [Xi,Yi,badchansdata(:,t)] = griddata(ygood, xgood , double(tmpdata(datachans,t)'),...
                                                    ybad, xbad, method); % interpolate data                                            
        end
        fprintf('\n');
    end
    
    % Local patch: allocate in the class of the data (was double, which promoted every
    % single-precision recording) and keep it 2-D until the very end, so the spherical
    % branch can fill it chunk by chunk.
    % (stock relied on tmpdata growing by out-of-range assignment when bad_elec is a
    % numeric index list; size it up front instead so the reshape below has a number)
    if isstruct(bad_elec)
        nchanout = length(bad_elec);
    else
        nchanout = max([origoodchans(:); badchans(:)]);
    end
    tmpdata                  = zeros(nchanout, EEG.pnts*EEG.trials, 'like', EEG.data);
    tmpdata(origoodchans, :) = EEG.data(:,:);

    if ~isempty(splineW)
        % The good channels now live in tmpdata, so the working copy inside this
        % function can go before the interpolated channels are computed - that is what
        % keeps the peak at roughly two copies of the recording.
        EEG.data = [];
        splineW  = cast(splineW, 'like', tmpdata);
        nsamples = size(tmpdata,2);
        if isa(tmpdata, 'double'), bytesper = 8; else, bytesper = 4; end
        step     = max(1, floor(256*2^20 / (max(1,length(goodchans))*bytesper)));
        for i1 = 1:step:nsamples
            i2 = min(i1+step-1, nsamples);
            tmpdata(badchans, i1:i2) = splineW * tmpdata(goodchans, i1:i2);
        end
    else
        tmpdata(badchans, :) = reshape(badchansdata, length(badchans), []);
    end

    if ~isempty(t_range)
        t_start = t_range(1);
        t_end   = t_range(2);
        if t_start~=ORIEEG.xmin || t_end~=ORIEEG.xmax
            if EEG.trials == 1
                times_2b_ignored = [1:floor(t_start*ORIEEG.srate), floor(t_end*ORIEEG.srate):floor(ORIEEG.xmax*ORIEEG.srate)];
                tmpdata(badchans, times_2b_ignored) = ORIEEG.data(badchans, times_2b_ignored);
            end
        end
    end

    EEG.data = reshape(tmpdata, nchanout, EEG.pnts, EEG.trials);
    clear tmpdata;
    EEG.nbchan = size(EEG.data,1);
    if ~strcmpi(method, 'sphericalfast')
        EEG = eeg_checkset(EEG);
    end
    
% get data channels
% -----------------
function datachans = getdatachans(goodchans, badchans);
      datachans = goodchans;
      badchans  = sort(badchans);
      for index = length(badchans):-1:1
          datachans(find(datachans > badchans(index))) = datachans(find(datachans > badchans(index)))-1;
      end
        
% -----------------
% spherical splines
% -----------------
% function [x, y, z, Res] = spheric_spline_old( xelec, yelec, zelec, values)
% 
% SPHERERES = 20;
% [x,y,z] = sphere(SPHERERES);
% x(1:(length(x)-1)/2,:) = []; x = [ x(:)' ];
% y(1:(length(y)-1)/2,:) = []; y = [ y(:)' ];
% z(1:(length(z)-1)/2,:) = []; z = [ z(:)' ];
% 
% Gelec = computeg(xelec,yelec,zelec,xelec,yelec,zelec);
% Gsph  = computeg(x,y,z,xelec,yelec,zelec);
% 
% % equations are 
% % Gelec*C + C0  = Potential (C unknown)
% % Sum(c_i) = 0
% % so 
% %             [c_1]
% %      *      [c_2]
% %             [c_ ]
% %    xelec    [c_n]
% % [x x x x x]         [potential_1]
% % [x x x x x]         [potential_ ]
% % [x x x x x]       = [potential_ ]
% % [x x x x x]         [potential_4]
% % [1 1 1 1 1]         [0]
% 
% % compute solution for parameters C
% % ---------------------------------
% meanvalues = mean(values); 
% values = values - meanvalues; % make mean zero
% C = pinv([Gelec;ones(1,length(Gelec))]) * [values(:);0];
% 
% % apply results
% % -------------
% Res = zeros(1,size(Gsph,1));
% for j = 1:size(Gsph,1)
%     Res(j) = sum(C .* Gsph(j,:)');
% end
% Res = Res + meanvalues;
% Res = reshape(Res, length(x(:)),1);

function W = spheric_spline_matrix( xelec, yelec, zelec, xbad, ybad, zbad, params)
% Local patch (see header): the nBad x nGood operator that the stock spheric_spline
% applies one full-length temporary at a time.
%
% With V the good-channel data (nGood x nSamples), mu = mean(V,1) and
% P = pinv([Gelec+lambda*I ; ones(1,nGood)]), the stock function computes
%
%   C      = P * [V - 1*mu ; 0]  = P(:,1:nGood) * (V - 1*mu)
%   allres = Gsph * C + 1*mu
%
% Every step is linear in V, and mu = (ones(1,nGood)/nGood) * V, so the whole thing
% is allres = W*V with W built once from the coordinates alone. Nothing here scales
% with the length of the recording.

nelec = length(xelec);
nbad  = length(xbad);

Gelec = computeg(xelec,yelec,zelec,xelec,yelec,zelec,params);
Gsph  = computeg(xbad ,ybad ,zbad ,xelec,yelec,zelec,params);

lambda = params(1);
P = pinv([Gelec+eye(nelec)*lambda; ones(1,nelec)]);

M = Gsph * P(:,1:nelec);                                  % applied to the mean-free data
W = M - (sum(M,2)/nelec)*ones(1,nelec) ...                % ... minus what the mean removed
      + ones(nbad,nelec)/nelec;                           % ... plus the mean added back

% compute G function
% ------------------
function g = computeg(x,y,z,xelec,yelec,zelec, params)

unitmat = ones(length(x(:)),length(xelec));
EI = unitmat - sqrt((repmat(x(:),1,length(xelec)) - repmat(xelec,length(x(:)),1)).^2 +... 
                (repmat(y(:),1,length(xelec)) - repmat(yelec,length(x(:)),1)).^2 +...
                (repmat(z(:),1,length(xelec)) - repmat(zelec,length(x(:)),1)).^2);

g = zeros(length(x(:)),length(xelec));
m = params(2);
maxn = params(3);

for n = 1:maxn % 200
    if ismatlab
        L = legendre(n,EI);
    else % Octave legendre function cannot process 2-D matrices
        for icol = 1:size(EI,2)
            tmpL = legendre(n,EI(:,icol));
            if icol == 1, L = zeros([ size(tmpL) size(EI,2)]); end
            L(:,:,icol) = tmpL;
        end
    end
    g = g + ((2*n+1)/(n^m*(n+1)^m))*squeeze(L(1,:,:));
end
g = g/(4*pi);    

