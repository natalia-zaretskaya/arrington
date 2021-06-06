function varargout = vpxPreprocess(varargin)
% analyzes eyetracker data with given settings
%
% use as:
% dataout = vpxPreprocess(datain, settings, plot_flag)
% or:
% dataout = vpxPreprocess
% with defaults
%
% input:
% e: a cell/structure array from vpxReadFile 
%   of size [nfiles, 1] for mono
%       and [nfiles, 2] for stereo tracking
% s: structure with fields
%     s.validMarker:     0 par-glint; 1 par
%     s.dpu:              degrees per eyetracker unit, [x y]
%     s.highpassFlag:    0/1, highpass filter the data to reduce e.g. motion
%       artifacts
%     s.highpassCutoffSecs:          lowpass cutoff in seconds
%     s.blinkThreshStd: blink threshold in standard deviations from the mean
%     s.s.removeBeforeBlinkSecs:            n samples to remove before the blink
%     s.s.removeAfterBlinkSecs:           n samples to remove after the blink
%     s.runningAverageWindowSecs:  window length in seconds for the running average filter

% PLOT_FLAG: 1 or 0

% output:
% DATAOUT: same as datain, but whit fields containing analysis results
% Natalia 17.06.2010
% updated 31.08.2011
% updated 20.03.2013 par width added
% updated 19.05.2014 various improvements, different input structure

set(0,'DefaultTextInterpreter', 'none')


if (nargin < 1)
    
    e                          = vpxReadFile;
    s.validMarker              = mode(e{1}.valid); % 0 par-glint; 1 par
    s.dpu                      = [30, 20]; % this is arbitrary
    s.highpassFlag             = 0;% careful, sometimes generates big artifacts in the timeseries (better substract the mean of each trial)
    s.blinkThreshStd           = 5; % median-based standard deviation of par width
    s.saccThreshStd            = 5; % median-based standard deviation of speed
    s.removeBeforeBlinkSecs    = 0.1;  % period in seconds to remove before the blink
    s.removeAfterBlinkSecs     = 0.1; % period in seconds samples to remove after the blink
    s.highpassCutoffSecs        = 10;% low cutoff in seconds
    s.runningAverageWindowSecs = 0.2; % running average window in seconds for smoothing x and y
    
    plot_flag = 1;
    
else
    e         = varargin{1};
    s         = varargin{2};
    plot_flag = varargin{3};
end


% --- internal eyetracker settings --- %
scr_eyetracker = [1 1]; % internal coordinate system of the eyetracker

% invalid sample ids:
% 6 - valid, but offscreen values
% 7 - blinks
% 5 - ? eyetracker invalid

% preprocessing
for i = 1:size(e,1)s
    for j = 1:size(e,2)
        if isempty(e{i,j});
            warning('no data in this session ');
            break
        end
        
        fprintf('Preprocessing: %s\n', e{i,j}.fname);
        
        s.sr = 1/mode(diff(e{i,j}.time));
        fprintf('determined sampling rate: %f Hz \n', s.sr);
        len = round(s.sr*s.runningAverageWindowSecs);
        
        % additional entries in the output structure
        e{i,j}.settings = s;
        e{i,j}.ix = []; % interpolated x
        e{i,j}.iy = []; % interpolated y
        e{i,j}.dist = []; % interpolated eucledian distance from the meadian
        e{i,j}.vx = []; % velocity in x dir
        e{i,j}.vy = []; %
        e{i,j}.v = []; % velocity
        e{i,j}.ipar = []; % par aspect ratio
        e{i,j}.ipw = []; % par width
        e{i,j}.blink_onsets = []; % indices of blinks
        e{i,j}.blink_offsets = []; % indices of blinks
        e{i,j}.sacc_onsets = []; % indices of detected saccades
        e{i,j}.sacc_offsets = []; % indices of detected saccades
        
        % - find valid data points (eyetracker found the par)
        valid = logical(e{i,j}.valid==s.validMarker);
        
        % 1. ----------- clean up samples ------------ %
        e{i,j}.par(isinf(e{i,j}.par)) = 1;
        
        % sometimes needed for good filtering:
        e{i,j}.x([1:5]) = nanmedian(e{i,j}.x(6:10));
        e{i,j}.y([1:5]) = nanmedian(e{i,j}.y(6:10));
        e{i,j}.x([end-5:end]) = nanmedian(e{i,j}.x(end-10:end-6));
        e{i,j}.y([end-5:end]) = nanmedian(e{i,j}.y(end-10:end-6));
        e{i,j}.par([1:5  end-5:end]) = nanmedian(e{i,j}.par);
        e{i,j}.pw([1:5  end-5:end]) = nanmedian(e{i,j}.pw);
        e{i,j}.valid([1:5  end-5:end]) = s.validMarker;
        
        % - center around zero
        e{i,j}.x = e{i,j}.x - nanmedian(e{i,j}.x(valid));
        e{i,j}.y = e{i,j}.y - nanmedian(e{i,j}.y(valid));
        
        % - find valid but impossible data points (outside of the screen)
        e{i,j}.valid(abs(e{i,j}.x)>scr_eyetracker(1)/2) = 6;
        e{i,j}.valid(abs(e{i,j}.y)>scr_eyetracker(2)/2) = 6;
        
        %  clean up eyetracker's dropped frames
%         tmp = strfind(e{i,j}.marker, 'Lost');
%         tmp = ~cellfun('isempty', tmp);
        e{i,j}.time = 0: length(e{i,j}.x)-1;
        e{i,j}.time= e{i,j}.time/s.sr;
        
        % 2. ---------- determine blinks: ----------- %
        %  median-based standard deviataion (less affected by outliers)
        %  formula: msdx = sqrt( median(x.^2) - (median(x))^2 );
        e{i,j}.par = e{i,j}.par - nanmedian(e{i,j}.par);
        msd = sqrt(   median(e{i,j}.par.^2) - (median(e{i,j}.par))^2   ) ;
        blink_thresh = msd*s.blinkThreshStd;
        
        % find blinks
        pupildev = abs(e{i,j}.par-nanmedian(e{i,j}.par));
        e{i,j}.blink_onsets = find(diff(pupildev>blink_thresh)>0); % find the starting indices of blinks
        e{i,j}.blink_offsets = find(diff(pupildev>blink_thresh)<0);  % find end indices of blinks
        
        % clean blinks
        [e{i,j}.blink_onsets, e{i,j}.blink_offsets] = cleanEvents(e{i,j}.blink_onsets, e{i,j}.blink_offsets);
        
        % include some samples before and after the blink
        e{i,j}.blink_onsets = e{i,j}.blink_onsets - round(s.removeBeforeBlinkSecs*s.sr);
        e{i,j}.blink_offsets = e{i,j}.blink_offsets + round(s.removeAfterBlinkSecs*s.sr);
        
        % after substracting s.removeBeforeBlinkSecs*s.sr or adding s.removeAfterBlinkSecs*s.sr there might be
        % non-existing indices, clean them:
        e{i,j}.blink_onsets(e{i,j}.blink_onsets<=0) = 1;
        e{i,j}.blink_offsets(e{i,j}.blink_offsets>length(e{i,j}.x)) = length(e{i,j}.x);
        
        % set data points between ons and offs of a blink to invalid
        if ~isempty(e{i,j}.blink_onsets)
            for blink_i = 1: length(e{i,j}.blink_onsets)
                e{i,j}.valid(e{i,j}.blink_onsets(blink_i):e{i,j}.blink_offsets(blink_i)) = 7;
            end
        end
        
        % define new onsets and offsets of blinks
        e{i,j}.blink_onsets = find(diff(e{i,j}.valid==7)>0); % find the starting indices of blinks
        e{i,j}.blink_offsets = find(diff(e{i,j}.valid==7)<0); % find the starting indices of blinks
        
        % clean again
        [e{i,j}.blink_onsets, e{i,j}.blink_offsets] = cleanEvents(e{i,j}.blink_onsets, e{i,j}.blink_offsets);
        
        % ----- interpolate blinks, extrapolate start and end values --- %
        % redefinde valid samples
        valid = logical(e{i,j}.valid==s.validMarker);
        
        % interpolate
        itmp = interp1( find(valid), e{i,j}.x(valid), find(~valid),  'linear', 'extrap'  );
        e{i,j}.ix = e{i,j}.x;
        e{i,j}.ix(~valid) = itmp;
        
        itmp = interp1( find(valid), e{i,j}.y(valid), find(~valid),  'linear', 'extrap'  );
        e{i,j}.iy = e{i,j}.y;
        e{i,j}.iy(~valid) = itmp;
        
        itmp = interp1( find(valid), e{i,j}.par(valid), find(~valid)',  'linear', 'extrap' );
        e{i,j}.ipar = e{i,j}.par;
        e{i,j}.ipar(~valid) = itmp;
        
        itmp = interp1( find(valid), e{i,j}.pw(valid), find(~valid)',  'linear', 'extrap' );
        e{i,j}.ipw = e{i,j}.pw;
        e{i,j}.ipw(~valid) = itmp;
        
        % 3. ------------- further preprocessing ---------------- %
        % center data around zero
        e{i,j}.ix = detrend(e{i,j}.ix);
        e{i,j}.iy = detrend(e{i,j}.iy);
        e{i,j}.ix = e{i,j}.ix- nanmedian(e{i,j}.ix(valid));
        e{i,j}.iy = e{i,j}.iy- nanmedian(e{i,j}.iy(valid));
        
        % smooth with running average to suppress noise
        e{i,j}.ix = filter(ones(1,len)/len,1, e{i,j}.ix);
        e{i,j}.iy = filter(ones(1,len)/len,1, e{i,j}.iy);
        
        % high-pass filter
        if s.highpassFlag == 1% better if only par was tracked
            norm_lcutoff = (1/s.highpassCutoffSecs)/(sr/2);
            [b,a] = butter(1, norm_lcutoff, 'high');
            e{i,j}.ix  = filtfilt(b,a,e{i,j}.ix);
            e{i,j}.iy  = filtfilt(b,a,e{i,j}.iy);
        end
        
        % transfer fom eyetracker units to degrees
        e{i,j}.ix = e{i,j}.ix.*s.dpu(1);
        e{i,j}.iy = e{i,j}.iy.*s.dpu(2);
        
        % velocity on interpolated data
        [ e{i,j}.vx]=gradient(e{i,j}.ix, 1/s.sr);
        [ e{i,j}.vy]=gradient(e{i,j}.iy, 1/s.sr);
        
        % smooth velocity to suppress noise
        e{i,j}.vx = filter(ones(1,len)/len,1, e{i,j}.vx);
        e{i,j}.vy = filter(ones(1,len)/len,1, e{i,j}.vy);
        
        e{i,j}.v=sqrt(e{i,j}.vx.^2+e{i,j}.vy.^2);
        
        % physical distance from fixation
        e{i,j}.dist = sqrt((e{i,j}.ix.^2 + e{i,j}.iy.^2));
        
        % -------------saccade detection  ----------------------- %
        % (based on the idea of Engbert, R. & Kliegl, R. (2002) )
        % 1. Computing thresholds
        msdx = sqrt( median(e{i,j}.vx.^2) - (median(e{i,j}.vx))^2 );
        msdy = sqrt( median(e{i,j}.vy.^2) - (median(e{i,j}.vy))^2 );
        radiusx = s.saccThreshStd*msdx;
        radiusy = s.saccThreshStd*msdy;
        
        % 2. detect saccades
        test = (e{i,j}.vx/radiusx).^2 + (e{i,j}.vy/radiusy).^2;
        e{i,j}.sacc_onsets = find(diff(test>1)>0);
        e{i,j}.sacc_offsets = find(diff(test>1)<0);
        
        % 3. clean saccades
        [e{i,j}.sacc_onsets, e{i,j}.sacc_offsets] = cleanEvents(e{i,j}.sacc_onsets, e{i,j}.sacc_offsets);
        
        % 4. get saccade amplitude and duration
        if ~isempty(e{i,j}.sacc_onsets) && ~isempty(e{i,j}.sacc_offsets)
            e{i,j}.sacc_amp = zeros(size(e{i,j}.sacc_onsets));
            e{i,j}.sacc_dur = zeros(size(e{i,j}.sacc_onsets));
            for k = 1:length(e{i,j}.sacc_onsets)
                e{i,j}.sacc_amp(k) = max(e{i,j}.dist(e{i,j}.sacc_onsets(k):e{i,j}.sacc_offsets(k)));
                e{i,j}.sacc_dur(k) = e{i,j}.time(e{i,j}.sacc_offsets(k)) - e{i,j}.time(e{i,j}.sacc_onsets(k));
            end
        else
            fprintf('no saccades found in this file')
        end
        
        % ---- extract user-specified conditions, if any ---- %
        e{i,j}.cond_ons = zeros(size(e{i,j}.time));
        idx = ~cellfun(@isempty, e{i,j}.marker); % get indices of non-empty cells
        condnames = unique(e{i,j}.marker(idx));
        for c = 1:length(condnames)
            foundstring = strfind(e{i,j}.marker, condnames{c});
            c_idx = ~cellfun(@isempty, foundstring);
            e{i,j}.cond_ons(c_idx) = c;
        end
        e{i,j}.condnames = condnames;
        
    end % for j
end % for i

varargout{1} = e;

%%
if plot_flag
    for i = 1:size(e,1)
        for j = 1:size(e,2)
            
            % ------------------------- PLOTS ------------------------- %
            c = 1;
            
            valid = logical(e{i,j}.valid==s.validMarker);
            figure
            
            % data quality and experimental conditions
            h = subplot(7,2,1:2);
            
            plot(e{i,j}.time, e{i,j}.valid, [ 'k.']); %
            hold on
            
            plot(e{i,j}.time(logical(e{i,j}.cond_ons)),...
                e{i,j}.cond_ons(logical(e{i,j}.cond_ons))+10, ['r.'])
            
            legend(h, {'quality', 'conditions'}, 1)
            ylabel('events')
            title(e{i,j}.fname);
            
            % raw x and y
            h = subplot(7,2,3:4);
            plot(e{i,j}.time, e{i,j}.x, 'g-', e{i,j}.time, e{i,j}.y, 'b-');
            hold on
            axis tight
            ylabel('raw x & y')
            legend(h, 'raw x', 'raw y', 1)
            
            % preprocessed x and y
            subplot(7,2,5:6)
            plot(e{i,j}.time, e{i,j}.ix, 'g-');
            hold on
            plot(e{i,j}.time, e{i,j}.iy, 'b-');
            axis tight
            ylabel('preproc x & y deg')
            legend(h, 'x', 'y', 1)
            axis tight
            ylim([-5 5])
            
            % velocity
            subplot(7,2,7:8)
            plot(e{i,j}.time, e{i,j}.v);
            hold on;
            for k = 1:length(e{i,j}.sacc_onsets)
                plot(e{i,j}.time(e{i,j}.sacc_onsets(k):e{i,j}.sacc_offsets(k)),...
                    e{i,j}.v(e{i,j}.sacc_onsets(k):e{i,j}.sacc_offsets(k)) , 'r')
            end
            ylabel('velocity')
            axis tight
            
            % par width
            subplot(7,2,9:10)
            plot(e{i,j}.time, e{i,j}.ipw, '-');
            hold on
            ylabel('pupil width')
            xlabel('sec')
            axis tight
            
            % par aspect ratio with interpolated blinks and blink thresholds
            subplot(7,2,11:12)
            plot(e{i,j}.time, e{i,j}.par, '-');
            hold on
            plot(e{i,j}.time, e{i,j}.ipar, 'g');
            ylabel('pupil aspect ratio')
            xlabel('Time (sec)')
            axis tight
            refline(0, nanmedian(e{i,j}.par)+blink_thresh)
            refline(0, nanmedian(e{i,j}.par)-blink_thresh)
            
            % resulting histogram of valid samples
            subplot(7,2,13)
            n_bins_x = ceil(length(e{i,j}.x(valid))^(1/3)*range(e{i,j}.x(valid))/(2*iqr(e{i,j}.x(valid)))); % Freedman-Diaconis rule (see Wikipedia)
            hist(e{i,j}.ix(valid), n_bins_x)
            hold on
            h = findobj(gca,'Type','patch');
            title('x (deg)')
            
            subplot(7,2,14)
            n_bins_y = ceil(length(e{i,j}.y(valid))^(1/3)*range(e{i,j}.y(valid))/(2*iqr(e{i,j}.y(valid)))); % Freedman-Diaconis rule (see Wikipedia)
            hist(e{i,j}.iy(valid), n_bins_y)
            hold on
            h = findobj(gca,'Type','patch');
            xlabel('screen position')
            title('y (deg)')
            
            % ---  --- %
            c = c+1;
        end % for j
    end % for i
end % if plot_flag


end


function [onsets, offsets] = cleanEvents(onsets, offsets)

if ~isempty(onsets) && ~isempty(offsets)
    % make sure the first is  onset and the last is offset
    if onsets(1)>offsets(1)
        offsets(1) = [];
    end
    if onsets(end)>offsets(end)
        onsets(end) = [];
    end
else
    disp('no events found' )
end

end