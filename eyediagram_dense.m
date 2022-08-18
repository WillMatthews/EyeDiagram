function [fh, scr_out] = eyediagram_dense(varargin)
% EYEDIAGRAM_DENSE Creates a more useful eye diagram compared to MATLAB's default.
% Uses different arguments as matlab's eyediagram function, but it's
% made to be close to a drop-in replacement.
%
% NOTE: Auto x axis label may be poor if your system suffers severe ISI.
% eyediagram_dense(signal, npoints)
%
% INPUTS:
%
% eyediagram_dense(signal, npoints, "time", time_vector)
%     makes auto x label use a time vector (same length as signal) to plot
%     the real time on the eye diagram
%
% eyediagram_dense(signal, npoints, "midbit", true)
%     forces auto x label to use middle of the bit (defined as midpoint of
%     zero crossings). Boolean input. Default false.
%
% eyediagram_dense(signal, npoints, "histogram", true)
%     plots histograms for x crossings and signal y level. Default false.
%
%
% RETURNS:
% [fh, scr_out] = eyediagram_dense(signal, npoints "time", time, "midbit")
%     where fh is a figure handle
%     scr_out is a struct with fields:
%        scr_out.x        vector the x axis scale (one el per el in matrix)
%        scr_out.y        vector the y axis scale (one el per el in matrix)
%        scr_out.image    image matrix after processing
%        scr_out.raw_data raw binned signal matrix

%% process args

signal = varargin{1};
rep    = varargin{2};
time   = [];
midbit = false;
plothist = false;

for i = 1:nargin
    if strcmpi("time", varargin{i})
        time = varargin{i+1};
    end
    if strcmpi("midbit", varargin{i})
        midbit = varargin{i+1};
    end
    if strcmpi("histogram", varargin{i})
        plothist = varargin{i+1};
    end
end



%% setup key variables
scr = 500; % Set screen size in x

num2delete = rep*mod(numel(signal)/rep,1);
num2delete = round(num2delete);
signal = signal((num2delete + 1):end); % Truncate signal so we don't have any half-traces on plot
sig_sliced = reshape(signal, rep, numel(signal)/rep)'; % slice up signal into 'rep' long sections
clear num2delete


screen_size = round(scr*[1,0.776018099547511]); % set screen dimensions (default aspect ratio)
screen = zeros(screen_size(1), screen_size(2)); % init image

% y (VOLT) axis bins
bins = linspace(min(signal),max(signal), screen_size(2));
yscale = bins;

% x (TIME) axis scale (not used - only for user)
if isempty(time)
    xscale = [];
else
    xscale = linspace(0, time(size(sig_sliced,2))-time(1), screen_size(1));
end

scr_out.x = xscale;
scr_out.y = yscale; % combine for output

%% volt binning
for i = 1:size(sig_sliced,2)
    slice = sig_sliced(i,:);
    % interpolate the signal slice - screen may be bigger or smaller than rep
    sliceinterp = interp1(linspace(0,1, rep),slice,linspace(0,1, screen_size(1)));
    for j = 1:numel(sliceinterp)
        [~,v] = min(abs(sliceinterp(j)-yscale)); % find closest voltage bin in slice
        screen(j,v) = screen(j,v) + 1;
    end
end

screen = screen'; % transpose the image
screen = flipud(screen);
scr_out.raw_data = screen;

%% image processing
image = screen/max(max(screen)); % scale to 0-1
image = sqrt(image); % gamma correction, gamma = 0.5
image = image/max(max(image)); % rescale to 0-1
%image = imgaussfilt(image, 1);
clear screen

%% plot code
fh = figure();
hEye = axes(fh);
if plothist
    hYHist = axes(fh);
    hXHist = axes(fh);
end
axes(hEye);


imshow(image);


scr_out.image = image;

%% auto-place y axis markers at BIT_ZERO and BIT_ONE
try
    % Histogram to find BIT_ZERO and BIT_ONE
    [vhist_counts, vhist_bin_edges] = histcounts(signal,screen_size(2)/2);
    vhist_bin_width = vhist_bin_edges(2)-vhist_bin_edges(1);
    vhist_bin_centres = vhist_bin_edges + vhist_bin_width/2;
    
    % Take second derivative of the low passed histogram counts
    % Minima in 2nd derivitive will be maxima on the histogram
    [~, idx1] = min(gradient(gradient(lowpass(vhist_counts,0.1))));
    v1 = vhist_bin_centres(idx1);
    halfway = round(screen_size(2)/4);
    
    % We are unsure if the peak in the histogram is a BIT_ZERO or a BIT_ONE
    % We should therefore truncate the half of the histogram we worked in, and
    % look for a maxima on the other half
    if idx1 < halfway
        vhist_counts = vhist_counts(halfway:end);
        ofs = halfway;
    else
        vhist_counts = vhist_counts(1:halfway);
        ofs = 0;
    end
    % Repeat process for other bit
    [~, idx2] = min(gradient(gradient(lowpass(vhist_counts,0.1))));
    v2 = vhist_bin_centres(idx2+ofs);
    idxpair = abs([idx1, idx2+ofs] - screen_size(2)/2);
catch
    warning("Error Encountered in Locating BIT ZERO and BIT ONE");
end

[~, lpFilt] = lowpass(sig_sliced(1,:), 0.1);
try
    %% Automatically find y=0 crossings for plot labelling
    zerocrossing_idxs = [];
    lp_crop_pt = ceil(length(lpFilt.Coefficients)/2);
    for i = 1:size(sig_sliced,1)
        % if we do not have a zero crossing - keep looking
        if ~((max(sig_sliced(i,:)) > 0) && (min(sig_sliced(i,:)) < 0))
            continue
        end
        % low pass signal - we have a zero crossing
        % add levels before and after to make length correct and remove
        % edge effects
        lpsig = filter(lpFilt, [repelem(sig_sliced(i,1), lp_crop_pt), ...
            sig_sliced(i,:), repelem(sig_sliced(i,end), lp_crop_pt)]);
        lpsig = lpsig((lp_crop_pt*2):end);
        
        if sig_sliced(i,1) > 0
            zerocrossing_idx = find(lpsig<0, 1); % locate position of zero cross pos to neg
        else
            zerocrossing_idx = find(lpsig>0, 1); % locate position of zero cross neg to pos
        end
        zerocrossing_idxs = [zerocrossing_idxs, zerocrossing_idx];
    end
    clear lpsig zerocrossing_idx
    
    % find distinct crossing points
    zerocrossing_idxs = sort(zerocrossing_idxs);
    gr = diff(zerocrossing_idxs);
    tpoints = find(gr > 25); % Detect discontinuities in crossing
    crossings = zeros(1, length(tpoints)+1);
    tp = [1, tpoints, length(zerocrossing_idxs)];
    for i = 1:length(crossings)
        crossings(i) = mean(zerocrossing_idxs((tp(i)+1):(tp(i+1)-1)));
    end
    clear gr tpoints tp i
    
    %% x (TIME) axis labelling
    deform_ratio_x = scr/size(sig_sliced,2); % screen is deformed - rescale
    xcross = crossings*deform_ratio_x;
    xstep = mean(diff(xcross)); % find step size
    if midbit
        xcross = xcross + xstep/2; % shift to the middle of the bit
        xcross = [xcross(1)-2*xstep,xcross(1)-xstep,xcross,xcross(end)+xstep]; % add extra step points
        xcross = xcross(xcross > 0 & xcross < screen_size(1)); % preserve if on screen
    end
    xcross_i = round(xcross/deform_ratio_x);
    xcross = round(xcross);
    xticks(xcross); % Assign x ticks
    
    
    if ~isempty(time)
        time_crossing = time(xcross_i); % assign times to ticks
        time_crossing = time_crossing-time_crossing(1);
    else
        time_crossing = 0:length(xcross); % otherwise assign bit number
    end
    
    [xprefix, xpoweroften] = get_prefix(time_crossing(end));
    ticklabels = round(time_crossing*10^(-xpoweroften),1); % tick text
    xticklabels(ticklabels);
    if isempty(time)
        xlabel("Bit Times");
    else
        xlabel(strcat("Time [",xprefix, "S]"));
    end
    
catch
    warning("Error Encountered in Locating Zero Crossings");
end

% image 0-388 y
% image 0-500 x

%% y (VOLT) axis labelling
maxabssig = max(abs([max(signal), min(signal)])); % use maximum abs value of signal as reference for our scale
[yprefix, ypoweroften] = get_prefix(maxabssig);
[sy, si] = sort(idxpair*2);
yticks(sy); % 0 is highest
sv = [v1, v2]*10^(-ypoweroften);
ticklabels = round(sv(si),1);
yticklabels(ticklabels);
ylabel(strcat("Voltage [",yprefix,"V]"));

axis on;


if plothist
    deform_ratio_x = scr/size(sig_sliced,2);
    axes(hXHist);
    h0 = histogram(zerocrossing_idxs,screen_size(1)/4);
    h0.EdgeColor = 'none';
    h0.FaceColor = [0,0,0];
    h0.FaceAlpha = 1;
    axis([0,screen_size(1)/deform_ratio_x, 0, inf]);
    
    axes(hYHist);
    h0 = histogram(signal,screen_size(2)/4);
    h0.EdgeColor = 'none';
    h0.FaceColor = [0,0,0];
    h0.FaceAlpha = 1;
    hYHist.View=[90,-90];
    axis([min(signal), max(signal), 0, inf]);
    
    
    
    pause(0.1);
    hEye.OuterPosition = [0,0,0.85,0.85];
    
    hXHist.OuterPosition = [0.041,0.8,0.805,0.2];
    hYHist.OuterPosition = [0.756,0.0418,0.2443,0.7945];
    
    hXHist.Position(1) = hEye.Position(1);
    hXHist.Position(3) = hEye.Position(3);
    
    hYHist.Position(2) = hEye.Position(2);
    hYHist.Position(4) = hEye.Position(4);
    
    
    axis(hEye, 'tight');
    axis(hXHist,'off');
    axis(hYHist,'off');
end
set(hEye,'FontSize',12,'FontWeight','bold');

    function [prefix, poweroften] = get_prefix(quant)
        pot = log10(quant);
        prefixes = ["p", "n", "u", "m","", "k", "M", "G", "T"];
        powers = -9:3:9;
        prefix_id = find(pot < powers,1);
        if prefix_id < 1 || prefix_id > 9
            warning("Could not assign unit prefix");
        end
        prefix = prefixes(prefix_id);
        poweroften = powers(prefix_id-1);
    end
end



