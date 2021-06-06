function vpxCalibrate(npoints, varargin)
% %--------------------------------------------------------------------------
% vpxCalibrate(npoints, varargin)
% presents calibration stimulus on the screen. It performs calibration by
% communicating with the Arrington ViewPoint Eyetracker.
% If the number of stimulus points is not mentioned the default value is 12.
%
% REQUIREMENTS:
% - Either you are running your stimuli on the same PC as the ViewPoint
%   Eyetracker or you have clinet link to the host PC running ViewPoint eyetracker
% - You have ViewPoint matlabe toolbox in your matlab path
%
% USAGE:
%   vpxCalibrate(npoints)
%   or
%   vpxCalibrate(npoints, 'Parameter', 'Value')
%
% INPUT:
% npoints: [default = 12] number of calibration points
%
% OPTIONAL PARAMETERS:
% setupId: [default = 0] scalar matching the setup definiiton
%   (including path to your Eyetracker SDK). See examples in the script.
%   Default is dummy mode without eyetracker linked.
%
% backGroundColor [default = 127 127 127] rgb of the background. It is
%   advised to match the luminance of the calibration background to the
%   luminance in your experiment
%
% doGammaCorrection = [default = 1] 1/0 to either do gammag correction
%   (make monitor linear) or not. If you are using gamma correction in your
%   simulus script it is useful to do so during calibraiton as well.
%   Otherwise the luminance will not be fully matched.
%
% based on the original ViewPoint toolbox script
% extended Natalia Zaretskaya Feb 2010
% updated to stereo mode
%--------------------------------------------------------------------------

% set defaults
setupId = 0;
backgroundColor = [127 127 127];
doGammaCorrection = 1;
stereoMode = 0;

if ~isempty(varargin)
    for i = 1:length(varargin)
        if strcmp(varargin{i}, 'setupId')
            setupId = varargin{i+1};
        elseif strcmp(varargin{i}, 'backgroundColor')
            backgroundColor = varargin{i+1};
        elseif strcmp(varargin{i}, 'doGammaCorrection')
            doGammaCorrection = varargin{i+1};
        elseif strcmp(varargin{i}, 'stereoMode')
            stereoMode = varargin{i+1};
        end
    end
end

% other settings
textColor = [255 0 0 ];
bigRectangleColor = [0 255 0];
smallRectangleColor = [255 0 0];
gammaLUTfile = '';

% configure eye tracker
if setupId == 1 % mri
    gammaLUTfile = 'MPI_3T_GammaTable.txt';
    pathofdll='D:\Arrington Eye Tracker\ViewPoint (2.8.4 beta 33d)\VPX_InterApp.dll' ;
    pathofheader1='D:\Arrington Eye Tracker\ViewPoint (2.8.4 beta 33d)\SDK\vptoolbox.h';
    pathofheader2='D:\Arrington Eye Tracker\ViewPoint (2.8.4 beta 33d)\SDK\vpx.h';
    % vpx_initialize(PathOfDll,PathOfHe ader1, PathOfHeader2)
    vpx_initialize(pathofdll,pathofheader1, pathofheader2)
    str = sprintf('%s', 'VideoMode Precision'); % 30 Hz
%     str = sprintf('%s', 'VideoMode Speed'); % 60 Hz
    VPX_SendCommandString(str);
    
elseif setupId == 2 % tms
    gammaLUTfile = 'gammaLUT_DELL10031148.txt';
    % initialize
    PathOfDll='D:\ViewPoint USB 220\Interfaces\Windows_XP_Vista\ViewPointClient Ethernet Interface\VPX_InterApp.dll' ;
    PathOfHeader1='D:\ViewPoint USB 220\SDK\vptoolbox.h';
    PathOfHeader2='D:\ViewPoint USB 220\SDK\vpx.h';
    vpx_initialize(PathOfDll,PathOfHeader1, PathOfHeader2)
    
elseif setupId == 3 % tms new viewpoint version
    gammaLUTfile = 'gammaLUT_DELL10031148 .txt';
    % initialize
    PathOfDll='D:\ViewPoint2.8.6.18\VPX_InterApp.dll' ;
    PathOfHeader1='D:\ViewPoint2.8.6.18\SDK\vptoolbox.h';
    PathOfHeader2='D:\ViewPoint2.8.6.18\SDK\vpx.h';
    vpx_initialize(PathOfDll,PathOfHeader1, PathOfHeader2)
    
    
elseif setupId == 4 % tms new viewpoint version STEREO
    gammaLUTfile = 'gammaLUT_samsung_syncmaster2443.txt';
    % initialize
    PathOfDll='D:\Stereo_ViewPoint2.8.6.18\VPX_InterApp.dll' ;
    PathOfHeader1='D:\Stereo_ViewPoint2.8.6.18\SDK\vptoolbox.h';
    PathOfHeader2='D:\Stereo_ViewPoint2.8.6.18\SDK\vpx.h';
    vpx_initialize(PathOfDll,PathOfHeader1, PathOfHeader2)
    
else
    fprintf('No eyetracker setup specified. Not initializing the library. \n')
end

if setupId > 0
    % set recording mode to precision (not sure if this is available in newer versions)
    str = sprintf('%s', 'VideoMode Precision');
    VPX_SendCommandString(str);
end

% Psychtoolbox settings
warning('off','MATLAB:dispatcher:InexactCaseMatch')
ListenChar(2)
HideCursor
KbName('UnifyKeyNames'); % same key indices for Windows and Mac
Screen('Preference', 'SkipSyncTests', 1); % skip sinchrinization tests
screenId = max(Screen('Screens'));
[scrWidth, scrHeight]=Screen('WindowSize', screenId);
if setupId <=0 % open a small window for testing
    scrWidth = scrWidth/3;
    scrHeight = scrHeight/3;
end


% ganmma correction
if doGammaCorrection
    oldGamma = Screen('ReadNormalizedGammaTable', screenId);
    if isempty(gammaLUTfile)
        newGamma = makeArtificialGamma(2.2);
    else
        newGamma = load(gammaLUTfile);
    end
    Screen('LoadNormalizedGammaTable', screenId, newGamma);
end

PsychImaging('PrepareConfiguration');
% [w,scrRect]=Screen('OpenWindow',screenId,0, [0  0 scrWidth scrHeight]);
[w, scrRect] = PsychImaging('OpenWindow', screenId, backgroundColor, [], [], [], stereoMode, []);
% [ptb.w.id, ptb.w.rect] = Screen('OpenWindow', slaveScreen, BlackIndex(slaveScreen), [], [], [], ptb.stereo.mode, ptb.numMultiSamples);

escape = KbName('escape');


if(nargin==0)
    npoints=12;
end

bigRectangle = scrWidth/25; % big rectangle size
smallRectangle = scrWidth/240; % small rectangle size


if setupId
    % calibration settings
    VPX_SendCommandString ('calibrationRealRect 0.2 0.2 0.8 0.8');
    strfinal=strcat(['calibration_points',blanks(1),num2str(npoints)]);
    vpx_SendCommandString(strfinal);
    vpx_SendCommandString('calibration_snapMode ON');
    vpx_SendCommandString('calibration_autoIncrement ON');
end


% draw the initial screen
Screen('FillRect', w, backgroundColor);
for i = 1:npoints
    if setupId
        [x,y]=vpx_GetCalibrationStimulusPoint(i);
    else
        x= 0.5; y=0.5;
    end
    
    
    %     currentRectangle = bigRectangle;
    rect2=[(x*scrRect(3))-smallRectangle,(y*scrRect(4))-smallRectangle,(x*scrRect(3))+smallRectangle,(y*scrRect(4))+smallRectangle];
    
    % buffer 1
    Screen('SelectStereoDrawBuffer', w, 0);
%     Screen('FillRect', w, smallRectangleColor, [rect2]');
    DrawFormattedText(w,...
        'eyetracker calibration: \n Follow receeding green squares with your gaze \n Try not to blink when the square gets very small',...
        'center', 'center', textColor);
    DrawFormattedText(w, num2str(i), x, y, textColor);
    
    % buffer 2
    Screen('SelectStereoDrawBuffer', w, 1);
%     Screen('FillRect', w, smallRectangleColor, [rect2]');
    DrawFormattedText(w,...
        'eyetracker calibration: \n Follow receeding green squares with your gaze \n Try not to blink when the square gets very small',...
        'center', 'center', textColor);
    DrawFormattedText(w, num2str(i), x, y, textColor);
    
    Screen('DrawingFinished', w);
end


Screen('Flip', w);
WaitSecs(1);
KbWait;
breakFlag = 0;
% start calibration
for i = randperm(npoints) % 1 : npoints
    
    if breakFlag
        break
    end
    
    if setupId
        [x,y]=vpx_GetCalibrationStimulusPoint(i);
    else x = 0.5; y = 0.5;
    end
    
    currentRectangle = bigRectangle;
    
    while currentRectangle >= smallRectangle
        
        [~, ~, keyCode, ~] = KbCheck();
        if keyCode(escape)
            breakFlag = 1;
        end
        
        rect=[(x*scrRect(3))-currentRectangle,(y*scrRect(4))-currentRectangle,...
            (x*scrRect(3))+currentRectangle,(y*scrRect(4))+currentRectangle];
        
        rect2=[(x*scrRect(3))-3,(y*scrRect(4))-3,...
            (x*scrRect(3))+3,(y*scrRect(4))+3];
        
        
        % buffer 1
        Screen('SelectStereoDrawBuffer', w, 0);
        Screen('FillRect', w, backgroundColor);
        Screen('FillRect', w,[bigRectangleColor; smallRectangleColor]', [rect; rect2]');
        currentRectangle = currentRectangle - 1;
        
        
        % buffer 2
        Screen('SelectStereoDrawBuffer', w, 1);
        Screen('FillRect', w, backgroundColor);
        Screen('FillRect', w,[bigRectangleColor; smallRectangleColor]', [rect; rect2]');
        currentRectangle = currentRectangle - 1;
        
        Screen('DrawingFinished', w);
        Screen('Flip', w);
    end
    
    if setupId
        strfinal1=strcat(['calibration_snap',blanks(1),num2str(i)]);
        vpx_SendCommandString(strfinal1);
    end
    
    WaitSecs(1);
    
end



ListenChar(0)
Screen('CloseAll');

if doGammaCorrection
    Screen('LoadNormalizedGammaTable', screenId, oldGamma);
end

if setupId
    vpx_unload
end
ShowCursor
end



function newGamma = makeArtificialGamma(desiredGammaExponent)
% Create artificial lookup table for the monitor if luminance was not
% measured. This is not ideal, but will at least approximate the gamma of
% most monitors.

% create arifical luminance measurment of each gun
they=linspace(0,1,256).^desiredGammaExponent;
they = round(255*repmat(they',1,3));
rgb_io = they;

% normalize the measured luminance
for i=1:min(size(rgb_io)),
    rgb_io(:,i)=rgb_io(:,i)-min(rgb_io(:,i)); % remove min from each
    rgb_io(:,i)=rgb_io(:,i)/max(rgb_io(:,i)); % normalize max : now rgb_io has range : [0,1].
end;

% Cheat: Make rgb_io monotonically increasing:
%		at start and/or end gamma-function may saturate and random fluctuations may make it non-monotonic.
for i=1:min(size(rgb_io)),
    for k=1:size(rgb_io,1)-1,
        if rgb_io(k,i)>=rgb_io(k+1,i), % if preceding value is >= to next value
            rgb_io(k+1,i)=rgb_io(k,i)+0.0001; % increase by 1/10000.
        end;
    end;
end;

thex=linspace(0, 255, 256); % range of values in LUT
they=linspace(0,1,256).^1; % [0,1] (desired, (not necessarily linear) rgb output function).

for i=1:min(size(rgb_io)),
    LUT(:,i)=round(interp1(rgb_io(:,i)',thex, they,'spline'))';
end;

newGamma = LUT/256; % lookup table to make monitor linear, [0 1]

end


