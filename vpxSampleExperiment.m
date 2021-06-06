function vpxSampleExperiment()
%
% this is an illustration of how to control the ViewPoint eyetracker form
% your stimulus script
%
% REQUIREMENTS:
% - Either you are running your stimuli on the same PC as the ViewPoint
%   Eyetracker or you have clinet link to the host PC running ViewPoint eyetracker
% - You have ViewPoint matlabe toolbox in your matlab path
%
% Based on the Arrington ViewPoint User's Manual
% Natalia Zaretskaya Feb 2010

% initialize eye tracker (change the paths according to your location of SDK)
PathOfDll='D:\Stereo_ViewPoint2.8.6.18\VPX_InterApp.dll' ;
PathOfHeader1='D:\Stereo_ViewPoint2.8.6.18\SDK\vptoolbox.h';
PathOfHeader2='D:\Stereo_ViewPoint2.8.6.18\SDK\vpx.h';
vpx_initialize(PathOfDll,PathOfHeader1, PathOfHeader2)

% set recording mode to speed for less quality but higher sampling rate
str = sprintf('%s', 'VideoMode Speed');
VPX_SendCommandString(str);

% pause writing into the data file
str = sprintf('%s', 'dataFile_Pause');
VPX_SendCommandString(str);

% configure synchronous string insertion
str = sprintf('%s', 'dataFile_asynchStringData No');
VPX_SendCommandString(str);

% open a file with name demo.txt
thestr_vpx = '"demo.txt"';
str = sprintf('dataFile_NewName %s', thestr_vpx );
vpx_SendCommandString (str);

% ----------------------------------------------------------------------- %
% do any preparatory phase of your experiment where you do not need to
% record the eye data, because the recording is paused
% ----------------------------------------------------------------------- %

% resume data recording and tell ViewPoint that the experiment starts
str = sprintf('%s', 'dataFile_Resume');
VPX_SendCommandString(str);
% insert marker idicating trial begin
str = sprintf('%s', 'dataFile_InsertString START');
VPX_SendCommandString(str);

for i = 1:nTrials
    
    % insert a marker indicating trial start (or anything else you like to 
    % tell the ViewPoint)
    str = sprintf('%s', 'dataFile_InsertString TRIAL', num2str(i));
    VPX_SendCommandString(str);
    WaitSecs(1)
end

% insert marker idicating experiment end
str = sprintf('%s', 'dataFile_InsertString END');
VPX_SendCommandString(str);
% close data file
str = sprintf('%s', 'dataFile_Close');
VPX_SendCommandString(str);
% unload the library
vpx_unload;

end