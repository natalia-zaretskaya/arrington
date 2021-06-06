function varargout = vpxReadFile (varargin)
%
% function varargout = vpxReadFile (varargin)
% reads Viewpoint eyetracker text file
% USAGE:
% vpxReadFile()
% or
% vpxReadFile(fileName(s), filePath)
%
% INPUTS (optional):
% fileName - can be a character or sell array of strings with file name(s)
% filePath - path to where the files are stored
%
% if no input is provided the user will be prompted to select files with a
% ui window
%
% OUTPUT (optional):
% e - cell array of stricture arrays (one for each file) with fields:
% .time
%. x
% .y
% .pw (pupil width)
% .par (pupil aspect ratio)
% .valid (whether the sample is valid or not)
% .marker (user marker)
% .fname (name of file from which the data is taken)
%
% Natalia Zaretskaya 18.06.10
% extended to stereo mode 03.2012

if (nargin < 1)
    
    [fileName,filePath] = uigetfile('*.txt','choose file(s)', 'MultiSelect', 'on');
    
    
elseif nargin < 2
    
    error('please specify file name and path')
    
elseif nargin < 3
    
    fileName = varargin{1};
    filePath = varargin{2};
    
end

if ischar(fileName)
    fileName = {fileName};
end

e = cell(length(fileName),1);
for i = 1:length(fileName)
    
    fprintf('Reading: %s \n', fileName{i});
    
    
    % --- determine version --- %
    fid = fopen(fullfile (filePath, deblank(fileName{i})));
    [tmp] = textscan(fid,'%s %s %s %s %s %s  %s  %s  %s  %s %s %s %s %s %s %s %s %s %s %s %s %s','headerlines',0, 'returnOnError',0);
    fclose(fid);
    
    % determine what is in the colums; this infomration is written in row with code 5
    columnNameLine = find(strcmp(tmp{1}, '5'));
    columnNames = cellfun(@(x) x(columnNameLine(1),:), tmp);
    fprintf('Found following column names:\n');
    fprintf('  %s \n', columnNames{~cellfun('isempty', columnNames)})
    
    if ~isempty(columnNames{end}),
        stereo_flag = 1;
    else
        stereo_flag = 0;
    end
    
    
    
    % read the file starting from the first data line
    dataLine = find(strcmp(tmp{1}, '10'));  
    start_position = dataLine(1);
    
    fid = fopen(fullfile (filePath, deblank(fileName{i})));
    if strcmp(tmp{4}{1}, 'PC60')
        
        disp([ 'Product Version: '   tmp{4}{1}])
        
        fprintf('Using following fields: \n')
        fprintf(' %s \n', columnNames{[2 4 5 7 8 9 12]})
        
        % read gaze data
        % 1-fieldcode 2-time 3-dtime 4-xgaze 5-ygaze 6-region 7-pupilwidth 8-pupilaspect 9-quality 10-fixation 11-count 12-marker 13-string
        [tmp] = textscan(fid,'%*s %f %*s %f %f %*s  %f  %f  %f  %*s %*s %s %*s','headerlines',start_position-1, 'returnOnError',0);
        % 2-1 4-2 5-3 8-4 9-5 12-6
        fclose(fid);
        
        fprintf('data length: %u samples \n', length(tmp{2}));
        
        % raw values
        e{i}.time = [tmp{1}];
        
        e{i}.x = [tmp{2}];
        
        e{i}.y = [tmp{3}];
        
        e{i}.pw = [tmp{4}]; % pupil width ratio
        
        e{i}.par = [tmp{5}]; % pupil aspect ratio
        
        e{i}.valid = tmp{6}; % 0:pupil-glint both ok 1:pupil ok 2:pupil-glint pupil ok 3,4,5:pupil bad
        
        e{i}.marker = tmp{7}; % user marker
        
        e{i}.fname = fileName{i};
        %     out_var{i}.e = e;
        
    elseif strcmp(tmp{4}{1}, 'USB-220')
        
        disp([ 'Product Version: '   tmp{4}{1}])
        
        % read gaze data
        if stereo_flag
            fprintf('stereo tracking \n')
            
            % colums:
            % 1-fieldcode 2-'TotalTime' 3-'DeltaTime' 4-'X_Gaze' 5-'Y_Gaze' 6-'Region' 7-'PupilWidth' 8-'PupilHeight' 9-'Quality' 10-'Fixation'
            % 11-'TotalTime' 12-'DeltaTime' 13-'X_Gaze' 14-'Y_Gaze' 15-'Region' 16-'PupilWidth' 17-'PupilHeight' 18-'Quality' 19-'Fixation' 20-'Count' 21-'Marker' 22-'String'
            % need:
            
            fprintf('Using following fields: \n')
            fprintf(' %s \n', columnNames{[2 4 5 7 8 9 11 13 14 16 17 18 21]})
            
            fprintf('data length: %u samples \n', length(tmp{2}));
            
            %                          1   2  3   4  5  6    7  8    9   10  11  12 13 14 15  16 17  18 19  20  21  22%
            [tmp] = textscan(fid,'%*s %f %*s %f %f %*s  %f  %f  %f  %*s %f  %*s %f %f %*s %f %f  %f %*s %*s %s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s','headerlines',start_position, 'returnOnError',1);
            fclose(fid);
            
            % raw values
            fprintf('data length: %f \n', length(tmp{2}));

            e{i,1}.time = tmp{1};
            e{i,1}.x = tmp{2};
            e{i,1}.y = tmp{3};
            e{i,1}.pw = tmp{4};
            e{i,1}.par = tmp{4}./tmp{5};
            e{i,1}.valid = tmp{6};
            e{i,1}.marker = tmp{13}; % user marker
            e{i,1}.fname = fileName{i};
            
            e{i,2}.time = tmp{7};
            e{i,2}.x = tmp{8};
            e{i,2}.y = tmp{9};
            e{i,2}.pw = tmp{10};
            e{i,2}.par = tmp{10}./tmp{11};
            e{i,2}.valid = tmp{12};
            e{i,2}.marker = tmp{13}; % user marker
            e{i,2}.fname = fileName{i};
            
        else
            
            fprintf('Using following fields: \n')
            fprintf(' %s \n', columnNames{[2 4 5 7 8 9 12]})
            % colums:
            % 1-fieldcode 2-time 3-dtime 4-xgaze 5-ygaze 6-region 7-pupilwidth 8-pupilHEIGHT 9-quality 10-fixation 11-count 12-marker 13-string
            [tmp] = textscan(fid,'%*s %f %*s %f %f %*s  %f  %f  %f  %*s %*s %s %*s','headerlines',start_position-1, 'returnOnError',1);
            % 2-1 4-2 5-3 7-4 8-5 9-6 12-7
            fclose(fid);
            
            fprintf('data length: %u samples \n', length(tmp{2}));
            
            
            % raw values
            e{i}.time = [tmp{1}];
            
            e{i}.x = [tmp{2}];
            
            e{i}.y = [tmp{3}];
            
            e{i}.pw = tmp{4}; % par width            
                        
            e{i}.par = [tmp{5}./tmp{4}]; % par aspect ratio
            
            e{i}.valid = tmp{6}; % 0:par-glint both ok 1:par ok 2:par-glint par ok 3,4,5:par bad
            
            e{i}.marker = tmp{7}; % user marker
            
            e{i}.fname = fileName{i};
            
        end % if stereo)flag
        
        
    else error('failed to determine file version')
        
    end
    
end % for i

varargout{1} = e;

end % function


