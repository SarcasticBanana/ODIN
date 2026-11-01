function [t, x, y, T_duration] = load_event_data(fname, nx, ny, zeroIndexed)
    % LOAD_EVENT_DATA Ingests and formats raw asynchronous event telemetry.
    
    % Read the text file, checking if the first line is a header to skip
    fid = fopen(fname,'r'); firstLine = fgetl(fid); fclose(fid);
    hasHeader = isnan(str2double(strtok(firstLine,', ')));
    dat = readmatrix(fname,'NumHeaderLines',double(hasHeader));
    
    % Strip to the core 4 columns: [x, y, polarity, timestamp_in_microseconds]
    dat = dat(:,1:4);
    
    % Reorder to [time, x, y, polarity]. Doing this before unique() ensures 
    % that identical events are evaluated chronologically.
    dat = dat(:,[4 1 2 3]);            
    
    % Event camera hardware occasionally glitches and outputs the exact same 
    % event trigger twice. This permanently deletes those exact duplicates.
    dat = unique(dat,'rows');          
    
    % Convert hardware microsecond timestamps to seconds and set the very 
    % first recorded event as t = 0 to establish a relative timeline.
    t = (dat(:,1) - dat(1,1))/1e6;     
    
    % MATLAB requires 1-based indexing for matrices. If the sensor software 
    % outputs 0-based coordinates, this shifts them up by 1.
    x = dat(:,2) + zeroIndexed;        
    y = dat(:,3) + zeroIndexed;
    T_duration = t(end);
    
    % Sanity check to ensure no events are mapped outside the sensor resolution
    assert(min(x)>=1 && max(x)<=nx && min(y)>=1 && max(y)<=ny, ...
        'Pixel coordinates fall outside nx-by-ny. Check nx, ny and zeroIndexed.');
    fprintf('Loaded %d events, Duration: %.3f s\n', numel(t), T_duration);
end