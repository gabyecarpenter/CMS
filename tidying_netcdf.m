%%Convert netcdf data to structure
trajlist = dir(fullfile(pwd,'/traj*'));

% Initialize the structure
bigstruct_0 = struct();

% Loop over the range of NetCDF files
for i = 1:length(trajlist)
    % Construct the filename
    filename = sprintf('traj_file_%02d.nc', i);  % Adjust format for two digits
    
    % Check if the file exists
    if isfile(filename)
        % Read the variables from the NetCDF file
        bigstruct_0(i).time = ncread(filename, 'time');
        bigstruct_0(i).location = ncread(filename, 'location');
        bigstruct_0(i).lon = ncread(filename, 'lon');
        bigstruct_0(i).lat = ncread(filename, 'lat');
        bigstruct_0(i).depth = ncread(filename, 'depth');
        bigstruct_0(i).distance = ncread(filename, 'distance');
        bigstruct_0(i).exitcode = ncread(filename, 'exitcode');
        bigstruct_0(i).releasedate = ncread(filename, 'releasedate');
    else
        warning('File %s does not exist.', filename);
    end
end

% bigstruct = bigstruct_0 
%Can create bigstruct immediatly if the number of release lines = number of
%traj files (not the case most of the time, so use the next loop)

%Depending on how many nodes/cores on the HPC, may need to adjust
%structure. This makes it so that each release line is a line in the
%structure. 

bigstruct = struct();
% bigstruct = struct('time', {}, 'location', {}, 'lon', {},'lat', {},'depth', {}, 'distance', {}, 'exitcode', {}, 'releasedate', {});

row_no = 0;
for i = 1:length(bigstruct_0)
    locs_in_file = unique(bigstruct_0(i).location);
    for j = 1:length(locs_in_file)
        row_no = row_no + 1
        r_index = bigstruct_0(i).location == locs_in_file(j);
        bigstruct(row_no).time = bigstruct_0(i).time;
        bigstruct(row_no).location = bigstruct_0(i).location(r_index);
        bigstruct(row_no).lon = bigstruct_0(i).lon(:,r_index);
        bigstruct(row_no).lat = bigstruct_0(i).lat(:,r_index);
        bigstruct(row_no).depth = bigstruct_0(i).depth(:,r_index);
        bigstruct(row_no).distance = bigstruct_0(i).distance(:,r_index);
        bigstruct(row_no).exitcode = bigstruct_0(i).exitcode(r_index);
        bigstruct(row_no).releasedate = bigstruct_0(i).releasedate(r_index);
    end
end



clear bigstruct_0

save traj_all.mat bigstruct -v7.3 %just saves the structure

%Quick figure to visualize trajectory files

figure()
hold on;
for i = 1:length(bigstruct) %count through each row of the structure and plot location
    plot([bigstruct(i).lon]-360, [bigstruct(i).lat]); 
end
axis equal
hold off;

% Let's bring in some land so we know what we're looking at.

S = shaperead('coastL1.shp');
for i = 1:length(S) %this for loop extracts the americas from a larger file
    Xloc = S(i).X;
    Yloc = S(i).Y;
    keepIndex = ~isnan(Xloc) & ~isnan(Yloc);
    Xloc = Xloc(keepIndex);
    Yloc = Yloc(keepIndex);
    fill(Xloc,Yloc,[.3 .3 .3]); hold on
end
for i = 1:length(bigstruct) %this plots our data
    plot([bigstruct(i).lon]-360, [bigstruct(i).lat]); 
end
axis equal
axis([-100,-80,14,32]) %this cuts it down