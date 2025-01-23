% Process ASCII output - SKIP IF YOUR OUTPUT IS IN NETCDF FORMAT
%The goal is to understand how each particle moves from its starting
%location. Leaving this code here in case there is a need to ouput in
%ASCII

trajlist = dir('traj_file*'); %create list of all trajectory file in the directory, pwd = what is my current location
    %can change to /output/traj if want to work in parent folder
tic %big chunk is going to run through all trajectory files and bring all of them into one structure
bigstruct = struct('loc',{},'lon',{},'lat',{},'dep',{},'flg',{},'tim',{});
for k = 1:length(trajlist)
    traj = trajlist(k).name;
    trajfile = load(strcat(trajlist(k).folder,'/',traj),'-ascii');

    locs = unique(trajfile(:,1)); % Would help to have the vector of poly IDs/release lines from the release file here
    outputrearrange = struct();
    for l = 1:length(locs)
        disp(k)
        disp(l)
        subtraj = trajfile(trajfile(:,1) == locs(l),:); %subset traj file by location
        outputrearrange(l).loc = locs(l);
        partcls = unique(subtraj(:,2));
        times = unique(subtraj(:,3));
        extstep = max(diff(times)); % This is likely the intended external time step
        baddies = [];
        baddies = find(rem(subtraj(:,3),extstep) > 0);
        % Check flags
        if (isempty(baddies))
            disp("No bad timesteps")
        elseif (any(subtraj(baddies,8) < 0))
            disp('There are flags at uneven timesteps')
            disp('Removing flagged lines from subtraj')
            disp('Remaking particles and unique times')
            subtraj(baddies,:) = [];
            partcls = unique(subtraj(:,2)); % Just being safe
            times = unique(subtraj(:,3));
            disp('Finishing removing baddies')
            if ~all(diff(times) == extstep)
                disp('Times still unequal STOP!')
            end
        else 
            disp('You have a bigger problem times are bad flags are good')
        end
        newlon = nan(length(times),length(partcls));
        newlat = nan(length(times),length(partcls));
        newdep = nan(length(times),length(partcls));
        newflg = nan(length(times),length(partcls));
        newtim = nan(length(times),length(partcls));
        for p = 1:length(partcls)
            subp = subtraj(subtraj(:,2) == partcls(p),:); %subset again by particle
            [~, I] = sort(subp(:,3));
            subp = subp(I,:);
            newlon(1:size(subp,1),p) = subp(:,4);
            newlat(1:size(subp,1),p) = subp(:,5);
            newdep(1:size(subp,1),p) = subp(:,6);
            newflg(1:size(subp,1),p) = subp(:,8); % DH added 1/31/24
            newtim(1:size(subp,1),p) = subp(:,3); % DH added 1/31/24
        end
        outputrearrange(l).lon = newlon;
        outputrearrange(l).lat = newlat;
        outputrearrange(l).dep = newdep;
        outputrearrange(l).flg = newflg;
        outputrearrange(l).tim = newtim;
    end
    bigstruct = [bigstruct,outputrearrange];
end
toc %bigstruct is the important one - have each of the locations (500) and a matrix for lat, long, and depth, (241 time stemps). 

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