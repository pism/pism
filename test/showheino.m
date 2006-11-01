function showheino(datdir,datprefix,runname);
% SHOWHEINO  Reads and displays time series data produced by
% ISMIP-HEINO run (i.e. a run like
%    "pisms -ismip H -datprefix MYGROUP -run ST ...").
% Example:
%    >> showheino('/home/user/pismdir/','MYGROUP','ST')
% will read files like /home/user/pismdir/MYGROUP_ST_ts_iv.dat
% If last year recorded in these files is 200k then will determine
% t1, t2, t3, t4.
%
% See ISMIP-HEINO documentation at
%    http://www.pik-potsdam.de/~calov/heino.html
% for meaning of titles, etc., and locations of points, and the
% meaning of "t1" etc.; see the "description of the model setup").
% [ELB 10/30/06]

filestart = [datdir datprefix '_' runname '_'];

figure
for j=1:2
    if j==1
        varname = 'iv';
        ylab = 'icevolume (10^6 km^3)';
    else % j==2
        varname = 'tba';
        ylab = 'temperate basal area (10^6 km^2)';
    end
    A = readHeino([filestart 'ts_' varname '.dat']);
    subplot(2,1,j), plot(A(1,:)/1000,A(2,:)), ylabel(ylab)
    if j==1, title('GLOBAL TIME SERIES'), end
    if j==2, xlabel('10^3 years'), end
end

figure
for j=1:3
    if j==1
        varname = 'ait';
        ylab = 'average ice thickness (km)';
    elseif j==2
        varname = 'ahbt';
        ylab = 'average homologous basal temp (K)';
    else % j==3
        varname = 'tba';
        ylab = 'temperate basal area (10^6 km^2)';
    end
    A = readHeino([filestart 'tss_' varname '.dat']);
    if j == 1
        A_ait = A;
    elseif j == 2
        A_ahbt = A;
    else % j==3
        A_tba = A;
    end
    subplot(3,1,j), plot(A(1,:)/1000,A(2,:)), ylabel(ylab)
    if j==1, title('SOFT SEDIMENT TIME SERIES'), end
    if j==3, xlabel('10^3 years'), end
end

for j=1:3
    if j==1
        varname = 'it';
        jtitle = 'ICE THICKNESS';
        ylab = 'thickness (km)';
    elseif j==2
        varname = 'hbt';
        jtitle = 'HOMOLOGOUS BASAL TEMPERATURE';
        ylab = 'homologous temperature (K)';
    else % j==3
        varname = 'bfh';
        jtitle = 'BASAL FRICTION HEATING';
        ylab = 'heating (W m^{-2})';
    end
    B = zeros(0,0);
    for k=1:7
        A = readHeino([filestart 'tsp' num2str(k) '_' varname '.dat']);
        if k == 1
            smallest = size(A,2);
        else
            smallest = min(smallest,size(A,2));
        end
        B(1:smallest,k) = A(2,1:smallest)';
    end
    figure
    plot(A(1,1:size(B,1))/1000,B)
    legend('P_1','P_2','P_3','P_4','P_5','P_6','P_7');
    title([jtitle ' TIME SERIES (AT POINTS P_1,...,P_7)'])
    xlabel('10^3 years'), ylabel(ylab)
end

disp('five figures created')

if max(A_ait(1,:)) >= 200000  % if run was completed
    A = A_ait(:, (A_ait(1,:) >= 150000)); % only consider last 50k yrs
    [yy,ii] = max(A(2,:));
    t1 = A(1,ii);
    [yy,ii] = min(A(2,:));
    t2 = A(1,ii)
else
    disp(['WARNING: run not completed?  last recorded (_tss_) year = ' ...
          num2str( max([max(A_ait(1,:)) max(A_ahbt(1,:)) max(A_tba(1,:))]) ) ...
          ' years'])
end
if max(A_ahbt(1,:)) >= 200000  % if run was completed
    A = A_ahbt(:, (A_ahbt(1,:) >= 150000)); % only consider last 50k yrs
    [yy,ii] = min(A(2,:));
    t3 = A(1,ii)
end
if max(A_tba(1,:)) >= 200000  % if run was completed
    A = A_tba(:, (A_tba(1,:) >= 150000)); % only consider last 50k yrs
    [yy,ii] = max(A(2,:));
    t4 = A(1,ii)
end

function A = readHeino(filename);
fileID = fopen(filename,'rt');
[A, count] = fscanf(fileID, '%f', [2, inf]);
disp(['  [ rows read from ' filename ': ' num2str(count/2) ' ]'])


