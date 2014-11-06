addpath ~/Documents/Programs/Thor/trunk/

% booleans to turn on/off plots
watch_distro = 1;
watch_pPosit = 1;
watch_eFact  = 1;

exInfo = ncinfo('ex.nc');

%% enhancement factors
eFactInfo = ncinfo('ex.nc','enhancement_factor');
eFactLoad = ncread('ex.nc','enhancement_factor');

nZ     = size(eFactLoad, 1); nZmid = nZ - floor(nZ/2);
nY     = size(eFactLoad, 2); nYmid = nY - floor(nY/2);
nX     = size(eFactLoad, 3); nXmid = nX - floor(nX/2);
nTimes = size(eFactLoad, 4);

if (watch_eFact)
    for ii = 1:nTimes
        for jj = nYmid %1:2:nY
            contourf( squeeze(eFactLoad(:,jj,:,ii)))
            title(['Y = ',num2str(jj), ', t_i = ', num2str(ii)])
            caxis([0,10])
            pause(0.2)
        end
    end
end

%% Distros
distroInfo = ncinfo('ex.nc', 'distributions');
distroLoad = ncread('ex.nc', 'distributions');

nParams   = size(distroLoad, 1); % [x,y,z,D,rho,To,Do]
nCrystals = size(distroLoad, 2);
nDistros  = size(distroLoad, 3);
nTimes    = size(distroLoad, 4);

if (watch_distro)
    % watch distro # d_n
    d_n = nXmid + nZmid*(nZ-1);
    W = zeros(nCrystals,3,nTimes);
    for ii = 1:nTimes
        W(:,:,ii) = squeeze(distroLoad(1:3,:,d_n,ii))';
        Sims.Analyze.equalAreaPoint(W(:,:,ii))
        title(num2str(ii))
        pause(0.2)
    end
end

%% particle positions
pxInfo = ncinfo('ex.nc', 'p_x');
pxLoad = ncread('ex.nc', 'p_x');
pzInfo = ncinfo('ex.nc', 'p_z');
pzLoad = ncread('ex.nc', 'p_z');

if (watch_pPosit)
    for ii = 1:nTimes
        scatter( pxLoad(:,ii),pzLoad(:,ii) )
        title(num2str(ii))
        pause(0.2)
    end
end



