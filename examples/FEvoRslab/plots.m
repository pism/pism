addpath ~/Documents/Programs/Thor/trunk/

% booleans to turn on/off plots
watch_distro = 0;
watch_pPosit = 0;
watch_eFact  = 0;

exInfo = ncinfo('ex.nc');

%% Distros
distroInfo = ncinfo('ex.nc', 'distributions');
distroLoad = ncread('ex.nc', 'distributions');

nParams   = size(distroLoad, 1); % [x,y,z,D,rho,To,Do]
nCrystals = size(distroLoad, 2);
nDistros  = size(distroLoad, 3);
nTimes    = size(distroLoad, 4);

if (watch_distro)
    % watch distro # d_n
    d_n = 3;
    W = zeros(nCrystals,3,nTimes);
    for ii = 1:nTimes
        W(:,:,ii) = squeeze(distroLoad(1:3,:,d_n,ii))';
        Sims.Analyze.equalAreaPoint(W(:,:,ii))
        title(num2str(ii))
        pause(0.5)
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
        pause(0.5)
    end
end


%% enhancement factors
eFactInfo = ncinfo('ex.nc','enhancement_factor');
eFactLoad = ncread('ex.nc','enhancement_factor');

nZ     = size(eFactLoad, 1);
nY     = size(eFactLoad, 2); nYmid = nY - floor(nY/2);
nX     = size(eFactLoad, 3);
nTimes = size(eFactLoad, 4);

if (watch_eFact)
    for ii = 1:nTimes
        for jj = 1:2:nY
            contourf( squeeze(eFactLoad(:,jj,:,ii)))
            title(['Y = ',num2str(jj), ', t_i = ', num2str(ii)])
            caxis([0,10])
            pause(0.2)
        end
    end
end
