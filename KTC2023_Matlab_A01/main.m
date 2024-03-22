function main(inputFolder,outputFolder,categoryNbr)
path(path,'MiscCodes/')

%% General parameters
PAUSE_FOR_EACH_TARGET = 1;%set if you want to pause after reconstructing each target. if 0, do not plot anything
Nel = 32; %number of electrodes
z = 1e-6*ones(Nel,1); %contact impedances
load([inputFolder '/ref.mat']) %reference data from water chamber
load('Mesh_sparse.mat')
sigma0 = ones(length(Mesh.g),1); %linearization point

%% Undersampling of the data for the difficulty levels
vincl = true(31,76);
rmind = 1:2*(categoryNbr - 1);
for ii=1:76
    if any(Injref(rmind,ii) ~= 0)
        vincl(:,ii) = false;
    end
    vincl(rmind,ii) = false;
end
vincl = vincl(:);

%% Smoothness prior
corrlength = 0.115;
var_sigma = 0.05^2;
mean_sigma = sigma0;
smprior = SMprior(Mesh.g,corrlength,var_sigma,mean_sigma);

%% Forward solver
solver = EITFEM(Mesh, Mesh2, Injref, Mpat, vincl, []);

%% Noise model
noise_std1 = 0.05; %standard deviation for first noise component (relative to each voltage measurement)
noise_std2 = 0.01; %standard deviation for second noise component (relative to the largest voltage measurement)
solver.SetInvGamma(noise_std1,noise_std2,Uelref)


%% Loop for each file

files=dir(fullfile(inputFolder,'*.mat'));
%numbered .mat files are from targets with inclusions
%ref.mat is a reference data (the last file in alphabetical order)
for objectno = 1:length(files)-1
    load([inputFolder '/data' num2str(objectno) '.mat'])
    deltaU = Uel - Uelref; %difference data

    if PAUSE_FOR_EACH_TARGET
        Usim = solver.SolveForward(sigma0,z);
        figure, plot(Uelref(vincl)), hold on, plot(Usim)
        set(gcf,'Units','normalized','OuterPosition',[0.0 0.2 0.3 0.4])
    end

    % Step 1: linear difference reconstruction for the conductivity change
    J = solver.Jacobian(sigma0,z);
    lambda = 20;
    deltareco = (J'*solver.InvGamma_n(vincl,vincl)*J + lambda*smprior.L'*smprior.L)\J'*solver.InvGamma_n(vincl,vincl)*deltaU(vincl);

    if PAUSE_FOR_EACH_TARGET
        sgplot = sigmaplotter(Mesh,99+objectno,'parula');
        sgplot.basic2Dplot(deltareco,{['linear difference reconstruction ' num2str(objectno)]});
    end

    %interpolate the reconstruction into a pixel image
    deltareco_pixgrid = interpolateRecoToPixGrid(deltareco,Mesh, 256);

    if PAUSE_FOR_EACH_TARGET
        t = figure;
        subplot(1,2,1)
        imagesc(deltareco_pixgrid), colorbar, axis square, title('Step 1 output')
        subplot(1,2,2)
        ordenado = sort(deltareco_pixgrid(:));
        plot(ordenado),  axis square, title('Step 1 output values: Ascending order')
    end


    % Step 2: Soft thresholding
    T = max(max(abs(deltareco_pixgrid)))*0.14; %thresh
    deltareco_pixgrid = sign(deltareco_pixgrid).*max(abs(deltareco_pixgrid)-T, 0); %soft thresholding

    if PAUSE_FOR_EACH_TARGET
        t = figure;
        subplot(1,2,1)
        imagesc(deltareco_pixgrid), colorbar, axis square, title('Step 2 output')
        subplot(1,2,2)
        ordenado = sort(deltareco_pixgrid(:));
        plot(ordenado),  axis square, title('Step 2 output values: Ascending order')
    end

    % Step 3: Post-processing using a CNN
    net = importKerasNetwork('ultimate_cnn1.h5'); %import CNN


    %Data normalization: range [0,1], where 0.5 represents no conductivity variation
    C = 2.5*max(max(abs(deltareco_pixgrid))); %the scale for each case is different
    S_max = C;
    S_min = -C;
    delta_S = S_max - S_min;

    S_max_x = max(max(deltareco_pixgrid));
    S_min_x = min(min(deltareco_pixgrid));
    valor_max_proj = (S_max_x - S_min) /   delta_S;
    valor_min_proj = (S_min_x- S_min) /   delta_S;
    deltareco_pixgrid = normalize(deltareco_pixgrid(:), 'range', [valor_min_proj, valor_max_proj]);
    deltareco_pixgrid(isnan(deltareco_pixgrid))=0.5;
    deltareco_pixgrid = reshape(deltareco_pixgrid, [256 256]);

    %CNN inference
    deltareco_pixgrid = imresize(deltareco_pixgrid, [64 64]); %the input to the CNN is a 64x64 image
    output = predict(net, deltareco_pixgrid);
    deltareco_pixgrid = double(output);

    %resize back to 256x256 and visualize
    deltareco_pixgrid = imresize(deltareco_pixgrid, [256, 256]);

    if PAUSE_FOR_EACH_TARGET
        t = figure;
        subplot(1,2,1)
        imagesc(deltareco_pixgrid), colorbar, axis square, title('Step 3 output')
        subplot(1,2,2)
        ordenado = sort(deltareco_pixgrid(:));
        plot(ordenado),  axis square, title('Step 3 output values: Ascending order')
    end

    % Step 4: Pre-segmentation
    c1 = (max(max(deltareco_pixgrid)) - 0.50)./4;
    c2 = (min(min(deltareco_pixgrid)) - 0.50)./4;
    deltareco_pixgrid(deltareco_pixgrid>0.5 + c1) = 1; %conductive
    deltareco_pixgrid(deltareco_pixgrid < 0.48 + c2) = -0; %resistive
    if PAUSE_FOR_EACH_TARGET
        t = figure;
        subplot(1,2,1)
        imagesc(deltareco_pixgrid), colorbar, axis square, title('Step 4 output')
        subplot(1,2,2)
        ordenado = sort(deltareco_pixgrid(:));
        plot(ordenado),  axis square, title('Step 4 output values: Ascending order')
    end

    % Step 5: Morphological filtering
    se = strel('disk',12);
    deltareco_pixgrid = imopen(deltareco_pixgrid, se);

    if PAUSE_FOR_EACH_TARGET
        t = figure;
        subplot(1,2,1)
        imagesc(deltareco_pixgrid), colorbar, axis square, title('Step 5 output')
        subplot(1,2,2)
        ordenado = sort(deltareco_pixgrid(:));
        plot(ordenado),  axis square, title('Step 5 output values: Ascending order')
    end


    %% Otsu's method for segmentation and final reconstructed image

    %treshold the image histogram using Otsu's method
    [level,x] = Otsu2(deltareco_pixgrid(:),256,7);

    reconstruction = zeros(size(deltareco_pixgrid));

    ind1 = find(deltareco_pixgrid < x(level(1)));
    ind2 = find(deltareco_pixgrid >= x(level(1)) & deltareco_pixgrid <= x(level(2)));
    ind3 = find(deltareco_pixgrid > x(level(2)));

    %Check which class is the background (assumed to be the one with the
    %most pixels in it).
    [bgnum,bgclass] = max([length(ind1) length(ind2) length(ind3)]);
    switch bgclass
        case 1 %background is the class with lowest values - assign other two classes as conductive inclusions
            reconstruction([ind2]) = 2;
            reconstruction([ind3]) = 2;
        case 2 %background is the middle class - assign the lower class as resistive inclusions and the higher class as conductive
            reconstruction([ind1]) = 1;
            reconstruction([ind3]) = 2;
        case 3 %background is the class with highest values - assign the other two classes as resistive inclusions
            reconstruction([ind1]) = 1;
            reconstruction([ind2]) = 1;
    end

    if PAUSE_FOR_EACH_TARGET
        t = figure;
        subplot(1,2,1)
        imagesc(reconstruction), colorbar, clim([0 2])
        axis square, title('Reconstructed image values')
        title('Final')
        subplot(1,2,2)
        ordenado = sort(reconstruction(:));
        plot(ordenado),  axis square, title('Reconstructed image values: Ascending order')
    end

    save([outputFolder '/' num2str(objectno) '.mat'],'reconstruction')


    if PAUSE_FOR_EACH_TARGET
        pause();
        close all
    end

end

end