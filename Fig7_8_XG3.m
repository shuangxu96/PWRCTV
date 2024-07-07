addpath(genpath('compared_methods'))
addpath(genpath('data'))
addpath(genpath('utils'))
%% load data
clc
clear

dataset = {'XG3_Beijing_x1200_y1600','XG3_Yulin_x0920_y1624'};
for index = 1:length(dataset)
    dataname = dataset{index};
    load(fullfile('data',dataname,[dataname,'.mat']))
    D = 4095; % Dynamic Range = 2^12-1; XG3 12-bit image

    Nhsi = double(hsi)/D;
    [M,N,B] = size(Nhsi);
    if size(pan,1)~=M || size(pan,2)~=N
        Pan  = double(imresize(pan,[M,N]))/D;
    else
        Pan  = double(pan)/D;
    end

    clear pan hsi

    savepath = fullfile('result',dataname);
    mkdir(savepath)

    elapsed_time = [];

    %% -------------------------- Denoising ------------------------------
    %% save Nhsi&Pan
    Nhsi = uint16(D*Nhsi);
    Pan = uint16(D*Pan);
    save([savepath,'\Noisy.mat'],'Nhsi','Pan')
    Nhsi = double(Nhsi)/D;
    Pan  = double(Pan)/D;

    %% LTHTV
    % parameters
    r = 4;
    tau = 0.04;
    lambda = 1000/sqrt(M*N);
    ten_rank = [ceil(M*0.7),ceil(N*0.7),r];
    % run algorithm
    disp('--------------- Run LTHTV --------------- ')
    tic
    output = LTHTV(Nhsi, tau,lambda,ten_rank,[1,1,0]);
    t1= toc; elapsed_time = [elapsed_time;t1];
    % save data
    output = uint16(D*output);
    save([savepath,'\LTHTV.mat'],'output')

    %% LRTV
    % parameters
    if size(Nhsi,3) > 100
        tau = 0.015;
        lambda = 20/sqrt(M*N);
        rank = 10;
    else
        tau = 0.01;
        lambda = 10/sqrt(M*N);
        rank = 5;
    end
    % run algorithm
    disp('--------------- Run LRTV --------------- ')
    tic
    output = LRTV(Nhsi, tau, lambda, rank);
    t1= toc; elapsed_time = [elapsed_time;t1];
    % save data
    output = uint16(D*output);
    save([savepath,'\LRTV.mat'],'output')

    %% RCTV
    % parameters
    beta = 200;
    lambda = 1;
    tau = 0.7*[1,1];
    r = 4;
    q = 10;
    % run algorithm
    disp('--------------- Run RCTV --------------- ')
    tic
    output = MoG_RBTV(Nhsi, beta, lambda, tau, r);
    t1= toc; elapsed_time = [elapsed_time;t1];
    % save data
    output = uint16(D*output);
    save([savepath,'\RCTV.mat'],'output')

    %% BALMF
    % parameters
    r = 4;
    % run algorithm
    disp('--------------- Run BALMF --------------- ')
    tic
    output = BALMF(Nhsi, r);
    t1= toc; elapsed_time = [elapsed_time;t1];
    % save data
    output = uint16(D*output);
    save([savepath,'\BALMF.mat'],'output')

    %% CTV
    % parameters
    opts = [];
    opts.rho = 1.5;
    % run algorithm
    disp('--------------- Run CTV --------------- ')
    tic
    output = ctv_rpca(Nhsi, opts);
    t1= toc; elapsed_time = [elapsed_time;t1];
    % save data
    output = uint16(D*output);
    save([savepath,'\CTV.mat'],'output')

    %% PWRCTV
    % parameters
    beta = 100;
    lambda = 1;
    tau = 0.5*[1,1];
    r = 7;
    q = 1;
    % run algorithm
    disp('--------------- Run PWRCTV --------------- ')
    tic
    output = PWRCTV(Nhsi, Pan, beta, lambda, tau, r, q);
    t1= toc; elapsed_time = [elapsed_time;t1];
    % save data
    output = uint16(D*output);
    save([savepath,'\PWRCTV.mat'],'output')

    %% Non-Blind Denoising Algorithm
    Ohsi = double(output)/D;

    %% NGMeet
    % parameters
    noiselevel = std(reshape(Ohsi-Nhsi, [M*N,B]));
    Par   = ParSetH(255*mean(noiselevel),B);
    % run algorithm
    disp('--------------- Run NGMeet --------------- ')
    tic
    output = NGmeet_DeNoising( 255*Nhsi, 255*Ohsi, Par);  %NGmeet denoisng function
    t1= toc; elapsed_time = [elapsed_time;t1];
    % save data
    output = uint16(output/255*D);
    save([savepath,'\NGMeet.mat'],'output')

    %% WNLRATV
    % parameters
    noise     = reshape(Nhsi - Ohsi, M*N,B);
    Sigma_ratio  = std(noise(:));
    initial_rank  = 3;
    Rank = 6;
    ModelPar.alpha = 30;
    ModelPar.belta = 1;
    ModelPar.gamma = 0.08;
    param   = SetParam_NWT(Nhsi, Sigma_ratio);
    param.initial_rank = initial_rank;
    param.maxiter = 15;
    param.patnum        = 200;
    param.lambda        = 2e-1;
    [prior, model] = InitialPara( param,0,B);
    % run algorithm
    disp('--------------- Run WNLRATV --------------- ')
    tic
    output = WNLRATV2(Nhsi,Ohsi, Rank,ModelPar, param, model, prior);
    t1= toc; elapsed_time = [elapsed_time;t1];
    % save data
    output = uint16(D*output);
    save([savepath,'\WNLRATV.mat'],'output')

    %% -------------- Visualization ------------------
    path = 'result\XG3_Beijing_x1200_y1600';
    % method = {'Noisy','LTHTV','LRTV','NGMeet','RCTV','WNLRATV','BALMF',...
    %     'CTV','RCILD', 'PWRCTV'};
    method = {'Noisy','LTHTV','LRTV','NGMeet','RCTV','WNLRATV','BALMF',...
        'CTV', 'PWRCTV'};

    data = {};
    err = [];
    for i = 1:length(method)
        load(fullfile(path,method{i}))
        if i==1
            data{i} = double(Nhsi);
        else
            data{i} = double(output)/D;
        end
    end

    band_index = [150,149,20];

    result = data;
    for i = 1:length(method)

        temp_image = result{i}(:,:,band_index);
        anotate_text = method{i};
        temp_image = rsshow(temp_image, 0.05);
        temp_image = insertText(temp_image,[1,1],anotate_text, ...
            'Font', 'Arial', ...
            'FontSize', 46, ...
            'BoxOpacity', 0.6, ...
            'BoxColor', 'black', ...
            'TextColor', 'white');

        imamp_rect = [110,210,50,50];
        imamp_linewidth = 2;
        imamp_scale = 5;
        imamp_location = 4;
        imamp_color = [1,0,0];
        imamp_alpha=1;
        colormap_ = jet(256);
        res_ratio = 3;
        temp_image = imamp(temp_image, imamp_rect, imamp_linewidth, ...
            imamp_scale, imamp_location, imamp_color, imamp_alpha);

        montage_image{i} = temp_image;

    end

    h=montage(montage_image', ...
        'Size', [2,5],...
        'BackgroundColor', 'white', ...
        'BorderSize', 1);
    output=h.CData;

    mkdir('visualization')
    imwrite(output, ['visualization\',dataname,'.jpg'], 'Quality', 100)
end