%% load data
clc
clear
dataset = {'Florence', 'Milan'};
for data_index = 1:length(dataset)
    dataname = dataset{data_index};
    load(fullfile('data',dataname,[dataname,'.mat']))
    Ohsi = imresize(im2double(hsi),[256,256]);
    [M,N,B] = size(Ohsi);
    if size(pan,1)~=M || size(pan,2)~=N
        Pan  = im2double(imresize(pan,[M,N]));
    else
        Pan  = im2double(pan);
    end

    clear pan hsi
    %% simulate noisy data
    for casenum = 1:5
        Nhsi = Ohsi;
        if casenum == 1
            NoiseG = true; gen_sig = 'sig = 10/255;';
            NoiseI = false;
            NoiseS = false;
        elseif casenum == 2
            NoiseG = true; gen_sig = 'sig = (rand(1)*25+5)/255;';
            NoiseI = false;
            NoiseS = false;
        elseif casenum == 3
            NoiseG = true; gen_sig = 'sig = (rand(1)*25+5)/255;';
            NoiseI = true; ipb = randperm(B); ipb = ipb(1:ceil(0.333*B)); gen_p = 'p = rand(1)*0.25+0.05;';
            NoiseS = false;
        elseif casenum == 4
            NoiseG = true; gen_sig = 'sig = (rand(1)*25+5)/255;';
            NoiseI = false;
            NoiseS = true; stb = randperm(B); stb = stb(1:ceil(0.333*B)); gen_s = 's = rand(1)*0.25+0.05;';
        elseif casenum == 5
            NoiseG = true; gen_sig = 'sig = (rand(1)*25+5)/255;';
            NoiseI = true; ipb = randperm(B); ipb = ipb(1:ceil(0.333*B)); gen_p = 'p = rand(1)*0.25+0.05;';
            NoiseS = true; stb = randperm(B); stb = stb(1:ceil(0.333*B)); gen_s = 's = rand(1)*0.25+0.05;';
        end

        for i=1:B
            % Gaussian noise
            if NoiseG
                eval(gen_sig)
                Nhsi(:,:,i) = Nhsi(:,:,i) + randn(M,N)*sig;
            end

            % Impulse noise
            if NoiseI && ismember(i,ipb)
                eval(gen_p)
                Nhsi(:,:,i) = imnoise(Nhsi(:,:,i),'salt & pepper', p);
            end

            % Deadlines
            if NoiseD && ismember(i,ddb)
                eval(gen_d)
                linenum   = ceil(N*d); % the number of deadlines
                linewidth = 1*ones(1,linenum); % the width of each deadline
                lineloc   = randi(N-max(linewidth),1,linenum); % the location of each deadline
                lineindex = [];
                for k=1:linenum
                    lineindex = [lineindex, lineloc(k):lineloc(k)+linewidth(k)];
                end
                Nhsi(:,lineindex,i) = 0;
            end

            % Stripes
            if NoiseS && ismember(i,stb)
                eval(gen_s)
                linenum = ceil(N*s);
                lineloc = ceil(N*rand(1,linenum));
                t = rand(1,length(lineloc))*0.5-0.25;
                Nhsi(:,lineloc,i) = bsxfun(@minus,Nhsi(:,lineloc,i),t);
            end
        end
        Nhsi(Nhsi>1) = 1; Nhsi(Nhsi<0)=0;
        Nhsi = uint8(Nhsi*255);

        %% set savepath
        savepath = fullfile('result',dataname,['case',num2str(casenum)]);
        mkdir(savepath)

        %% save simulations
        save([savepath,'\Noisy.mat'],'Nhsi','Ohsi')
        Nhsi = im2double(Nhsi);
        [m1, m2, m3, m4] = pwrctv_msqia(Ohsi,Nhsi);
        metrics = [m1, m2, m3, m4, 0];

        %% TDL
        % parameters
        noiselevel = std(reshape(Ohsi-Nhsi, [M*N,B]));
        vstbmtf_params.peak_value = 1;
        vstbmtf_params.nsigma = mean(noiselevel);
        % run algorithm
        disp('--------------- Run TDL --------------- ')
        tic
        output = TensorDL(Nhsi, vstbmtf_params);
        elapsed_time = toc;
        % save data
        output = uint8(255*output);
        save([savepath,'\TDL.mat'],'output')
        % metrics
        output = im2double(output);
        [m1, m2, m3, m4] = pwrctv_msqia(Ohsi,output);
        metrics = [metrics;m1, m2, m3, m4, elapsed_time];

        %% NGMeet
        % parameters
        noiselevel = std(reshape(Ohsi-Nhsi, [M*N,B]));
        Par   = ParSetH(255*mean(noiselevel),B);
        % run algorithm
        disp('--------------- Run NGMeet --------------- ')
        tic
        output = NGmeet_DeNoising( 255*Nhsi, 255*Ohsi, Par);  %NGmeet denoisng function
        elapsed_time = toc;
        % save data
        output = uint8(output);
        save([savepath,'\NGMeet.mat'],'output')
        % metrics
        output = im2double(output);
        [m1, m2, m3, m4] = pwrctv_msqia(Ohsi,output);
        metrics = [metrics;m1, m2, m3, m4, elapsed_time];

        %% BALMF
        % parameters
        r = 4;
        % run algorithm
        disp('--------------- Run BALMF --------------- ')
        tic
        output = BALMF(Nhsi, r);
        elapsed_time = toc;
        % save data
        output = uint8(255*output);
        save([savepath,'\BALMF.mat'],'output')
        % metrics
        output = im2double(output);
        [m1, m2, m3, m4] = pwrctv_msqia(Ohsi,output);
        metrics = [metrics;m1, m2, m3, m4, elapsed_time];

        %% TCTV
        % parameters
        opts = [];
        opts.rho = 1.25;
        opts.directions = [1,2,3];
        % run algorithm
        disp('--------------- Run TCTV --------------- ')
        tic
        output = TCTV_TRPCA(Nhsi, opts);
        elapsed_time = toc;
        % save data
        output = uint8(255*output);
        save([savepath,'\TCTV.mat'],'output')
        % metrics
        output = im2double(output);
        [m1, m2, m3, m4] = pwrctv_msqia(Ohsi,output);
        metrics = [metrics;m1, m2, m3, m4, elapsed_time];

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
        elapsed_time = toc;
        % save data
        output = uint8(255*output);
        save([savepath,'\LRTV.mat'],'output')
        % metrics
        output = im2double(output);
        [m1, m2, m3, m4] = pwrctv_msqia(Ohsi,output);
        metrics = [metrics;m1, m2, m3, m4, elapsed_time];

        %% LMHTV
        % parameters
        r = 4;
        tau = 0.0005;
        lambda = 10/sqrt(M*N);
        % run algorithm
        disp('--------------- Run LMHTV --------------- ')
        tic
        output = LMHTV(Nhsi, tau, lambda, r, [1,1,0]);
        elapsed_time = toc;
        % save data
        output = uint8(255*output);
        save([savepath,'\LMHTV.mat'],'output')
        % metrics
        output = im2double(output);
        [m1, m2, m3, m4] = pwrctv_msqia(Ohsi,output);
        metrics = [metrics;m1, m2, m3, m4, elapsed_time];

        %% LTHTV
        % parameters
        tau = 0.04;
        lambda = 1000/sqrt(M*N);
        ten_rank = [ceil(M*0.7),ceil(N*0.7),r];
        % run algorithm
        disp('--------------- Run LTHTV --------------- ')
        tic
        output = LTHTV(Nhsi, tau,lambda,ten_rank,[1,1,0]);
        elapsed_time = toc;
        % save data
        output = uint8(255*output);
        save([savepath,'\LTHTV.mat'],'output')
        % metrics
        output = im2double(output);
        [m1, m2, m3, m4] = pwrctv_msqia(Ohsi,output);
        metrics = [metrics;m1, m2, m3, m4, elapsed_time];

        %% CTV
        % parameters
        opts = [];
        opts.rho = 1.5;
        % run algorithm
        disp('--------------- Run CTV --------------- ')
        tic
        output = ctv_rpca(Nhsi, opts);
        elapsed_time = toc;
        % save data
        output = uint8(255*output);
        save([savepath,'\CTV.mat'],'output')
        % metrics
        output = im2double(output);
        [m1, m2, m3, m4] = pwrctv_msqia(Ohsi,output);
        metrics = [metrics;m1, m2, m3, m4, elapsed_time];

        %% RCTV
        % parameters
        beta = 100;
        lambda = 1;
        if casenum<=1
            tau = 0.4*[1,1];
        else
            tau = 0.7*[1,1];
        end
        r = 3;
        q = 5;
        % run algorithm
        disp('--------------- Run RCTV --------------- ')
        tic
        output = MoG_RBTV(Nhsi, beta, lambda, tau, r);
        elapsed_time = toc;
        % save data
        output = uint8(255*output);
        save([savepath,'\RCTV.mat'],'output')
        % metrics
        output = im2double(output);
        [m1, m2, m3, m4] = pwrctv_msqia(Ohsi,output);
        metrics = [metrics;m1, m2, m3, m4, elapsed_time];

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
        elapsed_time = toc;
        % save data
        output = uint8(255*output);
        save([savepath,'\WNLRATV.mat'],'output')
        % metrics
        output = im2double(output);
        [m1, m2, m3, m4] = pwrctv_msqia(Ohsi,output);
        metrics = [metrics;m1, m2, m3, m4, elapsed_time];

        %% PWRCTV
        % parameters
        beta = 100;
        lambda = 1;
        if casenum<=1
            tau = 0.4*[1,1];
            q = 10;
        else
            tau = 0.7*[1,1];
            q = 5;
        end
        r = 4;
        % run algorithm
        disp('--------------- Run PWRCTV --------------- ')
        tic
        output = PWRCTV(Nhsi, Pan, beta, lambda, tau, r, q);
        elapsed_time = toc;
        % save data
        output = uint8(255*output);
        save([savepath,'\PWRCTV.mat'],'output')
        % metrics
        output = im2double(output);
        [m1, m2, m3, m4] = pwrctv_msqia(Ohsi,output);
        metrics = [metrics;m1, m2, m3, m4, elapsed_time];

        save([savepath,'\metrics.mat'],'metrics')
    end
end

%% --------------------- Visualization ------------------------------------
casenum = 5;
dataset = {'Florence', 'Milan'};
for data_index = 1:length(dataset)
    dataname = dataset{data_index};
    path = ['result\',dataname,'\case',num2str(casenum)];
    % method = {'Noisy','Clean','LMHTV','LTHTV','LRTV','NGMeet','RCTV',...
    %     'WNLRATV','BALMF','CTV', 'RCILD', 'PWRCTV'};
    method = {'Noisy','Clean','LMHTV','LTHTV','LRTV','NGMeet','RCTV',...
        'WNLRATV','BALMF','CTV', 'PWRCTV'};

    data = {};
    for i = 1:length(method)
        if i~=2
            load(fullfile(path,method{i}))
        end
        if i==1
            data{i} = im2double(Nhsi);
            Ohsi = im2double(Ohsi);
        elseif i==2
            data{i} = Ohsi;
        else
            data{i} = im2double(output);
        end
    end

    if data_index==1
        band_index = [58,35,16];
    elseif data_index==2
        band_index = [56,41,1];
    end

    result = data;
    for i = 1:length(method)

        temp_image = result{i}(:,:,band_index);
        gt_image = result{2}(:,:,band_index);
        res_image = temp_image-gt_image;
        if i~=2
            temp_psnr = psnr(temp_image, gt_image);
            anotate_text = num2str(temp_psnr, '%2.2f');
        else
            anotate_text = 'PSNR';
        end
        anotate_text = [method{i},'/',anotate_text];
        temp_image = rsshow(temp_image);
        temp_image = insertText(temp_image,[1,1],anotate_text, ...
            'Font', 'Arial', ...
            'FontSize', 26, ...
            'BoxOpacity', 0.6, ... %
            'BoxColor', 'black', ...
            'TextColor', 'white');

        if data_index==1
            imamp_rect = [130,50,50,50];
        elseif data_index==2
            imamp_rect = [150,50,50,50];
        end

        imamp_linewidth = 2;
        imamp_scale = 2;
        imamp_location = 3;
        imamp_color = [1,0,0];
        imamp_alpha=1;
        colormap_ = jet(256);
        res_ratio = 5;
        temp_image = imamp(temp_image, imamp_rect, imamp_linewidth, ...
            imamp_scale, imamp_location, imamp_color, imamp_alpha);
        temp_image = imamp(temp_image, imamp_rect, imamp_linewidth, ...
            imamp_scale, 4, imamp_color, imamp_alpha, res_image, res_ratio, colormap_);

        montage_image{i} = temp_image;
    end

    h=montage(montage_image', ...
        'Size', [2,6],...
        'BackgroundColor', 'white', ...
        'BorderSize', 1);
    output=h.CData;
    
    mkdir('visualization')
    out_image_name = ['visualization',dataset{data_index},'_case',num2str(casenum),'.jpg'];
    imwrite(output, out_image_name, 'Quality', 100')
end