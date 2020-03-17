clear

allExpirements = dir('config*.m');

rng('shuffle')
for iterNum = 1:1:1e0
    for expirementFile = 1:1:numel(allExpirements)
            clear config
            run([allExpirements(expirementFile).folder,filesep,allExpirements(expirementFile).name]);

            [config,gpuFunc] = preprocessConfig(config);
            
            C = 0;
            Cs1 = 0;
            xRep = 0;
            
            Nl = numel(config.focalPointsL.base);
            Nv = numel(config.focalPointsV.base);
            
            gridx = (-2e4:5e2:2e4)/2;
%             gridx = 0;
            [X,Y,Z] = ndgrid(gridx,gridx,-50);
            xVec = [X(:)';Y(:)';Z(:)'];
            
            smpXvec = 0 * xVec;
            
            corrMat_sum = zeros(size(X,1),size(X,2),Nl/2,'gpuArray');
            
            corrMat_up = corrMat_sum;
            corrMat_right = corrMat_sum;
            corrMat_axis = corrMat_sum;
            
            deltaPix = 30;
            deltaAxis = round(deltaPix/sqrt(2));
            centerPoint = (Nv+1)/2;
            
            Cs = 0;
            
            tic
            for corrIter = 1:1:size(xVec,2)
                [~,us,smpX] = run_renderingDebug(config,gpuFunc,xVec(:,corrIter));
%                 u = reshape(u,[Nl,Nv,Nv]); u = permute(u,[2,3,1]);
                us = reshape(us,[Nl,Nv,Nv]); us = permute(us,[2,3,1]);
                
%                 C = C + u(:,:,1) .* conj(u);
                Cs1 = us(:,:,(1:Nl/2)) .* conj(us(:,:,((1+Nl/2):end)));
                Cs = Cs + us(:,:,(1:Nl/2)) .* conj(us(:,:,((1+Nl/2):end)));
                smpXvec(:,corrIter) = smpX;
                
                [I,J] = ind2sub(size(X),corrIter);
                corrMat_sum(I,J,:) = squeeze(sum(sum(abs(Cs1).^2,1),2));
%                 corrMat_up(I,J,:) = squeeze(abs(Cs(centerPoint+deltaPix,centerPoint,:)).^2);
%                 corrMat_right(I,J,:) = squeeze(abs(Cs(centerPoint,centerPoint+deltaPix,:)).^2);
%                 corrMat_axis(I,J,:) = squeeze(abs(Cs(centerPoint+deltaAxis,centerPoint+deltaAxis,:)).^2);
                corrMat_up(I,J,:) = squeeze(Cs1(centerPoint+deltaPix,centerPoint,:));
                corrMat_right(I,J,:) = squeeze(Cs1(centerPoint,centerPoint+deltaPix,:));
                corrMat_axis(I,J,:) = squeeze(Cs1(centerPoint+deltaAxis,centerPoint+deltaAxis,:));
            end
            t = toc
            numOfIterations = 1e4;
            
            corrMat_sum = gather(corrMat_sum);
            corrMat_up = gather(corrMat_up);
            corrMat_right = gather(corrMat_right);
            corrMat_axis = gather(corrMat_axis);

            T = datetime('now','Format','ddMM_HHmmss_SSS');
            save(['res',filesep,config.projectName,'_',num2str(iterNum),'_',char(T),'.mat'], ...
                'C','Cs','config','t','xRep','numOfIterations');
    end
end
