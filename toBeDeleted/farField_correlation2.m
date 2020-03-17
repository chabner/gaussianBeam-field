function [C_ff,Cs_ff] = farField_correlation2(config,l_base,v_base,itersNum) 
%% GB code

[v_x,v_y] = ndgrid(v_base,v_base);
v_x = v_x(:).'; v_y = v_y(:).';

N_l = numel(l_base);
N_v = numel(v_base).^2;

C_ff = zeros(sqrt(N_v),sqrt(N_v),N_l);
Cs_ff = C_ff;

for lNum = 1:1:numel(l_base)
    l_x = l_base(lNum);
%     l = [-l_x,l_x;0,0;sqrt(1-l_x.^2),sqrt(1-l_x.^2)];
    l = [0,l_x;0,0;1,sqrt(1-l_x.^2)];
    
%     v1 = [(v_x-l_base(lNum));v_y;sqrt(1-(v_x-l_base(lNum)).^2-v_y.^2)];
    v1 = [(v_x);v_y;sqrt(1-(v_x).^2-v_y.^2)];
    v2 = v1 + l(:,2) - l(:,1);
    v2(3,:) = sqrt(1 - v2(1,:).^2 - v2(2,:).^2);
    v = [v1,v2];
    
    af_ang_vl = zeros(2*N_v, 2);
    if config.sct_type>1
        for j=1:2
            af_ang_vl(:,j)=evalampfunc_general(l(:,j)'*v,config.sct_type,config.ampfunc,config.dimNum);
        end
    else
        af_ang_vl=evalampfunc_general(0,config.sct_type,config.ampfunc,config.dimNum);
    end

    for iterNum = 1:1:itersNum
        [C,Cs] = MCfieldOnWave_correlation( ...
            af_ang_vl, ...
            [] , ...
            l,    ... dirl
            v,    ... dirv
            [1/config.scattgMFP,1/config.attMFP] ,  ... sigt
            1,                                      ... albedo
            config.box_min,                         ... box_min
            config.box_max,                         ... box_max
            l,                                      ... l
            v,                                      ... v
            true,                                   ... is_ff_l
            true,                                   ... is_ff_v
            config.iterationsRender,                ... maxItr
            config.wavelenght,                      ... lambda
            false,                                  ... doCBS
            config.sampleFlag,                      ... smpFlg
            config.sct_type,                        ... sct_type
            config.ampfunc,                         ... ampfunc
            config.ampfunc0                         ... ampfunc0
        );
    
%         u_ff = reshape(u_ff.',2,sqrt(N_v),2*sqrt(N_v));
%         u_ff = permute(u_ff,[3,2,1]);
%         
%         us_ff = reshape(us_ff.',2,sqrt(N_v),2*sqrt(N_v));
%         us_ff = permute(us_ff,[3,2,1]);
%         
%         C_ff(:,:,lNum) = C_ff(:,:,lNum) + ...
%             u_ff(1:sqrt(N_v),:,1) .* conj(u_ff((1+sqrt(N_v)):end,:,2));
%         Cs_ff(:,:,lNum) = Cs_ff(:,:,lNum) + ...
%             us_ff(1:sqrt(N_v),:,1) .* conj(us_ff((1+sqrt(N_v)):end,:,2));
%         
%         C_ff = permute(C_ff,[2,1,3]);
%         Cs_ff = permute(Cs_ff,[2,1,3]);
        C_ff(:,:,lNum) = C_ff(:,:,lNum) + reshape(C,sqrt(N_v),sqrt(N_v));
        Cs_ff(:,:,lNum) = Cs_ff(:,:,lNum) + reshape(Cs,sqrt(N_v),sqrt(N_v));
    end
end

end