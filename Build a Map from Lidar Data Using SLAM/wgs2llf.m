%% =============== Function of wgs2llf =============== %%
function [LLF_NED_GT0] = wgs2llf(f_ref,initPOS)
    wgs84 = wgs84Ellipsoid;
    [N, E, D] = geodetic2ned(f_ref(:, 2), f_ref(:, 3), f_ref(:, 4), initPOS(1), initPOS(2), initPOS(3), wgs84);
    LLF_NED_GT0 = f_ref;
    LLF_NED_GT0(:, 2) = N;
    LLF_NED_GT0(:, 3) = E;
    LLF_NED_GT0(:, 4) = D;
end