%
function snalbedo = snowAlbedo_cryogrid(SNOWH_yday, SNOWH, AIRT, albei, melt_d)

if (AIRT <= 0)

    if (SNOWH_yday < SNOWH)

        snalbedo = 0.85;

    else

        snalbedo = albei - 0.008;
    
    end

else

    if (SNOWH_yday < SNOWH)

        snalbedo = 0.85;

    else

        snalbedo = 0.5 + exp(-0.24 * (melt_d)) * (0.85 - 0.5);
    
    end


end

if (snalbedo < 0.5)

    snalbedo = 0.5;

end

end

