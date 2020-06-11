function x = modified_zscore(x)

% median analog of the z-score
% https://www.ibm.com/support/knowledgecenter/SSEP7J_11.1.0/com.ibm.swg.ba.cognos.ug_ca_dshb.doc/modified_z.html

med = median(x);
mean_abs_dev = mad(x);
median_abs_dev = mad(x,1);

if median_abs_dev == 0 % this will happen if more then 50% of the data has identical values
    x = (x-med)/(1.253314*mean_abs_dev);
else
    x = (x-med)/(1.486*median_abs_dev);
end

end