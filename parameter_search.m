clear;

mix_wp_reg = 0:0.1:1;

for i=1:numel(mix_wp_reg)
    pars.dir_output                = ['/home/smajjk/Data/tmp-build-tpm-' num2str(mix_wp_reg(i))];
    pars.dat{1}.segment.mix_wp_reg = mix_wp_reg(i);
    build_template(pars);
end