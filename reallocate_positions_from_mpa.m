function [out_comp] = reallocate_positions_from_mpa(out_comp)
global GLO
if ~GLO.modify_positions
    return
end
effectors = {'reaches','saccades'};
par = {'tar_pos','cue_pos'};
for out_idx  = 1:numel(out_comp)
    for e=1:numel(effectors)
        euclideans  = GLO.euclideans_reach;
        fix=0;
        for i_fn = 1:numel(par)
            for i_trial  = 1:numel(out_comp{out_idx}.(effectors{e}))
                old_tar =out_comp{out_idx}.(effectors{e})(i_trial).(par{i_fn});
                old_fix =out_comp{out_idx}.(effectors{e})(i_trial).fix_pos;
                if isnan(old_tar)
                    new_tar=old_tar;
                    targets_diff=0;
                else
                    [~, i_loc] = min(abs(euclideans-old_tar+old_fix));
                    new_tar= euclideans(i_loc);
                    targets_diff=euclideans(i_loc)-old_tar;
                end
                out_comp{out_idx}.(effectors{e})(i_trial).(par{i_fn}) = new_tar;
                out_comp{out_idx}.(effectors{e})(i_trial).fix_pos = fix;
                if strcmp(par{i_fn},'tar_pos')
                    old_mov = out_comp{out_idx}.(effectors{e})(i_trial).endpos;
                    out_comp{out_idx}.(effectors{e})(i_trial).endpos = old_mov+targets_diff;
                end
            end
        end
    end
end