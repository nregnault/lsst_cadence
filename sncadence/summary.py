#multiplex cross_prod group_by "analyze[-1]"


import numpy as np

r = get_input()

print '*' * 72

to_merge = glob_parent('sim_history.npy', 'analyze')
print to_merge 

data = [np.load(nm) for nm in to_merge]

# merge the history ntuples into a big one 
d = np.hstack(data)

# extract just the final state from the NTuples
d_final_state = np.rec.fromrecords([dd[-1].tolist() for dd in data], 
                                   names = d.dtype.names)

# dump into the summary segment 
combined_sim_history_fn = get_data_fn('combined_sim_history.npy')
np.save(combined_sim_history_fn, d)

final_state_sim_history_fn = get_data_fn('final_state_sim_history.npy')
np.save(final_state_sim_history_fn, d_final_state)


# summary plots 
summary_plot_with_labels = get_data_fn('nsn_vs_z_labels.png')
cmd = ['cadence_plots.py', '-o', summary_plot_with_labels, '--title', '$N_{SN} vs. z$ [%s]' % r[0], final_state_sim_history_fn]
logged_subprocess(cmd)

summary_plot_no_labels = get_data_fn('nsn_vs_z_nolabels.png')
cmd = ['cadence_plots.py', '-o', summary_plot_no_labels, '--nolabels', '--title', '$N_{SN} vs. z$ [%s] % r[0]', final_state_sim_history_fn]
logged_subprocess(cmd)


seg_output = r
