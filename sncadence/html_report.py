#multiplex cross_prod group_by "'all'"

import os
import os.path as op

r = get_input()

summary_files = glob_parent('combined_sim_history.npy', 'summary')
summary_dir = op.dirname(op.dirname(summary_files[0]))

merge_results_files = glob_parent('*/*.txt', 'merge_results')
merged_results_dir = op.dirname(op.dirname(merge_results_files[0]))

#movie_files = glob_parent('*.mp4', 'movies')
#print movie_files

output_dir = op.dirname(get_data_fn('cadence_io_dir'))

# cadence_io 
cmd = ['cadence_io.py', '-O', output_dir, '--summary_dir', summary_dir, '--merged_results_dir', merged_results_dir]
logged_subprocess(cmd)

# html_report 
cmd = ['html_report.py', '-O', output_dir, output_dir + os.sep + 'cadences.pkl']
logged_subprocess(cmd)





