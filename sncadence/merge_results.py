#multiplex cross_prod where (movies[0]==global_metrics[0]) group_by '(movies[0], movies[-1])'

import os
import os.path as op


cadence, sn = get_input()

movies = glob_parent('*.mp4')
plots = glob_parent('*.png', 'global_metrics')
data = glob_parent('*.npz')
data += glob_parent('*.npy')

output_dir = op.dirname(get_data_fn('toto'))
for m in movies:
    os.symlink(m, output_dir + os.sep + 'video.mp4')
for p in plots:
    os.symlink(p, output_dir + os.sep + op.basename(p))
for d in data:
    os.symlink(d, output_dir + os.sep + op.basename(d))

# need to add this so that this segment has at least one real original product 
# otherwise, glob_parent() on this segment returns nothing
target = get_data_fn('original_product.txt')
with open(target, 'w') as f:
    f.write('target')



