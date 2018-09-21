#!/usr/bin/env python

import os
import os.path as op
import sys
import glob
import re
import logging
logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',
                    level=logging.INFO)
import Image



pattern = re.compile('(\d+)_(.+)_(\d{5})_0_(\d{5})_0_(\d+)')

def crop_last_frame(path, cadence_name, nside):
    filenames = glob.glob(path + os.sep + '0*.png')
    filenames.sort()
    fn = filenames[-1]
    im = Image.open(fn)
    nx, ny = im.size
    cropped = im.crop((0,0,nx,0.52*ny))
    dst = path + os.sep + cadence_name + '_%s_maps.png' % nside
    logging.info(' -> %s' % op.basename(dst))
    cropped.save(dst)
    
    
def main(root_dir='.', suffix=''):
    for dp, dn, fn in os.walk(root_dir):
        r = pattern.match(op.basename(dp))
        if r is None:
            continue
        _, cadence, _, _, nside = r.groups()
        logging.info(dp)
        crop_last_frame(dp, cadence, nside)
        
if __name__ == "__main__":
    main()
