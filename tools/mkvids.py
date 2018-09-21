#!/usr/bin/env python

import os
import os.path as op
import re
import glob
import subprocess
import argparse


pattern = re.compile('(\d+)_(.+)_(\d{5})_0_(\d{5})_0_(\d+)')

def main(root_dir='.', suffix='', snr_video=False):
    """
    """
    for dp, dn, fn in os.walk(root_dir):
        r = pattern.match(op.basename(dp))
        if r is None:
            continue
        _, cadence, _, _, nside = r.groups()
        cmd = ['ffmpeg', '-y', '-i', dp + os.sep + '%05d.png',
               '-pix_fmt', 'yuv420p',
               '-c:v', 'libx264', '-movflags', '+faststart',
               dp + os.sep + cadence + suffix + '.mp4']
        print ' '.join(cmd)
        subprocess.check_call(cmd)

        if snr_video:
            cmd = ['ffmpeg', '-y', '-i', dp + os.sep + 'snr%05d.png',
                   '-pix_fmt', 'yuv420p',
                   '-c:v', 'libx264', '-movflags', '+faststart',
                   dp + os.sep + cadence + suffix + '.snr.mp4']
            print ' '.join(cmd)
            subprocess.check_call(cmd)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='assemble the frames into one video')
    parser.add_argument('root_dir', nargs='+',
                        help='specify root directory (from where we search the frames)')
    parser.add_argument('--suffix',
                        default='', type=str,
                        help='suffix to append to the video name')
    parser.add_argument('--snr',
                        default=False,
                        action='store_true',
                        help='produce a video from the snr debug frame')
    args = parser.parse_args()
    print args
    
    for rd in args.root_dir:
        main(root_dir=rd, suffix=args.suffix, snr_video=args.snr)
