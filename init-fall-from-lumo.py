"""Generate fall_00000 from lumo_00000 for my phantom SPH APR test set 2ia-test-mini.

Run
    PYTHONPATH=/home/clmu/projects/20230201_MQ_RTinCEETest/src/sarracen python3 init-fall-from-lumo.py
in dump folder.
"""

import sarracen as sar

if __name__ == '__main__':

    FACTOR = 1./15.806274386

    sdf, sdf_sink = sar.read_phantom(f"lumo_00000")
    sdf.loc[:, ['vx', 'vy', 'vz']] *= FACTOR
    sdf_sink.loc[:, ['vx', 'vy', 'vz']] *= FACTOR
    sar.write_phantom(f"fall_00000.tmp", sdf, sdf_sink)
