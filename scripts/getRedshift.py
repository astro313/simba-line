import numpy as np

def get_z():

    import socket
    host = socket.gethostname()

    if 'flatironinstitute.org' or 'worker' in host:
        redshiftFile = '/mnt/ceph/users/daisyleung/simba/gizmo-extra/outputs_boxspace50.info'

        _, zs_table, snaps_table = np.loadtxt(redshiftFile, unpack=True)

    return zs_table, snaps_table


if __name__ == '__main__':
    zs_table, snaps_table = get_z()
    print(zip(zs_table, snaps_table))
