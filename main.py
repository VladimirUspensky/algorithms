import utils
from cfp import Algorithm
import os
import numpy as np


if __name__ == '__main__':

    kwargs = {'num_clusters': 2,
              'T0': 5,
              'Tf': 1e-3,
              'cooling_rate': 0.7,
              'n_iter': 4,
              'period': 6,
              'trap_stag_limit': 1,
              'max_cell_num': 20}

    for file in os.listdir('./datasets'):
        data = utils.CellsProductionData(os.path.join('./datasets', file))
        data()
        ann_sim_method = Algorithm(data)
        T, objective, borders, cell_num = ann_sim_method(**kwargs)
        with open(f'{file[:-4]}.sol', 'w') as file:
            m = np.arange(1, data.matr.shape[0] + 1)
            c = []
            for i in range(len(borders)):
                if i != len(borders) - 1:
                    l_border = borders[i][0]
                    h_border = borders[i + 1][0]
                else:
                    l_border = borders[i][0]
                    h_border = data.matr.shape[0]
                slice_size = range(l_border, h_border)
                for _ in slice_size:
                    c.append(i + 1)
            for m_idx, c_idx in zip(m, c):
                file.write(f'{c_idx} ')
            file.write('\n')
            parts = np.arange(1, data.matr.shape[1] + 1)
            clusters_for_parts = []
            for i in range(len(borders)):
                if i != len(borders) - 1:
                    l_border = borders[i][1]
                    h_border = borders[i + 1][1]
                else:
                    l_border = borders[i][1]
                    h_border = data.matr.shape[1]
                slice_size = range(l_border, h_border)
                for _ in slice_size:
                    clusters_for_parts.append(i + 1)
            for parts_id, parts_cluster_id in zip(parts, clusters_for_parts):
                file.write(f'{parts_cluster_id} ')
            file.write('\n')
