import numpy as np
from itertools import combinations
import scipy.cluster.hierarchy as spc
import random
from math import exp


'''np.asarray([[1,0,0,1,0],
                       [0,1,1,0,1],
                       [1,0,0,0,0],
                       [0,1,1,0,0],
                       [0,0,0,1,0]])'''


class Algorithm:
    def __init__(self, data):
        self.target_val = -1
        self.matr = data.matr
        self.c_b = None
        self.m_num = data.m_num
        self.p_num = data.p_num

    def similar_measure(self):
        new_matr = np.ones((self.p_num, self.p_num))
        pairs = {}
        for i, j in combinations(range(0, self.p_num), 2):
            a_ij, b_ij, c_ij = 0, 0, 0
            for m_num in range(self.m_num):
                if self.matr[m_num][i] and self.matr[m_num][j]:
                    a_ij += 1
                elif self.matr[m_num][i] and not self.matr[m_num][j]:
                    b_ij += 1
                elif not self.matr[m_num][i] and self.matr[m_num][j]:
                    c_ij += 1
            coef = a_ij / (a_ij + b_ij + c_ij)
            pairs[(i, j)] = coef
            pairs[(j, i)] = coef
            new_matr[i][j] = coef
            new_matr[j][i] = coef
        pairs = sorted(pairs.items(), key=lambda kv: -kv[1])
        return new_matr, pairs

    def get_target_val(self, matr=None, cell_borders=None):
        if matr is None:
            matr = self.matr
        if cell_borders is None:
            cell_borders = self.c_b
        n_1 = np.count_nonzero(matr)
        n_1_in = 0
        n_0_in = 0
        for j in range(len(cell_borders)):
            if j == len(cell_borders) - 1:
                l_border = cell_borders[j]
                h_border = [None, None]
            else:
                l_border = cell_borders[j]
                h_border = cell_borders[j + 1]
            cell = matr[l_border[0]:h_border[0], l_border[1]:h_border[1]]
            n_1_in += np.count_nonzero(cell)
            n_0_in += cell.size - np.count_nonzero(cell)
        n_1_out = n_1 - n_1_in
        return (n_1 - n_1_out) / (n_1 + n_0_in)

    def get_parts_solution(self, similarity_matrix, num_clusters):
        pdist = spc.distance.pdist(similarity_matrix)
        linkage = spc.linkage(pdist, method='complete')
        cell_clusters = spc.fcluster(linkage, num_clusters, 'maxclust')
        idx_sort = np.argsort(cell_clusters)
        sorted_cell_clusters = cell_clusters[idx_sort]
        _, cell_borders = np.unique(sorted_cell_clusters, return_index=True)
        self.matr = self.matr[..., idx_sort]
        return cell_borders

    def get_machines_solution(self, matr=None, cell_borders=None, save=True):
        if matr is None:
            matr = self.matr
        if cell_borders is None:
            cell_borders = self.c_b
        if len(cell_borders.shape) == 2:
            cell_borders = cell_borders[..., 1]
        machines_matrix = np.zeros((self.m_num, len(cell_borders)))
        for j in range(len(cell_borders)):
            if j == len(cell_borders) - 1:
                l_border = cell_borders[j]
                h_border = None
            else:
                l_border = cell_borders[j]
                h_border = cell_borders[j + 1]
            machine = matr[..., l_border:h_border]
            voids = np.count_nonzero(matr, axis=1) - np.count_nonzero(machine, axis=1)
            exceptional = machine.shape[1] - np.count_nonzero(machine, axis=1)
            machines_matrix[..., j] = voids + exceptional
        cell_clusters = np.argmin(machines_matrix, axis=1)
        idx_sort = np.argsort(cell_clusters)
        sorted_cell_clusters = cell_clusters[idx_sort]
        _, cell_borders_machines = np.unique(sorted_cell_clusters, return_index=True)
        while cell_borders_machines.shape != cell_borders.shape:
            cell_borders_machines = np.append(cell_borders_machines, self.matr.shape[0])
        if save:
            self.matr = matr[idx_sort, ...]
            self.c_b = np.column_stack((cell_borders_machines, cell_borders))
        return matr[idx_sort, ...], np.column_stack((cell_borders_machines, cell_borders))

    def single_move(self):
        matr = self.matr.copy()
        best_obj_val = None
        best_matr = None
        best_borders = None
        overload = False
        for i in range(len(self.c_b)):
            if i == len(self.c_b) - 1:
                l_border = self.c_b[i][1]
                if l_border == matr.shape[1]:
                    h_border = l_border + 1
                    overload = True
                else:
                    h_border = matr.shape[1]
                    overload = True
            else:
                l_border = self.c_b[i][1]
                h_border = self.c_b[i + 1][1]
            source_iter_range = range(l_border, h_border)
            target_pos = [(k, border) for k, border in enumerate(self.c_b)
                          if (border != self.c_b[i]).all()]
            for source_pos in source_iter_range:
                source_part = matr[..., source_pos]
                for k, target_position in target_pos:
                    temp_matr = np.insert(matr, target_position[1], source_part, axis=1)
                    temp_borders = self.c_b.copy()
                    if source_pos < target_position[1]:
                        try:
                            temp_borders[i + 1][1] -= 1
                        except Exception:
                            pass
                        if i + 1 != k:
                            temp_borders[k][1] -= 1
                        temp_matr = np.delete(temp_matr, source_pos, axis=1)
                    else:
                        temp_borders[i][1] += 1
                        if i != k + 1:
                            try:
                                temp_borders[k + 1][1] += 1
                            except Exception:
                                pass
                        temp_matr = np.delete(temp_matr, source_pos + 1, axis=1)
                    obj_val = self.get_target_val(temp_matr, temp_borders)
                    if overload:
                        temp_borders[-1, 1] -= 1
                    if best_obj_val is None:
                        best_obj_val = obj_val
                        best_matr = temp_matr
                        best_borders = temp_borders
                    elif (obj_val - self.target_val) > (best_obj_val - self.target_val):
                        best_obj_val = obj_val
                        best_matr = temp_matr
                        best_borders = temp_borders
        return best_matr, best_borders

    def exchange_move(self):
        matr = self.matr.copy()
        best_obj_val = None
        best_matr = None
        for i in range(len(self.c_b)):
            if i == len(self.c_b) - 1:
                l_border = self.c_b[i][1]
                if l_border == matr.shape[1]:
                    h_border = l_border + 1
                else:
                    h_border = matr.shape[1]
            else:
                l_border = self.c_b[i][1]
                h_border = self.c_b[i + 1][1]
            src_iter_range = range(l_border, h_border)
            target_cells_indexes = [n for n, border in enumerate(self.c_b) if
                                    (border != self.c_b[i]).all()]
            for src_pos in src_iter_range:
                for n in target_cells_indexes:
                    if n == len(self.c_b) - 1:
                        l_target_border = self.c_b[n][1]
                        if l_target_border == matr.shape[1]:
                            h_target_border = l_target_border + 1
                        else:
                            h_target_border = matr.shape[1]
                    else:
                        l_target_border = self.c_b[n][1]
                        h_target_border = self.c_b[n + 1][1]
                    target_iter_range = range(l_target_border, h_target_border)
                    for target_position in target_iter_range:
                        temp_matr = matr.copy()
                        temp_matr[..., [src_pos, target_position]] = temp_matr[..., [target_position, src_pos]]
                        obj_val = self.get_target_val(temp_matr, self.c_b)
                        if best_obj_val is None:
                            best_obj_val = obj_val
                            best_matr = temp_matr
                        elif (obj_val - self.target_val) > (best_obj_val - self.target_val):
                            best_obj_val = obj_val
                            best_matr = temp_matr
        return best_matr, self.c_b

    def get_init_solution(self, cell_number):
        sim_m, sim_p = self.similar_measure()
        cell_borders = self.get_parts_solution(sim_m, cell_number)
        self.get_machines_solution(cell_borders=cell_borders)
        self.target_val = self.get_target_val()
        best_current_solution = self.matr.copy()
        best_current_cell_borders = self.c_b.copy()
        return best_current_solution, best_current_cell_borders

    def __call__(self, *args, **kwargs):
        cell_number = kwargs['num_clusters']
        opt_cell_num = cell_number
        T_f = kwargs['Tf']
        cool_rate = kwargs['cooling_rate']
        n_iter = kwargs['n_iter']
        period = kwargs['period']
        trap_stag_limit = kwargs['trap_stag_limit']
        max_cell_num = kwargs['max_cell_num']
        T = kwargs['T0']
        cnt = 0
        cnt_MC = 0
        cnt_trap = 0
        cnt_stag = 0
        best_global_solution, best_at_all_borders = self.get_init_solution(cell_number)
        best_current_solution = best_global_solution.copy()
        best_current_cell_borders = best_at_all_borders.copy()
        best_at_all_target = self.get_target_val(best_global_solution, best_at_all_borders)
        initial_solution = True
        generate_solution = True
        while True:
            if initial_solution:
                initial_solution = False
            elif generate_solution:
                best_current_solution, best_current_cell_borders = self.get_init_solution(cell_number)
                generate_solution = False
            while cnt_MC < n_iter and cnt_trap < n_iter / 2:
                new_solution, new_borders = self.single_move()
                if cnt % period == 0:
                    self.exchange_move()
                neighbor_solution, neigbour_borders = self.get_machines_solution(new_solution, new_borders, save=False)
                neigbour_target_value = self.get_target_val(neighbor_solution, neigbour_borders)
                if neigbour_target_value > self.get_target_val(best_current_solution, best_current_cell_borders):
                    best_current_solution = neighbor_solution.copy()
                    best_current_cell_borders = neigbour_borders.copy()
                    self.matr = neighbor_solution.copy()
                    self.c_b = neigbour_borders.copy()
                    cnt_stag = 0
                    cnt_MC += 1
                elif neigbour_target_value == self.target_val:
                    self.matr = neighbor_solution.copy()
                    self.c_b = neigbour_borders.copy()
                    cnt_stag += 1
                    cnt_MC += 1
                else:
                    delta = neigbour_target_value - self.target_val
                    prob = random.random()
                    if exp(-delta / T) > prob:
                        self.matr = neighbor_solution.copy()
                        self.c_b = neigbour_borders.copy()
                        cnt_trap = 0
                    else:
                        cnt_trap += 1
                    cnt_MC += 1
            if T < T_f or cnt_stag >= trap_stag_limit:
                if self.get_target_val(best_current_solution, best_current_cell_borders) > self.get_target_val(
                        best_global_solution, best_at_all_borders) \
                        and cell_number < max_cell_num:
                    best_global_solution = best_current_solution.copy()
                    best_at_all_borders = best_current_cell_borders.copy()
                    best_at_all_target = self.get_target_val(best_global_solution, best_at_all_borders)
                    opt_cell_num = cell_number
                    cell_number += 1
                    T = kwargs['T0']
                    generate_solution = True
                    cnt = 0
                    cnt_MC = 0
                    cnt_trap = 0
                    cnt_stag = 0
                else:
                    return T, best_at_all_target, best_at_all_borders, opt_cell_num
            else:
                T *= cool_rate
                cnt_MC = 0
                cnt_MC += 1
                cnt += 1
