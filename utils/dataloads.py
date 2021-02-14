import numpy as np


class CellsProductionData:
    def __init__(self, path_to_data):
        self.path = path_to_data
        self.matr = None
        self.m_num = None
        self.p_num = None

    def __call__(self, *args, **kwargs):
        with open(self.path, "r") as f:
            m, n = map(int, f.readline().split())
            self.m_num = m
            self.p_num = n
            self.matr = np.zeros((self.m_num, self.p_num), dtype=np.bool)
            for i in range(self.m_num):
                for j in list(map(int, f.readline().split()))[1:]:
                    self.matr[i][j - 1] = 1

