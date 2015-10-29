# -*- coding: utf-8 -*-
__author__ = 'bod'

import copy
from simplexengine import *
from numpy import matrix


class SimplexTable:
    """
    Implements simplex method table
    """

    def __init__(self, func, **kwargs):
        """
        :param extrem_functions: detect what extremum type is max or min. Should contain select function
        :param kwargs: table elements, e.g. 'A' - vectors, 'B' - start basis, 'C' -
        coefficients beside x
        :return: new instance
        """
        self.__table = dict()
        self.__function = func
        for k, v in kwargs.items():
            self.__table[k] = v
        # initialize deltas as zero
        self.__table['D'] = [0] * len(self.__table['A'])
        self.__evaluate_deltas()

    def __evaluate_deltas(self):
        for i in xrange(len(self.__table['D'])):
            if i == 0:
                self.__evaluate_delta_null()
            else:
                self.__evaluate_delta(i)

    def __evaluate_delta_null(self):
        self.__table['D'][0] = 0
        for index, basis in enumerate(self.__table['B']):
            self.__table['D'][0] += self.__table['C'][basis] * self.__table['A'][0][index]

    def __evaluate_delta(self, index):
        self.__table['D'][index] = 0
        for i, basis in enumerate(self.__table['B']):
            self.__table['D'][index] += self.__table['C'][basis] * self.__table['A'][index][i]

        self.__table['D'][index] -= self.__table['C'][index]

    def __reevaluate(self):
        self.__table_back = copy.deepcopy(self.__table['A'])
        for i in xrange(len(self.__table['A'][0])):
            if i == self.__forward[0]:
                self.__reevaluate_forward_row()
            else:
                self.__reevaluate_row(i)

    def __reevaluate_forward_row(self):

        for i in xrange(len(self.__table['A'])):
            self.__table['A'][i][self.__forward[0]] /= float(self.__forward_value)

    def __reevaluate_row(self, index):
        a_ik = self.__table['A'][self.__forward[1]][index]
        for j in xrange(len(self.__table['A'])):
            a_rj = self.__table_back[j][self.__forward[0]]
            self.__table['A'][j][index] -= float((a_rj * a_ik)) / self.__forward_value

    def __find_forward_element(self):
        forward_row = self.__function.get_forward_row(self.__table['D'])
        min_element, min_index, i = -1, -1, 0
        for a0, ak in zip(self.__table['A'][0], self.__table['A'][forward_row]):
            if ak > 0:
                case = float(a0) / ak
                if min_element == -1:
                    min_index = i
                    min_element = case
                if self.__function.compare(case, min_element):
                    min_index = i
                    min_element = case
            i += 1
        self.__forward = [min_index, forward_row]
        self.__forward_value = self.__table['A'][forward_row][min_index]

    def next_iteration(self):
        """
        :return: code False if extremum element doesn't find, otherwise returns True
        :raises: Exception if a function is unbounded
        """
        self.__find_forward_element()
        if self.__function.constraint_on_unbounded(self.__table['A'][self.__forward[1]]):
            raise UnboundedException()
        if self.__function.constraint_on_find(self.__table['D']):
            return False
        self.__table['B'][self.__forward[0]] = self.__forward[1]
        self.__reevaluate()
        self.__evaluate_deltas()
        return not self.__function.constraint_on_find(self.__table['D'])

    def get_forward(self):
        return self.__forward

    def __getitem__(self, key):
        return self.__table[key]

    def get_plan(self):
        result = []
        for x in xrange(1, len(self.__table['A'])):
            if x in self.__table['B']:
                item = self.__table['A'][0][self.__table['B'].index(x)]
                result.append(round(item, 2))
            else:
                result.append(0)
        return result

if __name__ == '__main__':
    a = [[2, 6],
         [1, 2],
         [-1, 1],
         [1, 0],
         [0, 1]]
    b = [3, 4]
    c = [0, 3, 2, 0, 0]

    m_rank = ExtremumCondition.matrix_rank(matrix(a))
    a.append(b)
    m_expanded_rank = ExtremumCondition.matrix_rank(matrix(a))

    if m_rank != m_expanded_rank:
        print u'Система ограничений не совместима'
    else:
        a.pop()
        table = SimplexTable(MaxCondition(), A=a, B=b, C=c)
        try:
            while table.next_iteration():
                pass
            print table.get_plan()
        except UnboundedException:
            print U'Ошибка. Функция неограничена сверху'