__author__ = 'bod'
# -*- coding: utf-8 -*-
import abc
from simplexengine import *
from numpy.linalg import svd
from numpy import sum, where


class ExtremumCondition(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def get_forward_row(self, delta_list):
        pass

    @abc.abstractmethod
    def compare(self, case, value):
        pass

    @abc.abstractmethod
    def constraint_on_unbounded(self, row_items):
        pass

    @abc.abstractmethod
    def constraint_on_find(self, row_items):
        pass

    @staticmethod
    def matrix_rank(a, tol=1e-8):
        """:return matrix rank"""
        s = svd(a, compute_uv=False)
        return sum(where(s > tol, 1, 0))


class MaxCondition(ExtremumCondition):
    def get_forward_row(self, delta_list):
        return delta_list.index(min(delta_list))

    def constraint_on_unbounded(self, row_items):
        return len([x for x in row_items if x <= 0]) == len(row_items)

    def constraint_on_find(self, delta_list):
        return min(delta_list) >= 0

    def compare(self, case, value):
        return case < value


class UnboundedException(Exception):
    pass