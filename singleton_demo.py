# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 12:58:11 2020

@author: Duan Yutong (dyt@physics.bu.edu)
"""


class Singleton(type):
    _instances = {}

    def __call__(cls, *args, **kwargs):
        petal_id = kwargs.get('petal_id', None)
        if petal_id not in cls._instances:
            cls._instances[petal_id] = super(Singleton, cls).__call__(
                *args, **kwargs)
        return cls._instances[petal_id]


class PosMoveDB:
    def __init__(self, petal_id=None):
        print('PosMoveDB called, petal_id =', petal_id)
        self.petal_id = petal_id


class DBSingleton(PosMoveDB, metaclass=Singleton):
    pass


if __name__ == "__main__":
    for petal_id in [0, 0, 1, None]:
        db = DBSingleton(petal_id=petal_id)
        print(petal_id, db.petal_id)
