# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 13:37:56 2020

@author: Duan Yutong (dyt@physics.bu.edu)
"""
import psycopg2


class DBConn:
    conn = psycopg2.connect(
        host="desi-db", port="5442", database="desi_dev",
        user="desi_reader", password="reader")  # DB connection


# class Singleton(type):
#     _instances = {}

#     def __call__(cls, *args, **kwargs):
#         if not cls._instances:
#             cls._instances[None] = super(Singleton, cls).__call__(
#                 *args, **kwargs)
#         return cls._instances[None]


# class DBConnSingleton(DBConn, metaclass=Singleton):
#     pass
