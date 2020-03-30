# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 13:37:56 2020

@author: Duan Yutong (dyt@physics.bu.edu)
"""
import os
import psycopg2


class DBConn:
    host = os.getenv('$DOS_TELEMETRY_DB_HOST', 'desi-db')
    port = os.getenv('$DOS_TELEMETRY_DB_PORT', '5442')
    db = os.getenv('$DOS_TELEMETRY_DB_NAME', 'desi_dev')
    user = os.getenv('$DOS_TELEMETRY_DB_USER', 'desi_reader')
    pw = os.getenv('$DOS_TELEMETRY_DB_PASSWORD', 'reader')
    conn = psycopg2.connect(host=host, port=port, database=db,
                            user=user, password=pw)  # DB connection


# class Singleton(type):
#     _instances = {}

#     def __call__(cls, *args, **kwargs):
#         if not cls._instances:
#             cls._instances[None] = super(Singleton, cls).__call__(
#                 *args, **kwargs)
#         return cls._instances[None]


# class DBConnSingleton(DBConn, metaclass=Singleton):
#     pass
