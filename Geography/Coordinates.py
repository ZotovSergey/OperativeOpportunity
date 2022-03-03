from shapely.geometry import Point
from pyorbital import astronomy
import numpy as np
from math import pi


def geo_to_cart_coordinates(geo_coordinates, utc_time, earth_ellipsoid):
    """
    Метод переводит кординаты из географической системы координат (self.long, self.lat, self.alt) в
        геоцентрическую прямоугольную экваториальную систему с учетом эллипсоида Земли ellipsoid и времени в формате
        UTC utc_time. Формула взята из лекций http://lnfm1.sai.msu.ru/grav/russian/lecture/tfe/node3.html
    :param utc_time: время в формате UTC, показывающее, как повернута Земля в данный момент. Если его значение None,
        то вращение Земли не учитывается
    :param earth_ellipsoid: эллипсоид Земли. Его форма учитывается при переводе
    :return: объект CartesianCoordinates - координаты (self.long, self.lat, self.alt), переведенные в географическую
        систему координат
    """
    # a - большая полуось данного эллипсоид (в плоскости экватора)
    a = earth_ellipsoid.semi_major_axis
    # b - малая полуось данного эллипсоид (вдоль оси z)
    b = earth_ellipsoid.semi_minor_axis
    # Проверка того, задано ли время
    if utc_time is not None:
        # Если время задано, то определяется угол поворота Земли в данный момент времени
        earth_turn = astronomy.gmst(utc_time)
    else:
        earth_turn = 0
    # Перевод долготы и широты из градусов в радианы и из подвижной системы в неподвижную
    long_lat_alt = np.swapaxes(geo_coordinates, 0, 1)
    long = np.deg2rad(long_lat_alt[0]) + earth_turn
    lat = np.deg2rad(long_lat_alt[1])
    alt = long_lat_alt[2]
    # p вычисляется для удобства вычислений
    p = a ** 2 / np.sqrt(a ** 2 * np.cos(lat) ** 2 + b ** 2 * np.sin(lat) ** 2)
    # Вычисление прямоугольных координат
    x = (p + alt) * np.cos(lat) * np.cos(long)
    y = (p + alt) * np.cos(lat) * np.sin(long)
    z = ((b ** 2 / a ** 2) * p + alt) * np.sin(lat)
    # Возвращает объет CartesianCoordinates с переведенными координатами
    return np.swapaxes(np.array((x, y, z)), 0, 1)


def cart_to_geo_coordinates(cart_coordinates, utc_time, earth_ellipsoid):
    """
    Метод переводит кординаты из геоцентрической прямоугольной экваториальной системы (self.x, self.y, self.z) в
        географическую с учетом эллипсоида Земли earth_ellipsoid и времени в формате UTC utc_time. Формула взята из
        библиотеки pyorbital (https://github.com/pytroll/pyorbital/blob/master/pyorbital/orbital.py)
    :param utc_time: время в формате UTC, показывающее, как повернута Земля в данный момент. Если его значение None,
        то вращение Земли не учитывается
    :param earth_ellipsoid: эллипсоид Земли. Его форма учитывается при переводе
    :return: объект GeoCoordinates - координаты (self.x, self.y, self.z), переведенные в географическую систему
        координат
    """
    f = earth_ellipsoid.f
    a = earth_ellipsoid.semi_major_axis
    # Проверка того, задано ли время
    if utc_time is not None:
        # Если время задано, то определяется угол поворота Земли в данный момент времени
        earth_turn = astronomy.gmst(utc_time)
    else:
        earth_turn = 0
    # Вычисление долготы (в радианах) с учетом поворота Земли (в подвижной системе координат)
    long_lat_alt = np.swapaxes(cart_coordinates, 0, 1)
    x = long_lat_alt[0]
    y = long_lat_alt[1]
    z = long_lat_alt[2]
    long = ((np.arctan2(y, x) - earth_turn) % (2 * pi))
    # Проверка того, не выходит ли долгота за пределы диапазона (-pi, pi]
    long = np.where(long > pi, long - pi * 2, long)
    long = np.where(long <= -pi, long + pi * 2, long)
    # Расстояние до координат из начала отсчета в плоскости экватора
    r = np.sqrt(x ** 2 + y ** 2)
    # Вычисление широты в сферической системе
    lat = np.arctan2(z, r)
    # Вычисление эксцентриситета эллипсоида Земли в квадрате
    e2 = f * (2 - f)
    # Вычисление долготы (в радианах) в геодезической системе координат (нормальная проекция координаты на заданный
    #   эллипсоид) с точностью до 1e-10 км
    while True:
        lat_in_last_iter = lat
        sin_lat_in_last_iter = np.sin(lat_in_last_iter)
        c = 1 / (np.sqrt(1 - e2 * (sin_lat_in_last_iter ** 2)))
        lat = np.arctan2(z / a + c * e2 * sin_lat_in_last_iter, r / a)
        if np.all(abs(lat - lat_in_last_iter) < 1e-10):
            break
    # Вычисление наименьшенго расстояния до поверхности заданного эллипсоида Земли
    alt = r / np.cos(lat) - c * a
    # Возвращает объет GeoCoordinate с переведенными координатами, долгота и широта при этом переводится в градусы
    return np.swapaxes(np.array((np.rad2deg(long), np.rad2deg(lat), alt)), 0, 1)
