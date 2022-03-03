import numpy as np
from math import pi


class EarthEllipsoid:
    """
    @Описание:
        Класс EarthEllipsoid моделирует земной эллипсоид с заданной большой полуосью и сжатием. По умолчанию (если
            хотя бы один из аргументов semi_major_axis, f не задан) задаётся модель земного эллипсоида из
            World Geodetic System 1984 (WGS-84) с большим радиусом 6378137 м сжатием - 1 / 298.257223563. Также
            предусматривается возможность приблизительного вычисления расстояния между двумя точками с известными
            координатами на сфере с радиусом, соответсвующим большой полуоси моделируемого эллипсоида

    @Аргументы:
        semi_major_axis - большая полуось модели земного эллипсоида (км)
        f - сжатие модели земного эллипсоида

    @Константы
        SEMI_MAJOR_AXIS_WGS_84 - большая полуось модели эллипсоида Земли из WGS-84
        F_WGS_84 - сжатие модели эллипсоида Земли из WGS-84

    @Поля класса:
        semi_major_axis - большая полуось модели земного эллипсоида (км). Полю присваивается значение аргумента
            semi_major_axis, если аргумент semi_major_axis не задан (имеет значение None), полю semi_major_axis
            присваивается начение константы SEMI_MAJOR_AXIS_WGS_84
        f - сжатие модели земного эллипсоида. Полю присваивается значение аргумента f, если аргумент semi_major_axis
            не задан (имеет значение None), полю f присваивается начение константы F_WGS_84
        semi_major_axis - малая полуось модели земного эллипсоида (км). Значение вычисляется при инициализации объекта
            из аргументов semi_major_axis, если аргумент semi_major_axis не задан (имеет значение None), полю
            semi_major_axis присваивается начение константы SEMI_MAJOR_AXIS_WGS_84

    @Методы класса:
        dist_between_geo_coordinates - вычисляет расстояние по кратчайшей дуге между двумя географическими координатами
    """

    SEMI_MAJOR_AXIS_WGS_84 = 6378.137
    F_WGS_84 = 1 / 298.257223563

    def __init__(self, semi_major_axis=None, f=None):
        # Производится проверка, заданы ли аргументы
        if (semi_major_axis is None) or (f is None):
            # Если аргументы не заданны, то объект будет моделировать модель земного эллипсоида WGS-84
            self.semi_major_axis = self.SEMI_MAJOR_AXIS_WGS_84
            self.f = self.F_WGS_84
        else:
            self.semi_major_axis = semi_major_axis
            self.f = f
        # Вычисление малой полуоси модели эллипсоида Земли
        self.semi_minor_axis = self.semi_major_axis * (1 - self.f)

    def dist_between_geo_coordinates(self, geo_coordinates_1, geo_coordinates_2):
        """
        @Описание:
            Метод вычисляет расстояние по кратчайшей дуге между двумя географическими координатами, записанными в
            аргументы geo_coordinate_1 и geo_coordinate_2. При вычислении считаем, что Земля - шар с радиусом
            self.semi_major_axis.Вычисление производится по модифицированной формуле гаверсинусов
        :param geo_coordinates_1: объект класса Coordinates.GeoCoordinate, содержащий координаты некоторой точки
        :param geo_coordinates_2: объект класса Coordinates.GeoCoordinate, содержащий координаты некоторой точки
        :return: расстояние между точками с координатами geo_coordinate_1 и geo_coordinate_2 по кратчайшей дуге в
            километрах
        """
        mat_dimension = len(geo_coordinates_1.shape)
        long_lat_alt_1 = np.swapaxes(geo_coordinates_1, 0, mat_dimension - 1)
        long_lat_alt_2 = np.swapaxes(geo_coordinates_2, 0, mat_dimension - 1)
        # Перевод градусов широты в радианы
        lat1 = pi * long_lat_alt_1[1] / 180
        lat2 = pi * long_lat_alt_2[1] / 180
        # Вычисление разницы между градусами долготы первой и второй точек и перевод в радианы
        delta_long = pi * (long_lat_alt_1[0] - long_lat_alt_2[0]) / 180
        # Вычисления значений тригоометрических функций
        sin_lat1 = np.sin(lat1)
        cos_lat1 = np.cos(lat1)
        sin_lat2 = np.sin(lat2)
        cos_lat2 = np.cos(lat2)
        sin_delta_long = np.sin(delta_long)
        cos_delta_long = np.cos(delta_long)
        # Вучисление угла между точками по модифицированной формуле гаверсинусов
        delta_small_angle = np.arctan((((cos_lat2 * sin_delta_long) ** 2 +
                                        (cos_lat1 * sin_lat2 - sin_lat1 * cos_lat2 * cos_delta_long) ** 2) ** 0.5) /
                                      (sin_lat1 * sin_lat2 + cos_lat1 * cos_lat2 * cos_delta_long))
        # Приведение угла delta_small_angle к положительному значению
        delta_small_angle = np.where(delta_small_angle < 0, delta_small_angle + pi, delta_small_angle)
        return np.swapaxes(delta_small_angle, 0, mat_dimension - 2) * self.semi_major_axis

    def to_str(self, count_of_numerals_after_point_in_semi_major_axis=3, count_of_numerals_after_point_in_f=8):
        """
        @Описание:
            Вывод данных о эллипсоиде Земли - длину большой полуоси (км) и сжатие
        :param count_of_numerals_after_point_in_semi_major_axis: количество знаков после точки при выводе длины большой
            полуоси (в километрах). По умолчанию 3.
        :param count_of_numerals_after_point_in_f: количество знаков после точки при выводе сжатия. По умолчанию 8.
        :return:
        """
        return "".join(["Большая полуось:\t",
                        str(round(self.semi_major_axis, count_of_numerals_after_point_in_semi_major_axis)),
                        " км\nСжатие:\t\t\t\t",
                        str(round(self.f, count_of_numerals_after_point_in_f))])
