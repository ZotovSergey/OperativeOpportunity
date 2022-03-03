import numpy as np
import pickle
import re
from math import pi, modf
from copy import deepcopy
from datetime import datetime, timedelta
import pandas as pd
from pyorbital import astronomy, tlefile, orbital
from Ballistics.Orbit import Orbital
from AnalyticGeometry import AnalyticGeometry
from Geography.Coordinates import cart_to_geo_coordinates, geo_to_cart_coordinates
from Time.DaysNumber import to_determine_date_by_days_number_in_not_leap_year,\
    to_determine_days_number_in_not_leap_year
from Time.TimeCounter import to_get_time_counts_for_period
from Time.TimeUnits import to_get_unit_in_seconds
from Polygons.Geometry import polygon_contains_points

def to_make_tle(a, e, i, big_omega, omega, M0, epoch,
                satnumber='00000', classification='U', id_launch_year='00', id_launch_number='000',
                id_launch_piece='A  ', mean_motion_derivative=' .00000000', mean_motion_sec_derivative=' 00000-0',
                bstar=' 00000-0', element_number=' 000', loop_num='00000'):
    # a - большая полуось (км), float
    # e - эксцентриситет, float
    # i - наклонение (°), float
    # big_omega - долгота восходящего узла (°), float
    # omega - аргумент перицентра (°), float
    # M0 - средняя аномалия (°), float
    # epoch - эпоха спутника, дата и время с точностью до мкс в формате datetime
    # sat_num - номер спутника в системе NORAD (не имеет значения), 5 цифр, int или str, например: '00000'
    # classification - классификация спутника (не имеет значения), 1 символ, str, пример: 'U' (U=Unclassified — не секретный)
    # id_launch_year - идентификатор года запуска - две последние цифры года (не имеет значения), 2 цифры, int или str, например: '00'
    # id_launch_number - идентификатор номера запуска - номер запуска в году (не имеет значения), 3 цифры, str, например: '000'
    # id_launch_piece - идентификатор части запуска (не имеет значения), 3 символа, str пример: 'A  '
    # mean_motion_derivative - производная среднего движения, (должно быть нулевым при вычислении оперативности), 8 цифр после точки, float или str, например: '-.00000000'
    # mean_motion_sec_derivative - вторая производная среднего движения, (должно быть нулевым при вычислении оперативности), 5 цифр после точки и отрицательная степень, float или str, например: ' .00000-0'
    # bstar - rоэффициент торможения, (должно быть нулевым при вычислении оперативности), 5 цифр после точки и отрицательная степень, float или str, например: ' .00000-0'
    # element_number (не имеет значения), 4 символа, int или str, например: ' 000'
    # loop_num - номер витка (не имеет значения), 5 цифр, int или str, например: '00000'

    # геоцентрическая гравитационная постоянная (km ** 3 / s ** 2)
    MU_EARTH = 398600.4415888888

    def int_to_str(int_value, digits_num, fill_symbol='0'):
        str_value = str(int_value)
        while len(str_value) < digits_num:
            str_value = fill_symbol + str_value
        return str_value

    def float_fractional_part_to_str(float_value, digits_num_after_point):
        if type(float_value) != type(''):
            format_id = '.' + str(digits_num_after_point) + 'f'
            round_float_value = format(float_value, format_id)
            if float_value >= 0:
                str_value = ' .' + str(round_float_value)[2:]
            else:
                str_value = '-.' + str(round_float_value)[3:]
        else:
            str_value = float_value
        return str_value

    def float_to_str(float_value, whole_digit_num, digits_num_after_point):
        frac, whole = modf(float_value)
        whole_str = int_to_str(int(whole), whole_digit_num, fill_symbol=' ')
        frac_str = float_fractional_part_to_str(frac, digits_num_after_point)[1:]
        return whole_str + frac_str

    def sci_float_to_str(float_value, num_of_visible_digits_after_point):
        if type(float_value) != type(''):
            if float_value >= 0:
                str_value = ' '
            else:
                str_value = '-'
            minus_power = len(re.split(r'[1-9]', format(float_value, '.15f')[2:])[0])
            format_id = '.' + str(minus_power + num_of_visible_digits_after_point) + 'f'
            str_value += format(abs(float_value), format_id)[2 + minus_power:] + '-' + str(minus_power)
        else:
            str_value = float_value
        return str_value

    def to_calc_control_sum(line):
        control_sum = 0
        for symb in line:
            if str.isdigit(symb):
                control_sum += int(symb)
            if symb == '-':
                control_sum += 1
        return control_sum % 10
    # Вачисление частоты обращения (виток / день)
    n = (24 * 60 * 60) / (2 * pi * np.sqrt((a ** 3) / MU_EARTH))

    # Номера строк
    line1 = '1 '
    line2 = '2 '
    # Обозначение спутника
    line1 += int_to_str(satnumber, 5) + classification + ' ' + int_to_str(id_launch_year, 2) + int_to_str(id_launch_number, 3) + id_launch_piece + ' '
    line2 += int_to_str(satnumber, 5) + ' '
    # Эпоха
    epoch_year = str(int(epoch.year % 100 / 10) * 10 + epoch.year % 10)
    epoch_day = format((epoch - datetime(epoch.year, 1, 1, 0, 0, 0)).total_seconds() / 24 / 3600, '.8f')
    line1 += epoch_year + epoch_day + ' '
    # Параметры, определяющие изменения орбиты со временем. Для вычислений оперативности следует оставлять нулевыми
    line1 += float_fractional_part_to_str(mean_motion_derivative, 8) + ' ' + sci_float_to_str(mean_motion_sec_derivative, 5) + ' ' + \
             sci_float_to_str(bstar, 5) + ' 0 '
    # Прочее
    line1 += str(element_number)
    # Контрольная сумма
    line1 += str(to_calc_control_sum(line1))

    # Элементы орбиты
    line2 += float_to_str(i, 3, 4) + ' ' + float_to_str(big_omega, 3, 4) + ' ' + float_fractional_part_to_str(e, 7)[2:] + ' ' + \
             float_to_str(omega, 3, 4) + ' ' + float_to_str(M0, 3, 4) + ' ' + float_to_str(n, 2, 8)

    # Прочее
    line2 += int_to_str(loop_num, 5)
    # Контрольная сумма
    line2 += str(to_calc_control_sum(line2))

    return line1, line2

class Satellite:
    """
    @Описание:
        Класс Satellite моделирует работу спутника ДЗЗ. Предусматривается возможность двигаться по заданной орбите,
            вычислять полосу захвата, какая часть заданных тестовых полигонов попала в полосу захвата при соответствии
            внешних условий заданным.

    @Аргументы:
        sat_name - название спутника.
        tle_address - адресс текстового файла, содержащего данные TLE для спутниа с названием sat_name. Допустимо
            значение None
        angle_of_view - угол обзора гиперспектрометра, базирующегося на спутнике с названием sat_name (градусы).
        satellites_group - объект SatelliteGroup - группа спутников, в который входит моделируемый спутник.
            Предполагается, что создаваемый объект Satellite инициируется методом именно этого объекта.

    @Поля:
        sat_name - название спутника. Название должно совпадать с названием спутниа в системе NORAD. Присваивается
            значение аргумента sat_name при инициализации.
        orbit - объект Orbital, с помощью которого вычисляется положение моделируемого спутника на орбите Земли,
            заданной данными TLE, расположенных по адрессу из аргумента tle_address для спутника под названием sat_name
            в системе NORAD. Если значение аргумента sat_name - None, то данные TLE загружаются с веб-сайта
            celestrak.com для спутника под названием из аргумента sat_name. Если аргумент sat_name не задан или задан
            неправильно, то определение объекта Orbital невозможно. Создается при инициализации.
        satellite_coordinate_set - объект SatelliteCoordinateSet, содержащий важнейшие координаты моделируемого
            спутника. При инициалиции - None. Обновляется методом to_move_to_time.
        angle_of_view - половина угла обзора гиперспектрометра в радианах, базирующегося на моделируемом спутнике.
            Задается аргументом angle_of_view разделенным пополам при инициализации в градусах и переводится в радианы
            (радианы).
        satellites_group - объект SatelliteGroup, обозначающий спутниковую группировку, в которую входит моделируемый
            спутник. Задается аргументом satellite_group при инициализации.
        scanned_territory_for_last_step - список географических координат (GeoCoordinates) - вершины прямоугольника -
            часть полосы захвата, просканированные за последний шаг модельного времени. По умолчанию - пустой список.
            Определяется методом to_determine_close_polygons(self) и может очищаться методом to_act.
        close_polygons - список полигонов, достаточно близких к моделируемому спутнику, чтобы был шанс быть
            просканированными. По умолчанию - пустой список. Задаются методом to_determine_close_polygons.

    @Методы:
        to_act(self, next_simulation_time) - моделирует работу спутника от текущего модельного времени спутниковой
            группы до времени, заданного аргументом. Возвращает площадь просканированных тестовых полигонов за время
            моделирования.
        to_move_to_time(self, next_time) - определяет координаты спутника в заданное время и присваивает их значения
            моделируемому спутнику.
        to_determine_scan_area(self, previous_coordinate, current_coordinate) - вычисление точки пересечения эллипсоида
            Земли и прямой повернутой вокруг оси движения спутника (касательной к поверхности Земли в подспутниковой
            точке) с заданными координатами (параллельной Земле) на заданный угол.
        to_calculate_geo_coordinates_of_ground_point_in_field_of_view(self, satellite_coordinates_set,
                                                                      angle_of_rotation) - вычисляет точку пересечения
            эллипсоида Земли и прямой повернутой вокруг оси движения спутника (проекция вектора скорости спутника на
            касательную плоскость к поверхности Земли в подспутниковой точке) с координатами спутника на заданный угол.
        to_determine_close_polygons(self) - определяет, какие полигоны из заданных, то есть тех, которые должны быть
            просканированны для решения задачи, находятся достаточно близко к моделируемому спутнику, чтобы части их
            территории могли попасть в полосу захвата моделируемого спутника.
        to_scan(self, polygons_to_scan) - моделирует процесс съемки заданных тестовых полигонов за некоторое время.
            Производится проверка того, попадают ли сегменты близких полигонов в полосу захвата за некоторое время,
            допустимый ли в момент съемки зенитный угол Солнца, моделирует облачность или для всех полигонов, или для
            каждого отдельно, задает случайную облачность для каждого сегмента и проверяет, допустима ли облачность над
            полигоном для съемки и не закрыт ли сегмент облаками. Если все эти условия соблюдены, то сегмент считается
            просканированным и его площадь прибавляется к сумме, подаваемой на выход метода.
    """
    def __init__(self, sat_name, angle_of_view, infinite_orbit=True, tle_address=None, tle_line1=None, tle_line2=None):
        self.sat_name = sat_name
        # Извлечение данных TLE из файла по адресу tle_address или загрузка с celestrak.com для спутника sat_name, если
        if tle_address is not None:
            tle = tlefile.read(sat_name, tle_address)
            tle1 = tle.line1
            tle2 = tle.line2
        else:
            tle1 = tle_line1
            tle2 = tle_line2
        # Создание объекта Orbital из пакета pyorbital
        self.orbit = Orbital(self.sat_name, line1=tle1, line2=tle2)
        # Перевод в радианы
        self.half_angle_of_view = pi * angle_of_view / 360
        self.area_of_interest = None
        self.canon_ellipsoid = None
        self.initial_simulation_time = None
        self.final_simulation_time = None
        self.step = None
        self.initial_annual_observation_period = None
        self.final_annual_observation_period = None
        self.min_sun_angle_above_horizon = None
        self.captured_segments_of_polygons = None
        self.time_of_capturing = None
        self.interpreted_segments_of_polygons = None
        self.time_of_interpreted_capturing = None
        self.infinite_orbit = infinite_orbit

    def to_set_area_of_interest(self, area_of_interest):
        self.area_of_interest = area_of_interest
        self.canon_ellipsoid = AnalyticGeometry.CanonicalEllipsoid(
            self.area_of_interest.earth_ellipsoid.semi_major_axis,
            self.area_of_interest.earth_ellipsoid.semi_major_axis,
            self.area_of_interest.earth_ellipsoid.semi_minor_axis)
        self.captured_segments_of_polygons = []
        self.time_of_capturing = []

    def to_begin_calculation(self, initial_simulation_time, final_simulation_time, step=1,
                             initial_annual_observation_period=1, final_annual_observation_period=365,
                             one_calc_period_unit='days', one_calc_period_units_number=1, dat_aging_seasons_count=None,
                             track_note_file=None):
        self.captured_segments_of_polygons = []
        self.time_of_capturing = []
        self.initial_simulation_time = initial_simulation_time
        self.final_simulation_time = final_simulation_time
        self.step = step
        self.initial_annual_observation_period = initial_annual_observation_period
        self.final_annual_observation_period = final_annual_observation_period
        self.to_calculate(self.initial_simulation_time, self.final_simulation_time,
                          one_calc_period_unit=one_calc_period_unit,
                          one_calc_period_units_number=one_calc_period_units_number, track_note_file=track_note_file)

    def to_continue_calculation(self, final_simulation_time,
                                one_calc_period_unit='days', one_calc_period_units_number=1):
        new_initial_time = self.final_simulation_time
        # self.initial_simulation_time = new_initial_time
        self.final_simulation_time = final_simulation_time
        self.to_calculate(new_initial_time, self.final_simulation_time,
                          one_calc_period_unit=one_calc_period_unit,
                          one_calc_period_units_number=one_calc_period_units_number)

    def to_calculate(self, initial_simulation_time, final_simulation_time,
                     one_calc_period_unit='days', one_calc_period_units_number=1, track_note_file=None):
        observation_period_inside_one_year = self.final_annual_observation_period > \
                                             self.initial_annual_observation_period
        one_period_sec = one_calc_period_units_number * to_get_unit_in_seconds(one_calc_period_unit)
        #     Начало периода
        calc_periods_beginning = initial_simulation_time
        # Проверяется, входет ли начальное время в допустимый период наблюдения. Если нет, то это время переносится на
        #   начало следующего периода наблюдений и именно оно считается начальным модельным временем и записывается в
        #   отчетах и выходных данных
        # Подсчет номера дня в невисокосном году начального модельного времени спутниковой группировки
        days_number_in_year = to_determine_days_number_in_not_leap_year(calc_periods_beginning)
        if (observation_period_inside_one_year ^
            ((days_number_in_year >= self.initial_annual_observation_period) and
             (days_number_in_year <= self.final_annual_observation_period))) or \
                (days_number_in_year == self.initial_annual_observation_period) or \
                (days_number_in_year == self.final_annual_observation_period):
            new_year = initial_simulation_time.year
            if (observation_period_inside_one_year and
                (days_number_in_year > self.final_annual_observation_period)) or not observation_period_inside_one_year:
                new_year += 1
            calc_periods_beginning = to_determine_date_by_days_number_in_not_leap_year(
                self.initial_annual_observation_period, new_year)
        # Год конца периода
        annual_period_ending_year = calc_periods_beginning.year
        if observation_period_inside_one_year and \
                to_determine_days_number_in_not_leap_year(calc_periods_beginning) > \
                self.final_annual_observation_period:
            annual_period_ending_year += 1
        annual_periods_ending = to_determine_date_by_days_number_in_not_leap_year(
            self.final_annual_observation_period, annual_period_ending_year)

        while (calc_periods_beginning < final_simulation_time):
            # Создание массива отсчетов времени за заданный период
            calc_periods_ending = calc_periods_beginning + timedelta(seconds=one_period_sec - 1)
            simulation_time_counts = to_get_time_counts_for_period(calc_periods_beginning,
                                                                   calc_periods_ending,
                                                                   self.step)
            self.to_act(simulation_time_counts, track_note_file=track_note_file)
            calc_periods_beginning = calc_periods_beginning + timedelta(seconds=one_period_sec)

            if calc_periods_beginning >= annual_periods_ending:
                if observation_period_inside_one_year:
                    annual_periods_beginning = to_determine_date_by_days_number_in_not_leap_year(
                        self.initial_annual_observation_period, annual_period_ending_year + 1)
                else:
                    annual_periods_beginning = to_determine_date_by_days_number_in_not_leap_year(
                        self.initial_annual_observation_period, annual_period_ending_year)
                annual_period_ending_year = annual_periods_beginning.year
                if observation_period_inside_one_year and \
                        to_determine_days_number_in_not_leap_year(calc_periods_beginning) > \
                        self.final_annual_observation_period:
                    annual_period_ending_year += 1
                annual_periods_ending = to_determine_date_by_days_number_in_not_leap_year(
                    self.final_annual_observation_period, annual_period_ending_year)
                calc_periods_beginning = annual_periods_beginning + timedelta(seconds=self.step)

    def to_act(self, simulation_time_counts, track_note_file=None):
        """
        @Описание:
            Метод моделирует работу спутника от текущего модельного времени группы
                (self.satellite_group.simulation_time) до времени, заданного аргументом next_simulation_time. В процесс
                моделирования входит моделиросвание движения моделируемого спутника, его полосы захвата, процесса
                съемки. Возвращает площадь просканированных тестовых полигонов за время моделирования в кв. м.
        :param next_simulation_time: время, до которого проходит моделирование. Следует использовать время не далекое от
            модельного, так как спутник между двумя координатами, сответствующих начальному и конечному времени
            моделирования, будет двигаться прямолинейно.
        :return: площадь просканированных тестовых полигонов за время моделирования (кв. м).
        """
        # Вычисление координат спутника за весь заданный период
        (pos_x_counts, pos_y_counts, pos_z_counts), (vel_x_counts, vel_y_counts, vel_z_counts) = \
            self.orbit.get_position(simulation_time_counts, normalize=False, infinite_orbit=self.infinite_orbit)
        # Моделирование движения спутника
        track = self.to_make_track(np.array(np.swapaxes(np.array([pos_x_counts, pos_y_counts, pos_z_counts]), 0, 1)),
                                   np.array(np.swapaxes(np.array([vel_x_counts, vel_y_counts, vel_z_counts]), 0, 1)),
                                   simulation_time_counts)
        if track_note_file is not None:
            file = open(track_note_file, 'a')
            for coord in track.geo_coordinates[:, 0:2]:
                file.write(str(coord[0]) + '\t' + str(coord[1]) + '\n')
            file.close()
        tracks_counts_inx, close_polygons_for_tracks_counts_inx = self.to_determine_close_polygons(track)
        # Проверка, есть ли вблизи тестовые полигоны
        scanned_polygons = self.to_determine_scan_area(track, tracks_counts_inx)
        # file = open('D:\\Documents\\Конференция Леса Евразии\\Test2.txt', 'w')
        # for i, inx in enumerate(tracks_counts_inx):
        #     file.write(str(track.geo_coordinates[inx][0]) + '\t' + str(track.geo_coordinates[inx][1]) + '\t' +
        #                str(scanned_polygons[i][3][0]) + '\t' + str(scanned_polygons[i][3][1]) + '\t' +
        #                str(scanned_polygons[i][2][0]) + '\t' + str(scanned_polygons[i][2][1]) + '\t' +
        #                str(scanned_polygons[i][1][0]) + '\t' + str(scanned_polygons[i][1][1]) + '\t' +
        #                str(scanned_polygons[i][0][0]) + '\t' + str(scanned_polygons[i][0][1]) + '\n')
        # file.close()
          # Моделируется сканирование близких полигонов и возвращается площадь

        self.to_scan(scanned_polygons, close_polygons_for_tracks_counts_inx, track.utc_time[tracks_counts_inx])
        self.to_clear_interpreted_segments_of_polygons()

    def to_make_track(self, cart_coordinates, velocity_vectors, time):
        """
        @Описание:
            Метод определяет координаты спутника во время next_time и присваивает их значения моделируемому спутнику.
        :param next_time: время в формате UTC в которое определяются координаты моделируемого спутника.
        :return: координаты в объекте класса SatelliteCoordinatesSet записываются в self.satellite_coordinates_set
        """
        return SatellitesTrack(cart_coordinates, velocity_vectors, time, self.area_of_interest.earth_ellipsoid)

    def to_determine_subsatellite_track(self, track, indexes):
        return track.subsatellite_coordinates(indexes, self.area_of_interest.earth_ellipsoid)

    def to_determine_close_polygons(self, track):
        """
        @Описание:
            Метод определяет, какие полигоны из заданных, то есть тех, которые должны быть просканированны для решения
                задачи (определенных в объекте класса Task, записанных в списке
                self.satellites_group.task.polygons_group), находятся достаточно близко к моделируемому спутнику, чтобы
                части их территории могли попасть в полосу захвата моделируемого спутника.
        :return: близкие к подспутниковой точке полигоны записываются в список близких к моделируемуму спутнику
                 полигонов.
        """
        sat_geo_coordinates_matrix = np.swapaxes(
            np.array([track.geo_coordinates] * len(self.area_of_interest.splitted_polygons_list)), 0, 1)
        centers_geo_coordinates_matrix = np.array([self.area_of_interest.spitted_polygons_centers_geo_coordinates] *
                                                  len(track.geo_coordinates))
        distances_to_polygons_centers = self.area_of_interest.earth_ellipsoid.dist_between_geo_coordinates(
            sat_geo_coordinates_matrix, centers_geo_coordinates_matrix)

        radii_matrix = np.array([self.area_of_interest.radii_list] * len(track.geo_coordinates))
        swath_matrix = np.swapaxes(np.array([track.geo_coordinates[:, 2]] *
                                            len(self.area_of_interest.splitted_polygons_list)), 0, 1) * \
                       np.tan(self.half_angle_of_view)
        step_path_matrix = np.swapaxes(np.array([(np.sum(track.velocity_vector ** 2, axis=1)) ** 0.5 * self.step] *
                                                len(self.area_of_interest.splitted_polygons_list)), 0, 1)
        near_distance_matrix = (1.1 * radii_matrix + swath_matrix + step_path_matrix)
        tracks_counts_close_polygons_inx = np.where(distances_to_polygons_centers < near_distance_matrix)
        tracks_counts_for_close_polygons_inx = tracks_counts_close_polygons_inx[0]
        close_polygons_inx = tracks_counts_close_polygons_inx[1]
        # if len(tracks_counts_for_close_polygons_inx):
        tracks_counts_inx = np.unique(tracks_counts_for_close_polygons_inx)
        close_polygons_for_tracks_counts_inx = \
            np.array([close_polygons_inx[np.where(tracks_counts_for_close_polygons_inx == track_count_inx)] for
                      track_count_inx in tracks_counts_inx])
        return tracks_counts_inx, close_polygons_for_tracks_counts_inx

    def to_determine_scan_area(self, track, current_coordinates_inx):
        """
        @Описание:
            Определяет полосу захвата моделируемого спутника при прямолинейном движении спутника между точками с
                координатами current_coordinates и next_coordinates. Результат записывается в поле
                self.scanned_territory_for_last_step в виде списка геосоординат (GeoCoordinates).
        :param current_coordinates: координаты моделируемого спутника в начальной точке (SatelliteCoordinatesSet).
        :param next_coordinates: координаты моделируемого спутника в конечной точке (SatelliteCoordinatesSet).
        :return: полоса захвата записывается в поле self.scanned_territory_for_last_step в виде объекта
        shapely.geometry.Polygon.
        """
        counts_inx = []
        if len(current_coordinates_inx) > 0:
            counts_inx = [1]
            for i in range(1, len(current_coordinates_inx)):
                if current_coordinates_inx[i] - current_coordinates_inx[i - 1] <= 1:
                    counts_inx.append(i + 1)

        polygons_counts_inx = np.unique(np.append(current_coordinates_inx, current_coordinates_inx - 1))
        if -1 in polygons_counts_inx:
            track = self.track_with_previous_value(track)
            polygons_counts_inx += 1
            current_coordinates_inx += 1
        # Радиус-вектор спутника
        sat_cart_coordinates_under_polygons_counts = track.cartesian_coordinates[polygons_counts_inx]
        sat_velocity_vectors_under_polygons_counts = track.velocity_vector[polygons_counts_inx]
        time = track.utc_time[polygons_counts_inx]
        #   Вычисление вектора, указывающего направление надира для спутника
        subsatellite_geo_coordinates, subsatellite_cart_coordinates = \
            self.to_determine_subsatellite_track(track, polygons_counts_inx)
        sat_nadir_vector_from_polygons_counts = subsatellite_cart_coordinates - \
                                                sat_cart_coordinates_under_polygons_counts
        # Вычисление границ полосы захвата для текущего положения спутника
        left_swath_borders = self.to_calculate_geo_coordinates_of_view(
            sat_cart_coordinates_under_polygons_counts, sat_velocity_vectors_under_polygons_counts, time,
            sat_nadir_vector_from_polygons_counts, self.half_angle_of_view)
        right_swath_borders = self.to_calculate_geo_coordinates_of_view(
            sat_cart_coordinates_under_polygons_counts, sat_velocity_vectors_under_polygons_counts, time,
            sat_nadir_vector_from_polygons_counts, -self.half_angle_of_view)
        # Определение углов прямоугольной полосы захвата с помощью метода
        #   self.to_calculate_geo_coordinates_of_ground_point_in_field_of_view
        # counts_inx = []
        # for i in range(1, len(polygons_counts_inx)):
        #     if polygons_counts_inx[i] - polygons_counts_inx[i - 1] <= 1:
        #         counts_inx.append(i)
        scanned_polygons = np.array([[left_swath_borders[i - 1],
                                      right_swath_borders[i - 1],
                                      right_swath_borders[i],
                                      left_swath_borders[i]] for i in counts_inx])
        return scanned_polygons

    def to_calculate_geo_coordinates_of_view(self, sat_cart_coordinates, velocity, time, nadirs, angle_of_rotation):
        """
        @Описание:
            Вычисление координаты точки на Земле, в которую направлен вектор повернутый от надира спутника в сторону от
            его движения параллельно эллипсоиду Земли на заданный угол.
        :param sat_vector: координаты спутника (Coordinates.CartesianCoordinates)
        :param velocity_vector: вектор скорости спутника (AnalyticGeometry.Vector)
        :param nadir_vector: вектор-надир для спутника (Vector)
        :param angle_of_rotation: угол поворота в градусах (int, double)
        :return: Широту и долготу искомой точки (int, double) (в скобках для того, чтобы записать в объект Polygon)
        """
        # Вычисление координат вектора - оси вращения, вокруг которой будет вращаться вектор-надир, чтобы найти вектор
        #   движения спутника параллельно плоскости эллипсоида Земли
        sat_vectors = AnalyticGeometry.Vector(sat_cart_coordinates)
        nadirs_vectors = AnalyticGeometry.Vector(nadirs)
        velocity_vectors = AnalyticGeometry.Vector(velocity)
        rot_axises = nadirs_vectors * velocity_vectors
        # Вычисление координат вектора с помощью вращения движения спутника параллельно плоскости эллипсоида Земли
        parallel_mov_vectors = nadirs_vectors.to_rotate_vector(rot_axises, 90)
        # Вычисление вектора, направленного в искомую точку с помощью вращения
        vector_to_searched_point = nadirs_vectors.to_rotate_vector(parallel_mov_vectors, angle_of_rotation)
        # Вычисление координат искомой точки по пересечения vector_to_searched_point и
        #   self.satellites_group.canon_ellipsoid
        searched_point = AnalyticGeometry. \
            to_found_point_of_intersection_of_line_and_canonical_ellipsoid_nearest_to_stating_point_of_line(
                AnalyticGeometry.Line(sat_vectors, vector_to_searched_point), self.canon_ellipsoid)
        # # Перевод декартовых  координат в географческие и вывод в виде объекта geometry.Point
        geo_coordinates_of_search_point = cart_to_geo_coordinates(
            np.swapaxes(np.array([searched_point.x, searched_point.y, searched_point.z]), 0, 1),
            time, self.area_of_interest.earth_ellipsoid)
        return geo_coordinates_of_search_point

    def track_with_previous_value(self, track):
        previous_time = track.utc_time[0] - np.timedelta64(self.step, 's')
        (previous_coordinate_x, previous_coordinate_y, previous_coordinate_z),\
        (previous_velocity_x, previous_velocity_y, previous_velocity_z) = \
            self.orbit.get_position(previous_time, normalize=False)
        previous_value = self.to_make_track(
            np.array([[previous_coordinate_x, previous_coordinate_y, previous_coordinate_z]]),
            np.array([[previous_velocity_x, previous_velocity_y, previous_velocity_z]]),
            np.array([previous_time]))
        return previous_value + track

    def to_scan(self, scanned_areas, close_polygons_for_tracks_counts_inx, time):
        """
        @Описание:
            Метод моделирует процесс съемки близких полигонов self.close_polygons) за некоторое время. Производится
                проверка того, попадают ли сегменты близких полигонов (объекты из списка segments_list объекта Polygon)
                в полосу захвата за некоторое время self.scanned_territory_for_last_step, допустимый ли в момент съемки
                self.satellites_group.simulation_time зенитный угол Солнца (меньше
                self.satellites_group.task.max_zenith_angle), моделирует облачность или для всех полигонов, или для
                каждого отдельно, задает случайную облачность для каждого сегмента и проверяет, допустима ли облачность
                над полигоном для съемки (меньше self.satellites_group.task.max_cloud_score) и не закрыт ли сегмент
                облаками. Если все эти условия соблюденыЮ, то сегмент считается просканированным и его площадь
                прибавляется к сумме, подаваемой на выход метода.
        :return: площадь просканированной территории (кв. м)
        """
        # if len(self.time_of_capturing) > 0:
        #     print(self.time_of_capturing[-1])
        for i, scanned_area in enumerate(scanned_areas):
            captured_segments_of_polygons = [[]] * len(self.area_of_interest.splitted_polygons_list)
            for polygon_inx in close_polygons_for_tracks_counts_inx[i]:
                segments_of_polygon_captured = polygon_contains_points(scanned_area,
                                                                       np.array(self.area_of_interest.
                                                                                splitted_polygons_list[polygon_inx].
                                                                                segments_geo_coordinates_list))
                captured_segments_of_polygons[polygon_inx] = list(np.where(segments_of_polygon_captured)[0])
            if sum([len(i) for i in captured_segments_of_polygons]) > 0:
                self.captured_segments_of_polygons.append(captured_segments_of_polygons)
                self.time_of_capturing.append(time[i])

    def to_consider_solar_angle(self, min_solar_angle_above_horizon):
        for count_inx in range(0, len(self.interpreted_segments_of_polygons)):
            for pol_inx in range(0, len(self.interpreted_segments_of_polygons[count_inx])):
                seg_inx = 0
                while seg_inx < len(self.interpreted_segments_of_polygons[count_inx][pol_inx]):
                    if (90 - astronomy.sun_zenith_angle(
                            self.time_of_interpreted_capturing[count_inx],
                            self.area_of_interest.splitted_polygons_list[pol_inx].
                                    segments_geo_coordinates_list[seg_inx][0],
                            self.area_of_interest.splitted_polygons_list[pol_inx].
                                    segments_geo_coordinates_list[seg_inx][1])) < \
                            min_solar_angle_above_horizon:
                        del self.interpreted_segments_of_polygons[count_inx][pol_inx][seg_inx]
                    else:
                        seg_inx += 1
        count_inx = 0
        while count_inx < len(self.interpreted_segments_of_polygons):
            if sum([len(i) for i in self.interpreted_segments_of_polygons[count_inx]]) > 0:
                count_inx += 1
            else:
                del self.interpreted_segments_of_polygons[count_inx]
                del self.time_of_interpreted_capturing[count_inx]

    def to_consider_cloudiness(self, max_cloudiness_score, cloudiness_distribution):
        if type(cloudiness_distribution) != type([]):
            common_cloudiness = True
        else:
            common_cloudiness = False
        last_day = 0
        last_year = 0
        cloudiness_scores_for_polygons = []
        for count_inx in range(0, len(self.interpreted_segments_of_polygons)):
            timestamp = ((self.time_of_interpreted_capturing[count_inx] - np.datetime64('1970-01-01T00:00:00')) /
                         np.timedelta64(1, 's'))
            current_datetime = datetime.utcfromtimestamp(timestamp)
            current_day = to_determine_days_number_in_not_leap_year(current_datetime)
            current_year = current_datetime.year
            if (current_day != last_day) or (current_year != last_year):
                if common_cloudiness:
                    cloudiness_score = cloudiness_distribution.to_randomize_cloudiness(
                        self.time_of_interpreted_capturing[count_inx])
                    cloudiness_scores_for_polygons = [cloudiness_score] * len(self.area_of_interest.polygons_list)
                else:
                    cloudiness_scores_for_polygons = []
                    for cloud_dist in cloudiness_distribution:
                        cloudiness_scores_for_polygons.append(
                            cloud_dist.to_randomize_cloudiness(self.time_of_interpreted_capturing[count_inx]))
            for pol_inx in range(0, len(self.interpreted_segments_of_polygons[count_inx])):
                if cloudiness_scores_for_polygons[pol_inx] > max_cloudiness_score:
                    self.interpreted_segments_of_polygons[count_inx][pol_inx] = []
            last_day = current_day
            last_year = current_year
        count_inx = 0
        while count_inx < len(self.interpreted_segments_of_polygons):
            if sum([len(i) for i in self.interpreted_segments_of_polygons[count_inx]]) > 0:
                count_inx += 1
            else:
                del self.interpreted_segments_of_polygons[count_inx]
                del self.time_of_interpreted_capturing[count_inx]

    def to_clear_interpreted_segments_of_polygons(self):
        self.interpreted_segments_of_polygons = self.captured_segments_of_polygons
        self.time_of_interpreted_capturing = self.time_of_capturing

    def to_save(self, save_name, save_directory):
        """
        # @Описание:
        #     Метод сохраняет объект self на диск в виде файла с заданным именем в заданной директории. Объект может быть
        #         загружен с помощью метода to_load_data.
        # :param save_name: название создаваемого файла сохранения (String)
        # :param save_directory: адрес директории, в которую сохраняется файл сохранения (String)
        # :return: файл сохранения, содержащий данный объект, который можно загрузить с помощью метода to_load_data.
        #     Название айла - save_name и он записывается в директории по адресу save_directory
        """
        with open("".join([save_directory, '\\', save_name, '.file']), "wb") as file:
            pickle.dump(self, file, pickle.HIGHEST_PROTOCOL)

def to_load_satellite(save_address):
    """
    @ Описание:
        Метод возвращает объект Task с данными о задаче, пролетах, выполнении задач, сохраненный по адресу save_address.
    :param save_address: адресс файла, из которого загружается сохраненный объект OutputDataMaker (String).
    :return: объект класса OutputDataMaker, содержащий данные о задаче, пролетах, выполнении задач, прочитанные из файла
        по адресу save_address
    """
    with open(save_address, "rb") as file:
        loaded_data = pickle.load(file)
    return loaded_data


class SatellitesTrack:
    def __init__(self, cartesian_coordinates, velocity_vectors, utc_time, earth_ellipsoid):
        self.cartesian_coordinates = cartesian_coordinates
        # Вычисление координат спутника в географической системе координат путем перевода координат
        #   cartesian_coordinates из прямоугольной экватериальной системы координат в географическую
        self.geo_coordinates = cart_to_geo_coordinates(cartesian_coordinates, utc_time, earth_ellipsoid)
        self.velocity_vector = velocity_vectors
        self.utc_time = utc_time

        # file = open('D:\\Documents\\Конференция Леса Евразии\\Test5.txt', 'w')
        # for i, coord in enumerate(self.geo_coordinates):
        #     file.write(str(coord[0]) + '\t' + str(coord[1]) + '\t' + str(coord[2]) + '\t' + str(np.sqrt(self.velocity_vector[i][0] ** 2 + self.velocity_vector[i][1] ** 2 + self.velocity_vector[i][2] ** 2)) + '\n')
        # file.close()

    def __add__(self, other):
        res = deepcopy(self)
        res.cartesian_coordinates = np.append(res.cartesian_coordinates, other.cartesian_coordinates, axis=0)
        res.geo_coordinates = np.append(res.geo_coordinates, other.geo_coordinates, axis=0)
        res.velocity_vector = np.append(res.velocity_vector, other.velocity_vector, axis=0)
        res.utc_time = np.append(res.utc_time, other.utc_time)
        return res

    def subsatellite_coordinates(self, indexes, earth_ellipsoid):
        # Вычисление координат подспутниковой точки в географической системе координат путем обнуления высоты для поля
        #   self.geo_coordinates
        subsatellite_geo_coordinates = np.swapaxes(np.array([self.geo_coordinates[indexes, 0],
                                                             self.geo_coordinates[indexes, 1],
                                                             np.zeros(len(indexes))]), 0, 1)
        # Вычисление координат подспутниковой точки в прямоугольной экватериальной системе координат путем перевода
        #   координат self.subsatellite_geo_coordinates из географической системы координат в прямоугольную
        #   экватериальную
        subsatellite_cartesian_coordinates = geo_to_cart_coordinates(subsatellite_geo_coordinates,
                                                                     self.utc_time[indexes],
                                                                     earth_ellipsoid)
        return subsatellite_geo_coordinates, subsatellite_cartesian_coordinates
