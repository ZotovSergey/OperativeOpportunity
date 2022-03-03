import numpy as np
import pickle
from datetime import datetime, timedelta
from random import random
from itertools import product
from netCDF4 import Dataset
from Time.DaysNumber import to_determine_days_number_in_not_leap_year


def to_load_cloud_distr_dataset(save_address):
    """
    @ Описание:
        Метод возвращает объект CloudDistrDataset с данными о распределении балла облачности, сохраненный по адресу
            save_address.
    :param save_address: адресс файла, из которого загружается сохраненный объект CloudDistrDataset (String).
    :return: объект класса CloudDistrDataset, прочитанные из файла по адресу save_address
    """
    with open(save_address, "rb") as file:
        loaded_data = pickle.load(file)
    return loaded_data


class CloudDistrDataset:
    """
    @Описание:
        Объекты класса содержат информацию о распределениях облачности над некоторой областью для заданных в этом же
            объекте годовых интервалов времени. Объекты этого класса также содержат информацию о области, для которой
            определены распределения, и о данных, которые использовались для формирования этих распределений. Объекты
            класса используются для статистического определения облачности над полигонами, моделируемыми классом
            TaskCalculations.Polygons.Polygon. Для этого используются методы
            TaskCalculations.Polygons.Polygon.to_randomize_cloudiness_to_polygon и
            TaskCalculations.Polygons.to_randomize_cloudiness. Предусмотрена возможность сохранения в файл.
    @Поля:
        cloud_probability_by_intervals - двумерный список. Массив содержит значения дискретных распределений вероятности
            балла облачности для каждого годового интервала распределения, границы которых записаны в поле
            self.intervals. В столбцах - вероятность для каждого балла от 0 до некоторого максимального значения балла,
            который равен len(self.cloud_probability_by_intervals) - 1, в строках - распределение вероятности для
            каждого интервала между днями в self.intervals. Задаётся при инициализации аргументом
            cloud_probability_by_intervals. По умолчанию None. Если cloud_probability_by_intervals - None и
            self.intervals - None, то приравнивается к [[1, 0]], что означает, что облачности нет никогда.
        intervals - список, содержащий границы годовых интервалов (в номерах дней в невисокосном году), для которых
            определены распределения из поля self.cloud_probability_by_intervals. Всего интервалов len(self.intervals)
            - 1. Каждый интервал включает все дни между двумя соседними значениями в списке, включая значение правой
            границы. Задаётся при инициализации аргументом intervals. По умолчанию None. Если
            self.cloud_probability_by_intervals - None и intervals - None, то приравнивается к [0, 365], что означает,
            что распределение облачности одинакого весь год.
        north_lat - широта северной границы области, для которой определены распределения облачности. Задаётся при
            инициализации аргументом north_lat. По умолчанию None.
        south_lat - широта южной границы области, для которой определены распределения облачности. Задаётся при
            инициализации аргументом south_lat. По умолчанию None.
        west_long - долгота западной границы области, для которой определены распределения облачности. Задаётся при
            инициализации аргументом west_long. По умолчанию None.
        east_long - долгота восточной границы области, для которой определены распределения облачности. Задаётся
            при инициализации аргументом east_long. По умолчанию None.
        coordinates - список уникальных значений, содержащий значения координат точек, для которых была определена.
            статистика облачности в исходных данных, по которым считались распределения. Координаты записаны в формате
            (широта, долгота). Задаётся при инициализации аргументом coordinates. По умолчанию None.
        time - список уникальных значений, содержащий значения времени, в которое были определены статистические
            значения облачности в исходных данных, по которым считались распределения. Время записано в формате
            datetime. Задаётся при инициализации аргументом time. По умолчанию None.
    @Методы:
        to_save_data(self, save_name, save_directory) - сохраняет объект self по адрсу save_directory в файле с
            названием save_name и расширением .file.
    """
    def __init__(self, cloud_probability_by_intervals=None, intervals=None, coordinates=None, time=None):
        if cloud_probability_by_intervals is not None and intervals is not None:
            self.cloud_probability_by_intervals = cloud_probability_by_intervals
            self.intervals = intervals
        else:
            # если распределения или годовые интервалы для распределений не заданы, то задается одно распределение на
            #   весь год, по которому облачность всегда нулевая
            self.cloud_probability_by_intervals = [[1, 0]]
            self.intervals = [0, 365]
        # координаты границ области, для которой задано распределение
        self.north_lat = None
        self.south_lat = None
        self.west_long = None
        self.east_long = None
        self.coordinates = coordinates
        self.time = time

    def to_randomize_cloudiness(self, time):
        """
            @Описание:
                Метод вероятностным методом определяет балл облачности в соответвии с распределениями из cloud_distr для
                    заданных годовых интервалов времени для времени из аргумента time
            :param cloud_distr: двумерный массив, содержащий распределения вероятности выпадения некоторого
                балла облачности для годовых периодов, границы которых записаны в аргументе borders_of_ranges
            :param time: время, для которого определяется балл облачности (datetime)
            :return: случайный балл облачности
            """
        # Определение дня от начала года в невисокосном году
        timestamp = ((time - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's'))
        days_number_in_year = to_determine_days_number_in_not_leap_year(datetime.utcfromtimestamp(timestamp))
        # Проверяется, входит ли заданный день в заданные годовые периоды, для которых заданы распределения облачности.
        #   Если нет, то возвращается 0.
        if (days_number_in_year < self.intervals[0]) or (days_number_in_year > self.intervals[-1]):
            return 0
        # Определяется, в какой период входит время time
        i = 0
        while days_number_in_year <= self.intervals[i]:
            i += 1
        distribution = self.cloud_probability_by_intervals[i]
        # Вычисление случайного балла облачности в соответствии с распределением для определенного выше годового периода
        rand = random()
        j = 0
        sum_proportion = distribution[0]
        while rand >= sum_proportion:
            j += 1
            sum_proportion += distribution[j]
        return j

    def to_save_data(self, save_name, save_directory):
        """
        @Описание:
            Метод сохраняет объект на диск в виде файла с заданным именем в заданной директории. Объект может быть
                загружен с помощью метода to_load_data.
        :param save_name: название создаваемого файла сохранения (String)
        :param save_directory: адрес директории, в которую сохраняется файл сохранения (String)
        :return: файл сохранения, содержащий данный объект, который можно загрузить с помощью метода to_load_data.
            Название файла - save_name и он записывается в директории по адресу save_directory
        """
        with open("".join([save_directory, '\\', save_name, '.file']), "wb") as file:
            pickle.dump(self, file, pickle.HIGHEST_PROTOCOL)


class CloudDataset:
    """
    @Описание:
        Объекты класса содержат статистическую информацию об облачности (балл облачности) за некоторое время на
            некоторой сетке координат. Эти данные извлекаются из файла netcdf4, находящегося по адресу, заданному
            аргументом dataset_address. Также аргументом default_max_score задается предполагаемый максимальный балл
            облачности в этих данных.
    @Поля:
        lat - список широт сетки статистических данных по баллу облачности. При инициализации извлекается из
            файла netcdf4 по адресу, заданному аргументом dataset_address (float).
        long - список долготы сетки статистических данных по баллу облачности. При инициализации извлекается из файла
            netcdf4 по адресу, заданному аргументом dataset_address (float).
        time - список времени сетки статистических данных по баллу облачности. При инициализации извлекается из файла
            netcdf4 по адресу, заданному аргументом dataset_address (datetime).
        tcdc - трехмерный список. Сетка данных по баллу облачности. Сетка составлена на пересечении значений полей
            self.lat, self.long, self.time. При инициализации извлекается из файла netcdf4 по адресу, заданному
            аргументом dataset_address (float).
        max_score - Предполагаемый максимальный балл облачности, для данных в self.tcdc. Задается при инициализации
            аргументом default_max_score. По умолчанию - 100. (int)
    """
    def __init__(self, dataset_address, default_max_score = 100):
        dataset = Dataset(dataset_address)
        # извлечение данных из dataset
        self.lat = np.array(dataset['lat'][:])
        self.long = np.array(dataset['lon'][:])
        self.time = np.array(dataset['time'][:])
        self.tcdc = np.array(dataset['tcdc'][:])
        self.max_score = default_max_score


class CloudCalculator:
    """
    @Описание:
        Класс содержит методы, с помощью которых вычисляются распределения балла облачности в виде объектов
            CloudDistrDataset по наборам данных - объектам класса CloudDistrDataset
    @Поля:
        north_lat - северная граница области вычислений по широте. При инициализации - None. Вычисляется методом
            to_set_coordinates_of_calculation_area(float)
        south_lat - южная граница области вычислений по широте. При инициализации - None. Вычисляется методом
            to_set_coordinates_of_calculation_area (float)
        west_long - западная граница области вычислений по долготе. При инициализации - None. Вычисляется методом
            to_set_coordinates_of_calculation_area (float)
        east_long - восточная граница области вычислений по долготе. При инициализации - None. Вычисляется методом
            to_set_coordinates_of_calculation_area (float)
        max_score - максимальный балл облачности (при котором небо закрыто полностью). При инициализации - 10. Задается
            методом to_set_max_score (int)
        annual_intervals - список номеров дней в невисокосном году, ограничивающих периоды в году, для которых будет
            вычисляться средняя облачность. По умолчанию эти периоды - это 12 месяцев. Первое число день перед новым
            интервалом, второе - последнее число интервала. (int)
        cloud_datasets_list - список наборов данных, из которых берутся архивные данные по облачности (CloudDataset)
    @Методы:
        to_set_coordinates_of_calculation_area - вычисляет границы прямоугольной области, для которой вычисляется
            средняя облачность: северную и южную в градусах долготы, западную и восточную в градусах долготы.
        to_set_max_score - задает максимальный балл облачности (при котором небо закрыто полностью)
        to_set_annual_intervals - задает границы интервалов в невисокосном году в номерах дней в невисокосном году, для
            которых определяется средняя облачность.
        to_add_cloud_dataset - добавляет к данным, по которым будет вычисляться вероятностное распределение облачности
            еще один набор данных.
        to_calculate_cloud_distr - вычисляет дискретное распределение балла облачности для каждого заданного годового
            интервала по всем наборам данных, содержащихся в поле self.cloud_datasets_list и возвращает в виде объекта
            CloudDistrDataset.
        to_extract_cloud_data_from_one_dataset - вычисляет дискретное распределение балла облачности для каждого
            заданного годового интервала времени по одноу заданному аргументом набору данных, не нормированное, в виде
            объекта двумерного списка. Кроме того, возвращает списки уникальных координат точек, для которых были
            определены исходные данные, и времени исходных данных.
        to_extract_data_about_territory - извлекает из заданного набора данных данных для точек, входящих в область, для
            которой считаются распределения облачности.
        to_transfer_to_new_score - переводит быллы облачности с некоторым максимальным баллом old_max_score в баллы с
        максимально допустимым баллом self.max_score.
    @Константы:
        DEFAULT_INTERVALS_BY_MONTHS - границы интервалов в невисокосном году, для которых определяются средние значения
            облачности по умолчанию. Каждый интервал - месяц в невисокосном году.
        ZERO_TIME - время, которое в файлах netCDF4 вопринимается, как нулевое время
    """
    # границы интервалов в невисокосном году, для которых определяются средние значения облачности по умолчанию. Каждый
    #   интервал - месяц в невисокосном году.
    DEFAULT_INTERVALS_BY_MONTHS = [0, 31, 59, 90, 120, 151, 182, 212, 243, 273, 304, 334, 365]
    # время, которое в файлах netCDF4 вопринимается, как нулевое время
    ZERO_TIME = datetime(1800, 1, 1, 0)

    def __init__(self):
        # северная граница области вычислений по широте
        self.north_lat = None
        # южная граница области вычислений по широте
        self.south_lat = None
        # западная граница области вычислений по долготе
        self.west_long = None
        # восточная граница области вычислений по долготе
        self.east_long = None
        # максимальный балл облачности (при котором небо закрыто полностью)
        self.max_score = 10
        # список номеров дней в невисокосном году, ограничивающих периоды в году, для которых будет вычисляться средняя
        #   облачность. По умолчанию эти периоды - это 12 месяцев. Первое число день перед новым интервалом, вторя -
        #   последнее число интервала.
        self.annual_intervals = self.DEFAULT_INTERVALS_BY_MONTHS
        # список файлов, из которых берутся архивные данные по облачности
        self.cloud_datasets_list = []

    def to_set_coordinates_of_calculation_area(self, center_lat, center_long, lat_angle_diff, long_angle_diff):
        """
        @Описание:
            Метод вычисляет границы прямоугольной области, для которой вычисляется средняя облачность: северную и южную
                в градусах долготы, западную и восточную в градусах долготы. Для вычислений на вход подаются координаты
                центра предполагаемой области вычислений, её длина с юга на север в градусах широты и длина с запада на
                восток в градусах долготы.
        :param center_lat: широта центра области вычислений облачности (double)
        :param center_long: долгота центра области вычислений облачности (double)
        :param lat_angle_diff: длина области вычислений облачности с юга на север в градусах широты (double)
        :param long_angle_diff: длина области вычислений облачности с запада на восток в градусах долготы (double)
        :return: вычисляет self.north_lat, self.south_lat, self.west_long, self.east_long
        """
        # вычисление границ области вычисления облачности
        self.north_lat = center_lat + lat_angle_diff / 2
        self.south_lat = center_lat - lat_angle_diff / 2
        self.west_long = center_long - long_angle_diff / 2
        self.east_long = center_long + long_angle_diff / 2

    def to_set_max_score(self, max_score):
        """
        @Описание:
            Метод задает максимальный балл облачности (при котором небо закрыто полностью)
        :param max_score: максимальный балл облачности (int)
        :return: заполняет поле self.max_score
        """
        self.max_score = max_score

    def to_set_annual_intervals(self, annual_intervals):
        """
        @Описание:
            Метод задает границы интервалов в невисокосном году в днях, для которых определяется средняя облачность.
                Если на вход метода подается пустой список или None, то то задаются нтервалы по умолчанию - по месяцам.
                Если задается список с одним значением, то год елится на два интервала: от началагода до этого дня
                включительно и после заданного дня до конца года.
        :param annual_intervals: список границ интервалов в невисокосном году в днях (int)
        :return: заполняет поле self.annual_intervals
        """
        if annual_intervals is not None and len(annual_intervals) > 0:
            self.annual_intervals = annual_intervals
        elif annual_intervals is not None and len(annual_intervals) > 0:
            self.annual_intervals = annual_intervals
            self.annual_intervals.insert(0, self.DEFAULT_INTERVALS_BY_MONTHS[0])
            self.annual_intervals.appand(self.DEFAULT_INTERVALS_BY_MONTHS[-1])
        else:
            self.annual_intervals = self.DEFAULT_INTERVALS_BY_MONTHS

    def to_add_cloud_dataset(self, dataset_address, max_score=100):
        """
        @Описание:
            Метод добавляет к данным, по которым будет вычисляться вероятностное распределение облачности еще один набор
                данных. Предполагается, что каждый набор данных представляет архивные данные по облачности. Отдельно
                задается предполагаемый максимальный балл облачности.
        :param dataset_address: адрес нового набора данных (String).
        :param max_score: предполагаемый максимальный балл облачности, для данных, содержащихся в файле по адресу
            dataset_address. По умолчанию - 100. (int)
        :return: добавляет новый набор данных в поле self.datasets_list
        """
        self.cloud_datasets_list.append(CloudDataset(dataset_address, max_score))

    def to_calculate_cloud_distr(self):
        """
        @Описание:
            Метод вычисляет дискретное распределение балла облачности для каждого заданного годового интервала
                (заданного полем self.intervals) по всем наборам данных, содержащихся в поле self.cloud_datasets_list и
                возвращает в виде объекта CloudDistrDataset.
        :return: объект CloudDistrDataset, содержащий распределения балла облачности для каждого заданного годового
            интервала, вычисленный по наборам данных из self.cloud_datasets_list
        """
        # создание двумерного массива, заполненного нулевыми значениями. Массив будет содержать значения дискретных
        #   нормированных распределений вероятности балла облачности для каждого годового интервала распределения,
        #   полученного из всех заданных наборов данных. В столбцах - вероятность для каждого балла от 0 до max_score, в
        #   строках - распределение вероятности для каждого интервала между днями в annual_intervals
        cloud_distr_by_all_datasets = np.array([[0] * (self.max_score + 1)] * (len(self.annual_intervals) - 1))
        # список, в который будут записываться координаты точек, из которых берутся данные (только уникальные
        #   координаты) (tuple)
        data_coordinates = set()
        # список, в который будет записываться время данных, используемых в формировании распредения(только уникальные
        #   значения) (datetime)
        data_time = set()

        for k in range(0, len(self.cloud_datasets_list)):
            # вычисление распределения по одному набору данных
            cloud_data_from_one_dataset, coordinates_from_dataset, time_from_dataset = \
                self.to_extract_cloud_data_from_one_dataset(self.cloud_datasets_list[k])
            # добавление уникальных координат и времени данных в списки
            data_coordinates.update(coordinates_from_dataset)
            data_time.update(time_from_dataset)
            # добавление распределения по одному набору данных к общему распределению
            cloud_distr_by_all_datasets += cloud_data_from_one_dataset
        # распределение нормируется и возвращается
        sums_of_each_distr_for_annual_periods = np.sum(cloud_distr_by_all_datasets, axis=1)
        cloud_distribution = (cloud_distr_by_all_datasets.T / sums_of_each_distr_for_annual_periods).T
        return CloudDistrDataset(cloud_distribution, self.annual_intervals, data_coordinates, data_time)

    def to_extract_cloud_data_from_one_dataset(self, dataset):
        """
        @Описание:
            Метод вычисляет дискретное распределение балла облачности для каждого заданного годового интервала
                (заданного полем self.intervals) по набору данных dataset, не нормированное в виде объекта двумерного
                списка. Кроме того, возвращает списки уникальных координат точек, для которых были определены исходные
                данные, и времени исходных данных.
        :param dataset: набор данных, по которомы вычисляются распределения облачности (CloudDataset)
        :return:
            двумерный список, содержащий значения дискретных распределений вероятности балла облачности для каждого
                годового интервала распределения, границы которых записаны в поле self.annual_intervals. В столбцах -
                вероятность для каждого балла от 0 до self.max_score,, в строках - распределение вероятности для каждого
                интервала между днями в self.annual_intervals (int)
            список уникальных координат точек, для которых были определены исходные данные в виде: (широта, долгота)
                ((float, float))
            список уникальных времени исходных данных (datetime)
        """
        # создание двумерного массива, заполненного нулевыми значениями. Массив будет содержать значения дискретных
        #   не нормированных распределений вероятности балла облачности для каждого годового интервала распределения,
        #   полученного из набора данных dataset. В столбцах - вероятность для каждого балла от 0 до max_score, в
        #   строках - распределение вероятности для каждого интервала между днями в annual_intervals
        non_norm_cloud_distr_by_dataset = np.array([[0] * (self.max_score + 1)] * (len(self.annual_intervals) - 1))
        # извлечение данных по облачности на заданной территории из dataset, а также координат, в которых берутся данные
        data_for_territory, data_coordinates = self.to_extract_data_about_territory(dataset)
        # перевод баллов, в которых данные были в наборе данных dataset в балл, в которых будет распределение облачности
        data_for_territory_in_final_score = self.to_transfer_to_new_score(data_for_territory, dataset.max_score)
        # объявление массива, в который будет добавляться время данных (datetime)
        data_time = []
        # заполнение двумерного массива - не нормированных распределений вероятности балла облачности для каждого
        #   годового интервала распределения, полученного из набора данных dataset
        for i in range(0, len(dataset.time)):
            # извлечение из массива времени данных в часах от self.ZERO_TIME
            time_in_hours = dataset.time[i]
            # перевод времени данных из часов от self.ZERO_TIME в datetime
            time = self.ZERO_TIME + timedelta(hours=time_in_hours)
            # добавление time в список времени данных time_in_hours
            data_time.append(time)
            # перевод времени данных из datetime в дни в невисокосном году
            time_in_days = to_determine_days_number_in_not_leap_year(time)
            # данные по облачности за выбранное время
            score_of_time = data_for_territory_in_final_score[i]
            # расппределение всех измерений по годовым интервалам
            for j in range(0, len(self.annual_intervals) - 1):
                if self.annual_intervals[j] < time_in_days <= self.annual_intervals[j + 1]:
                    non_norm_cloud_distr_by_dataset[j][score_of_time] += 1
                    break
        return non_norm_cloud_distr_by_dataset, data_coordinates, data_time

    def to_extract_data_about_territory(self, dataset):
        """
        @Описание:
            Метод извлекает из набора данных dataset (CloudDataset) данных для точек, входящих в область, ограниченную
                координатами (self.north_lat, self.west_long) и (self.south_lat, self.east_long) и возвращает их в виде
                двумерного массива - списков значений балл облачности для каждого времени dataset.time, а также список
                уникальных координат точек, из которых берутся данные
        :param dataset: набор данных, из которого извлекаются данные (CloudDataset)
        :return:
            двумерный список - списки значений балла облачности для каждого времени dataset.time (float)
            список уникальных координат точек, из которых берутся данные в виде (широта, долгота) ((float, float))
        """
        # извлечение данных о координатах станций в данном наборе данных dataset
        lat_list = dataset.lat
        long_list = dataset.long
        # определение индексов значений сетки координат по долгоге и широте, которые входят в область, для которой
        #   считается распределение облачности
        lat_valid_indexes_list = np.where((lat_list >= self.south_lat) & (lat_list <= self.north_lat))[0]
        long_valid_indexes_list = np.where((long_list >= self.west_long) & (long_list <= self.east_long))[0]
        # извлечение данных о облачности для заданной территории
        tcdc = dataset.tcdc[:, lat_valid_indexes_list][:, :, long_valid_indexes_list]
        # составление списка координат
        data_coordinates = list(product(lat_valid_indexes_list, long_valid_indexes_list))
        # возвращение извлеченных данных и группировка всех значений облачности по времени, также возвращение всех
        #   координат
        return np.reshape(tcdc, (len(tcdc), len(tcdc[0]) * len(tcdc[0][0]))), data_coordinates

    def to_transfer_to_new_score(self, data_in_old_score, old_max_score):
        """
        @Описание:
            Метод переводит баллы облачности заданные аргументом data_in_old_score с максимальным баллом old_max_score
            в баллы с максимальным баллом self.max_score.
        :param data_in_old_score: баллы облачности с максимальным баллом old_max_score. Допустим numpy.array (float)
        :param old_max_score: максимальный балл облачности баллов из аргумента data_in_old_score
        :return: баллы облачности с максимально допустимым значением self.max_score (int)
        """
        # вычисление шага в старых баллах, соответствующего одному новому
        step_in_old_score = old_max_score / self.max_score
        # перевод в новые баллы и возвращение
        return (data_in_old_score / step_in_old_score).round().astype(int)
