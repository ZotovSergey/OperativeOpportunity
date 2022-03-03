import numpy as np
from math import pi
from itertools import product
from shapely import geometry
from shapefile import Reader
from Geography.EarthEllipsoid import EarthEllipsoid
from Polygons.Geometry import polygon_contains_points


class PolygonsGroup:
    """
    @Описание:
        Класс PolygonsGroup моделирует группу наземных тестовых полигонов, которые должны быть просканированны, для того
            чтобы поставленная задача была выполнена. Предусматривается возможность разбиение на мелкие сегменты всех
            полигонов в группе, подсчёта количества сканирований каждого сегмента в отдельности, подсчёта
            приблизительной площади всех полигонов, моделирования облачности над полигонами в общем или над каждым в
            отдельности.

    @Аргументы:
        earth_ellipsoid - объект класса EarthEllipsoid.

    @Поля:
        earth_ellipsoid - объект класса EarthEllipsoid, задающий параметры модели эллипсоида Земли, на котором
            расположена группа полигонов. При инициализации присваивается значение аргумента earth_ellipsoid
        polygons_list – список, состоящий из объектов класса Polygon, где каждый элемент обозначает тестовый полигон
            принадлежащий моделируемой группе полигонов. При инициализации - пустой список. Заполняется методом
            self.to_read_shape_file
        full_segments_list - список, содержащий все сегменты, на которые разделены все полигоны из моделируемой группы.
            При инициалиации задаётся пустым списком. Заполняется с помощью метода to_split_all_polygons
        full_area – общая площадь всех полигонов, входящих в моделируемую группу. Общая площадь вычисляется
            суммированием площадей всех полигонов, входящих в группу, то есть всех, перечисленных в списке (кв. м)
            self.polygons_list, а площадь каждого полигона – объекта класса Polygon записана в их поле polygons_area.
            При инициализации - 0. Вычисляется при применении метода self.to_split_all_polygons.
        common_cloudiness_distr - объект CloudDistrDataset содержащий распределения вероятности балла облачности в
            районе, где находится моделируемая группа полигонов для некоторых годовых периодов. По умолчанию облачность
            всегда отсутствует.
        common_current_cloudiness_in_score - балл облачности над моделируемой группой полигонов. При инициалиации
            присваивается 0. Обновляется методом to_randomize_common_cloudiness объекта класса PolygonGroup.
        common_calculations_of_cloudiness - если True, вычислять облачность над всей группой полигонов вместе, если
            False, то отдельно для каждого полигона в группе (boolean). По умолчанию True.

    @Методы:
        to_read_shape_file(self, shape_file_address) – читает из shape-файла информацию о тестовых полигонах, полигоны
            добавляются в список self.polygons_list.
        to_set_polygons_names(self, polygons_names) - задает имена полигонам, входящим в моделируемую группу.
        to_split_all_polygons(self, lat_fineness, long_fineness, earth_ellipsoid) – разбивает каждый полигон в группе на
            сегменты с мелкостью разбиения lat_fineness по широте и long_fineness по долготе, вычисляет приблизительную
            общую площадь всех полигонов в моделируемой группе.
        to_calc_percentages_of_grabbed_areas(self) – определяет ход выполнения тематической задачи – сколько процентов
            площади тестовых полигонов попало в поле зрения ГСК, сколько раз.
        to_clear_count_of_grabs(self) - очищает все сегменты всех полигонов от "сканирований". То есть после применения
            считается, что полигоны не сканировались.
        to_add_common_cloudiness_distr(self, common_cloud_distr) - добавляет распределения вероятности балла облачности
            в районе, где находится моделируемая группа полигонов для некоторых годовых периодов. Распределения
            записываются в виде объекта CloudDistrDataset.
        to_add_cloudiness_distr_to_each_polygon(self, cloudiness_distr_tables_list, cloudiness_distr_ranges_list) -
            добавляет распределения вероятности появления облачности некоторого балла для каждого полигона в
            моделируемой группе для годовых периодов, границы которых также записываются.
        to_normalize_distribution(distribution) - нормируетраспределение на единицу distribution.
        to_calculate_cloudiness_above_group(self, time) - определяет текущий балл облачности для всей моделируемой
            группы полигонов в соответвии с распределением из self.common_cloudiness_distr_table для времени из
            аргумента time.
    """

    def __init__(self, earth_ellipsoid=EarthEllipsoid()):
        self.earth_ellipsoid = earth_ellipsoid
        self.polygons_list = []
        self.centers_geo_coordinates = []
        self.radii_list = []
        self.full_segments_geo_coordinates_list = None
        self.splitted_polygons_list = []
        self.spitted_polygons_centers_geo_coordinates = []
        self.full_area = 0

    def to_read_shape_file(self, shape_file_address, polygons_names=None):
        """
        @Описание:
            Открывает для чтения shape-файл, в цикле обходится все геометрические объекты типа “Polygon” (в терминах из
                описания формата shape-файла), из записанных в открытый shape-файл. На каждом витке создаётся объект
                self.Polygon(shape, self.earth_ellipsoid), где shape - один из геометрических объектов типа “Polygon”,
                self.earth_ellipsoid - объект EarthEllipsoid, задающий параметры эллипсоида Земли, на котором находится
                полигон, и добавляется в self.polygonsList, задает стандартные имена.
        :param shape_file_address: адресс shape-файла, из которого читается информация о тестовых полигонах.
        :param polygons_names: список названий полигонов. Допустим неполный или пустой список, так как для оставшихся
            полигонов, которым не достается названия они прописываются автоматически. Допустимо None - заменяется на
            пустой список. По умолчанию None.
        """
        self.polygons_list = []
        self.centers_geo_coordinates = []
        self.radii_list = []
        if polygons_names is None:
            polygons_names = []
        # Чтение информации о полигонах из файла по адресу shape_file_address
        shape_file = Reader(shape_file_address)
        shape_file_list = list(shape_file.iterShapes())
        # Создание объектов self.Polygon(shape, self.earth_ellipsoid), где shape - один из объектов из списка, объект
        # EarthEllipsoid, задающий параметры эллипсоида Земли, на котором находится создаваемый полигон shape_file_list
        # и запись в self.polygons_list в цикле
        for shape in shape_file_list:
            self.polygons_list.append(Polygon(shape, self))
            self.centers_geo_coordinates.append(self.polygons_list[-1].center_geo_coordinate)
            # self.radii_list.append(self.polygons_list[-1].radius)
        # self.radii_list = np.array(self.radii_list)
        self.centers_geo_coordinates = np.array(self.centers_geo_coordinates)
        # Задание полигонам стандартных названий: "Полигон 1", "Полигон 2", "Полигон 3" ...
        self.to_set_polygons_names(polygons_names)

    def to_set_polygons_names(self, polygons_names):
        """
        Метод задает имена из списка polygons_names полигонам, входящим в моделируемую группу (self.polygons_list).
            Полигонам из списка self.polygons_list сообщаются имена из списка polygons_names в соответствии их номерам
            в списках. Полигонам, которым название не задано (None) или не хватает имен (если длина списка
            self.polygons_list больше polygons_names), сообщается стандартное название "Полигон 1", "Полигон 2",
            "Полигон 3" ...
        :param polygons_names: список имен, присвиваемых полигонам из списка self.polygons_list (в списке допустимы
            значения None)
        :return:
        """
        number_of_unnamed_polygon = 1
        i = 0
        while i < len(polygons_names):
            if polygons_names[i] is not None:
                self.polygons_list[i].to_set_polygons_name(polygons_names[i])
            else:
                self.polygons_list[i].to_set_polygons_name('Полигон ' + str(number_of_unnamed_polygon))
                number_of_unnamed_polygon += 1
            i += 1
        while i < len(self.polygons_list):
            self.polygons_list[i].to_set_polygons_name('Полигон ' + str(number_of_unnamed_polygon))
            number_of_unnamed_polygon += 1
            i += 1

    def to_split_polygons(self, lat_fineness, long_fineness):
        """
        @Описание:
            Метод разбивает каждый полигон в группе (из списка polygons_list) на сегменты с мелкостью разбиения
                lat_fineness по широте и long_fineness по долготе, вычисляет приблизительную общую площадь всех
                полигонов в моделируемой группе
        :param lat_fineness: мелкость разбиения сегментов по широте (км)
        :param long_fineness: мелкость разбиения сегментов по долготе (км)
        """
        self.splitted_polygons_list = []
        self.radii_list = []
        self.spitted_polygons_centers_geo_coordinates = []
        self.full_segments_geo_coordinates_list = np.empty((0, 3))
        self.full_area = 0
        # В цикле к каждому полигону применяется метод to_split_polygon(lat_fineness, long_fineness)
        for i, polygon in enumerate(self.polygons_list):
            splitted_pol = polygon.to_split_polygon(lat_fineness, long_fineness)
            # Сегменты, на которые разделился полигон из группы включаются в список всех сегментов моделируемой группы
            #   self.full_segments_list
            if splitted_pol is not None:
                self.spitted_polygons_centers_geo_coordinates.append(splitted_pol.center_geo_coordinate)
                self.radii_list.append(splitted_pol.radius)
                self.splitted_polygons_list.append(splitted_pol)
                self.full_segments_geo_coordinates_list = np.append(self.full_segments_geo_coordinates_list,
                                                                    splitted_pol.segments_geo_coordinates_list, axis=0)
                # Площадь полигона из группы прибавляется к общей площади полигонов из моделируемой группы
                self.full_area += splitted_pol.area
        self.spitted_polygons_centers_geo_coordinates = np.array(self.spitted_polygons_centers_geo_coordinates)
        self.radii_list = np.array(self.radii_list)

    def to_union_splitted_polygons(self, polygons_numbers=None, union_polygons_name=None):
        if polygons_numbers is None:
            polygons_numbers = range(0, len(self.splitted_polygons_list))

        if union_polygons_name is None:
            union_polygons_name = self.splitted_polygons_list[polygons_numbers[0]].name
            for pol_num in polygons_numbers[1:]:
                union_polygons_name = ", ".join([union_polygons_name,
                                                 self.splitted_polygons_list[polygons_numbers[pol_num]].name])

        union_segments_geo_coordinates = []
        union_segments_areas = []
        for pol_num in polygons_numbers:
            union_segments_geo_coordinates += list(self.splitted_polygons_list[pol_num].segments_geo_coordinates_list)
            union_segments_areas += list(self.splitted_polygons_list[pol_num].segments_area_list)
        union_segments_geo_coordinates = np.array(union_segments_geo_coordinates)
        union_segments_areas = np.array(union_segments_areas)

        union_polygon = SplittedPolygon(union_polygons_name, union_segments_geo_coordinates, union_segments_areas,
                                        self.earth_ellipsoid)

        spitted_polygons_centers_geo_coordinates = list(self.spitted_polygons_centers_geo_coordinates)
        radii_list = list(self.radii_list)
        for i in reversed(polygons_numbers):
            del self.splitted_polygons_list[i]
            del spitted_polygons_centers_geo_coordinates[i]
            del radii_list[i]
        self.splitted_polygons_list.append(union_polygon)
        spitted_polygons_centers_geo_coordinates.append(union_polygon.center_geo_coordinate)
        radii_list.append(union_polygon.radius)
        self.spitted_polygons_centers_geo_coordinates = spitted_polygons_centers_geo_coordinates
        self.radii_list = radii_list


class Polygon:
    """
    @Описание:
        Класс Polygon моделирует наземный полигон, который входит в группу полигонов, моделируемых PolygonsGroup, и
            должен быть просканирован для решени задачи. Предусматривается возможность разбиения полигона на мелкие
            сегменты, учёта количества сканирования каждого сегмента в отдельности, подсчёта приблизительной площади
            полигона, моделирования облачности над полигоном и каждым сегментом отдельно.

    @Аргументы
        shape - содержит геометрический объект типа “Polygon” (в терминах из описания shape-формата), прочитанный из
            shape-файла. Записанный геометрический объект определяет границы тестового полигона
        polygons_group - объект PolygonGroup - группа полигонов, к которой относится создаваемы объект

    @Поля:
        shape – содержит геометрический объект типа “Polygon” (в терминах из описания shape-формата), прочитанный из
            shape-файла. Записанный геометрический объект определяет границы тестового полигона. При инициализации
            объект в поле записывается значение аргумента shape.
        name - название полигона. При инициализации присваивается None. Задается методом to_set_polygons_name.
        segments_list – список, состоящий из объектов класса Segment, где каждый элемент обозначает один из сегментов на
            которые разбит моделируемый полигон. Заполняется при применении метода to_split_polygon.
        polygons_area – площадь моделируемого полигона с некоторой точностью, равна общей площади всех сегментов, на
            которые разбит тестовый полигон, то есть сумма полей segments_area объектов Segment, входящих в список
            self.segmentsList. При инициализации объекта приравнивается к нулю, значение обновляется при применении
            метода to_split_polygon.
        own_group - содержит объект PolygonsGroup - группу полигонов, к которой относится моделируемый полигон. В
            own_group содержится поле. Присваивается аргумент polygons_group при инициализации объекта.
            earth_ellipsoid, используещееся в вычислениях для моделируемого полигона.
        top_border_lat – широта самой северной точки на границе тестового полигона. Вычисляется при инициализации
        bot_border_lat – широта самой южной точки на границе тестового полигона. Вычисляется при инициализации
        left_border_long – долгота самой западной точки на границе тестового полигона. Вычисляется при инициализации
        right_border_long – долгота самой восточной точки на границе тестового полигона. Вычисляется при инициализации
        center – объект Coordinates.GeoCoordinatesAndPointSet, содержащий координаты точки на поверхности Земли. (высота
            над поверхностью Земли всегда 0) - центра полигона. Вычисляется при инициализации
        radius – расстояние от центра тестового полигона до самой дальней от него точки на границы
            полигона. Определяется при инициализации.
        cloudiness_distr - объект CloudDistrDataset содержащий распределения вероятности балла облачности в
            районе, где находится моделируемая группа полигонов для некоторых годовых периодов. По умолчанию облачность
            всегда отсутствует. Устанавливается методом to_add_cloudiness_distr_to_each_polygon объекта класса
            PolygonGroup.
        current_cloudiness_in_score - балл облачности над моделируемым полигоном. При инициалиации присваивается 0.
            Обновляется методом to_randomize_cloudiness_to_each_polygon объекта класса PolygonGroup.


    @Методы
        to_set_polygons_name(self, polygons_name) - задает название полигона.
        to_split_polygon(self, lat_fineness, long_fineness) - разделяет моделируемый полигон на сегменты с мелкостью
            разбиения по широте lat_fineness и по долготе long_fineness
        to_calculate_segments_area(self, lat_of_grids_nodes, longFineness) - метод вычисляет площади сегментов,
            на которые делится полигон в зависимости от широты их верхних и нижних границ lat_of_grids_nodes и мелкости
            разбиения по долготе long_fineness
        to_calculate_space_from_equator_to_lat(self, lat) - вычисляет площадь поверхности эллипсоида Земли
            self.own_group.earth_ellipsoid от экватора до заданной аргументом lat широты
        segment_is_hidden(self) – случайным образом определяет закрыт ли некоторый сегмент моделируемого полигона
            облаком или тенью от облака. Вероятность зависит от текущего балла облачности над полигоном
            self.current_cloudiness_in_score
        to_randomize_cloudiness_to_polygon(self, time) - добавляет распределения вероятности появления облачности
            некоторого балла моделируемого полигона для годовых периодов, границы которых также
            записываются
    """

    def __init__(self, shape, polygon_group):
        self.name = None
        self.own_group = polygon_group
        [self.left_long_border, self.bot_lat_border, self.right_long_border, self.top_lat_border] = shape.bbox
        self.center_geo_coordinate = np.array([(self.left_long_border + self.right_long_border) / 2,
                                               (self.top_lat_border + self.bot_lat_border) / 2,
                                               0])
        self.border_points = np.array(shape.points[:-1])
        # distances_to_point = self.own_group.earth_ellipsoid.dist_between_geo_coordinates(
        #     np.array([self.center_geo_coordinate] * len(self.border_points)), self.border_points)
        # self.radius = np.max(distances_to_point)

    def to_set_polygons_name(self, polygons_name):
        """
        Метод задает название полигона polygons_name.
        :param polygons_name: Название полигона.
        :return: поле self.name приравнивается к значению аргумента polygons_name.
        """
        self.name = polygons_name

    def to_split_polygon(self, lat_fineness, long_fineness):
        """
        @Описание:
            Разделяет моделируемый полигон на сегменты с мелкостью разбиения по широте lat_fineness и по долготе
                long_fineness
        :param lat_fineness: мелкость разбиения сегментов по широте (км)
        :param long_fineness: мелкость разбиения сегментов по долготе (км)
        """
        # Координаты центрального сегмента будут совпадать с центром полигона
        coordinates_of_central_segment = self.center_geo_coordinate
        # Определение координат центров сегментов в виде сетки -двумерного списка
        lat_of_segments = np.concatenate([np.flip(np.arange(
            coordinates_of_central_segment[1] - lat_fineness, self.bot_lat_border, -lat_fineness)),
            np.arange(coordinates_of_central_segment[1], self.top_lat_border, lat_fineness)])
        long_of_segments = np.concatenate([np.flip(np.arange(
            coordinates_of_central_segment[0] - long_fineness, self.left_long_border, -long_fineness)),
            np.arange(coordinates_of_central_segment[0], self.right_long_border, long_fineness)])
        # Определение широт границ сегментов с юга на север
        lat_of_grids_nodes = np.arange(lat_of_segments[0] - lat_fineness / 2,
                                       lat_of_segments[-1] + lat_fineness,
                                       lat_fineness)
        # Вычисление площадей сегментов в зависимости от их широты
        area_of_segments_of_lat = self.to_calculate_segments_area(lat_of_grids_nodes, long_fineness)

        potential_segments_geo_coordinates = np.array(list(product(long_of_segments, lat_of_segments, [0])))
        areas_grid = np.repeat(area_of_segments_of_lat, len(long_of_segments))

        segments_inx = np.where(self.contains_points(potential_segments_geo_coordinates))
        segments_geo_coordinates_list = potential_segments_geo_coordinates[segments_inx]
        segments_areas = areas_grid[segments_inx]

        if len(segments_geo_coordinates_list) > 0:
            return SplittedPolygon(self.name, segments_geo_coordinates_list, segments_areas,
                                   self.own_group.earth_ellipsoid)
        else:
            return None

    def to_calculate_segments_area(self, lat_of_grids_nodes, long_fineness):
        """
        @Описание:
            Метод вычисляет площади сегментов, на которые делится полигон в зависимости от широты их верхних и нижних
                границ lat_of_grids_nodes и мелкости разбиения по долготе long_fineness
        :param lat_of_grids_nodes: список границ сегментов полигона по широте (градусы)
        :param long_fineness: мелкость разбиения моделируемого полигона на сегменты по долготе (градусы)
        :return: список длиной в (len(lat_of_grids_nodes) - 1), содержащий площади сегментов (кв. км), на которые
            делится моделируемый полигон. В элементе с номером i содержит площадь сегмента, расположенного между
            широтами lat_of_grids_nodes[i] и lat_of_grids_nodes[i + 1]
        """
        # Создание пустого массива, в который будут записываться результаты
        area_of_segments_of_lat = []
        # Вычисление площадь эллипсоида self.own_group.earth_ellipsoid от экватора до нижней границы самого южного
        #   сегмента
        area_from_equator_to_segments = self.to_calculate_space_from_equator_to_lat(lat_of_grids_nodes)
        # В цикле вычислеяются площади сегментов от широты и добавляется в area_of_segments_of_lat
        for i in range(1, len(area_from_equator_to_segments)):
            # Вычисляется площадь некоторого сегмента и добавляется в список результатов
            area_of_segments_of_lat.append((area_from_equator_to_segments[i] - area_from_equator_to_segments[i - 1]) *
                                           long_fineness / 360)
        return np.array(area_of_segments_of_lat)

    def to_calculate_space_from_equator_to_lat(self, lat):
        """
        @Описание:
            Метод вычисляет площадь поверхности эллипсоида Земле self.own_group.earth_ellipsoid от экватора до заданной
                аргументом lat широты. Вычисление происходит методом интегрирования части эллипсоида по поверхности
        :param lat: широта, до которой вычесляется площадь (градусы)
        :return: площадь поверхности на эллипсоиде Земли self.own_group.earth_ellipsoid от экватора эллипсоида до
            заданной широты lat (кв. км)
        """
        semi_major_axis = self.own_group.earth_ellipsoid.semi_major_axis
        semi_minor_axis = self.own_group.earth_ellipsoid.semi_minor_axis
        y = semi_minor_axis * np.sin(lat * np.pi / 180)
        a = semi_major_axis * semi_major_axis - semi_minor_axis * semi_minor_axis
        b = a ** 0.5
        c = semi_minor_axis ** 4
        d = (c + y * y * a) ** 0.5
        return 2 * pi / semi_minor_axis * d * ((c * np.log(b * d + y * a)) / (b * d) + y)

    def contains_points(self, points):
        return polygon_contains_points(self.border_points, points)


class SplittedPolygon:
    def __init__(self, name, segments_geo_coordinates_list, segments_area_list, earth_ellipsoid):
        self.name = name
        self.segments_geo_coordinates_list = segments_geo_coordinates_list
        long_lat_seg = np.swapaxes(self.segments_geo_coordinates_list, 0, 1)
        long_seg = long_lat_seg[0]
        lat_seg = long_lat_seg[1]
        self.left_long_seg = min(long_seg)
        self.right_long_seg = max(long_seg)
        self.bot_lat_seg = max(lat_seg)
        self.top_lat_seg = min(lat_seg)
        self.center_geo_coordinate = np.array([(self.right_long_seg + self.left_long_seg) / 2,
                                               (self.bot_lat_seg + self.top_lat_seg) / 2])
        distances_to_segments_centers = earth_ellipsoid.dist_between_geo_coordinates(
            np.array([self.center_geo_coordinate] * len(self.segments_geo_coordinates_list)),
            self.segments_geo_coordinates_list)
        self.radius = np.max(distances_to_segments_centers)
        self.segments_area_list = segments_area_list
        self.area = np.sum(self.segments_area_list)
