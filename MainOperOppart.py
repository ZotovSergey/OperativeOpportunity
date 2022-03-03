# import numpy as np
# import math
# from datetime import datetime
# from Geography import EarthEllipsoid
# from Satellite.Satellite import Satellite, to_load_satellite, to_make_tle
# from Polygons import Polygons
# from Weather import Cloudiness
#
#
# if __name__ == '__main__':
#     MU_EARTH = 398600.4415888888
#     i = 51.6416
#     big_omega = 247.4627
#     e = 0.0006703
#     omega = 130.5360
#     M0 = 325.0288
#     n = 10
#     a = np.power((24 * 60 * 60 * np.sqrt(MU_EARTH)) / (2 * np.pi * n), 2 / 3)
#     line1, line2 = to_make_tle(a, e, i, big_omega, omega, M0, datetime(2018, 6, 5, 1, 0, 0, 0))
#     start = datetime.now()
#     earth_ellipsoid = EarthEllipsoid.EarthEllipsoid()
#     polygons_group = Polygons.PolygonsGroup(earth_ellipsoid)
#     #   Загрузка полигонов из файла
#     polygons_group.to_read_shape_file(
#         "D:\\Data\\Shapefiles\\Валуйское лесничество\\WGS_region (union).shp")
#     polygons_group.to_split_polygons(0.0025, 0.005)
#     polygons_group.to_union_splitted_polygons()
#     # sat = Satellite('ISS (ZARYA)', 3.5, tle_address='D:\\Data\\TLE\\tle.txt')
#     sat = Satellite('ISS (ZARYA)', 3.5, tle_line1=line1, tle_line2=line2)
#     sat.to_set_area_of_interest(polygons_group)
#     initial_annual_observation_period = 121
#     final_annual_observation_period = 273
#     initial_simulation_time = datetime(2018, 6, 5, 0, 0, 0)
#     final_simulation_time = datetime(2018, 6, 5, 1, 0, 0)
#     # file = open('D:\\Documents\\Конференция Леса Евразии\\Test.txt', 'w')
#     # for pol_inx in range(0, len(sat.area_of_interest.polygons_list)):
#     #     for seg_inx in range(0, len(sat.area_of_interest.polygons_list[pol_inx].segments_geo_coordinates_list)):
#     #         file.write(
#     #             str(sat.area_of_interest.polygons_list[pol_inx].segments_geo_coordinates_list[seg_inx][0]) + '\t' +
#     #             str(sat.area_of_interest.polygons_list[pol_inx].segments_geo_coordinates_list[seg_inx][1]) + '\n')
#     # file.close()
#     sat.to_begin_calculation(initial_simulation_time, final_simulation_time,
#                              initial_annual_observation_period=initial_annual_observation_period,
#                              final_annual_observation_period=final_annual_observation_period, track_note_file='D:\\Data\\file2.txt')
#     sat.to_save('test2118', 'D:\\Data')
#     # final_simulation_time = datetime(2318, 6, 5, 0, 0, 0)
#     # sat.to_continue_calculation(final_simulation_time)
#     # sat.to_save('test2318', 'D:\\results')
#     #
#     # final_simulation_time = datetime(2418, 6, 5, 0, 0, 0)
#     # sat.to_continue_calculation(final_simulation_time)
#     # sat.to_save('test2418', 'D:\\results')
#     #
#     # final_simulation_time = datetime(2518, 6, 5, 0, 0, 0)
#     # sat.to_continue_calculation(final_simulation_time)
#     # sat.to_save('test2518', 'D:\\results')
#     #
#     # final_simulation_time = datetime(2618, 6, 5, 0, 0, 0)
#     # sat.to_continue_calculation(final_simulation_time)
#     # sat.to_save('test2618', 'D:\\results')
#     #
#     # final_simulation_time = datetime(2718, 6, 5, 0, 0, 0)
#     # sat.to_continue_calculation(final_simulation_time)
#     # sat.to_save('test2718', 'D:\\results')
#     #
#     # final_simulation_time = datetime(2818, 6, 5, 0, 0, 0)
#     # sat.to_continue_calculation(final_simulation_time)
#     # sat.to_save('test2818', 'D:\\results')
#     #
#     # final_simulation_time = datetime(2918, 6, 5, 0, 0, 0)
#     # sat.to_continue_calculation(final_simulation_time)
#     # sat.to_save('test2918', 'D:\\results')
#     #
#     # final_simulation_time = datetime(3018, 6, 5, 0, 0, 0)
#     # sat.to_continue_calculation(final_simulation_time)
#     # sat.to_save('test3018', 'D:\\results')
#
#     #file = open('D:\\Documents\\Конференция Леса Евразии\\Test3.txt', 'w')
#     #for count_inx in range(0, len(sat.interpreted_segments_of_polygons)):
#     #    for pol_inx in range(0, len(sat.interpreted_segments_of_polygons[count_inx])):
#     #        for seg_inx in sat.interpreted_segments_of_polygons[count_inx][pol_inx]:
#     #            file.write(str(sat.area_of_interest.splitted_polygons_list[pol_inx].segments_geo_coordinates_list[seg_inx][0]) + '\t' +
#     #                       str(sat.area_of_interest.splitted_polygons_list[pol_inx].segments_geo_coordinates_list[seg_inx][1]) + '\n')
#     #file.close()
#
#     print(datetime.now() - start)

# Sentinel. Бронницкое лесничество
import numpy as np
import math
from datetime import datetime
from Geography import EarthEllipsoid
from Satellite.Satellite import Satellite, to_load_satellite, to_make_tle
from Polygons import Polygons
from Weather import Cloudiness


if __name__ == '__main__':
    start = datetime.now()
    earth_ellipsoid = EarthEllipsoid.EarthEllipsoid()
    polygons_group = Polygons.PolygonsGroup(earth_ellipsoid)
    #   Загрузка полигонов из файла
    # polygons_group.to_read_shape_file(
    #    'D:/Проекты/Оперативность/Бронницкое лесничество/Shapes/выпуклый контур Броннецкого лесничества.shp')
    polygons_group.to_read_shape_file(
       'D:/Проекты/Оперативность/Бронницкое лесничество/Shapes/выпуклый контур Савватьевского лесничества.shp')
    polygons_group.to_split_polygons(0.0025, 0.005)
    polygons_group.to_union_splitted_polygons()
    tle_path = 'D:/Проекты/Оперативность/Бронницкое лесничество/TLE/resurs_p3_tle.txt'
    #sat = Satellite('RESURS P3', 3.62, tle_address=tle_path)
    sat = Satellite('RESURS P3', 53.28, tle_address=tle_path)
    sat.to_set_area_of_interest(polygons_group)
    initial_annual_observation_period = 120
    final_annual_observation_period = 243
    initial_simulation_time = datetime(2021, 5, 1, 0, 0, 0)
    final_simulation_time = datetime(2031, 5, 1, 0, 0, 0)
    sat.to_begin_calculation(initial_simulation_time, final_simulation_time,
                             initial_annual_observation_period=initial_annual_observation_period,
                             final_annual_observation_period=final_annual_observation_period,
                             track_note_file='D:/Проекты/Оперативность/Бронницкое лесничество/Промежуточные данные/track_file.txt')
    sat.to_save('resurs_p3_mult_2031_s', 'D:/Проекты/Оперативность/Бронницкое лесничество/Промежуточные данные')
    # final_simulation_time = datetime(2318, 6, 5, 0, 0, 0)
    # sat.to_continue_calculation(final_simulation_time)
    # sat.to_save('test2318', 'D:\\results')
    #
    # final_simulation_time = datetime(2418, 6, 5, 0, 0, 0)
    # sat.to_continue_calculation(final_simulation_time)
    # sat.to_save('test2418', 'D:\\results')
    #
    # final_simulation_time = datetime(2518, 6, 5, 0, 0, 0)
    # sat.to_continue_calculation(final_simulation_time)
    # sat.to_save('test2518', 'D:\\results')
    #
    # final_simulation_time = datetime(2618, 6, 5, 0, 0, 0)
    # sat.to_continue_calculation(final_simulation_time)
    # sat.to_save('test2618', 'D:\\results')
    #
    # final_simulation_time = datetime(2718, 6, 5, 0, 0, 0)
    # sat.to_continue_calculation(final_simulation_time)
    # sat.to_save('test2718', 'D:\\results')
    #
    # final_simulation_time = datetime(2818, 6, 5, 0, 0, 0)
    # sat.to_continue_calculation(final_simulation_time)
    # sat.to_save('test2818', 'D:\\results')
    #
    # final_simulation_time = datetime(2918, 6, 5, 0, 0, 0)
    # sat.to_continue_calculation(final_simulation_time)
    # sat.to_save('test2918', 'D:\\results')
    #
    # final_simulation_time = datetime(3018, 6, 5, 0, 0, 0)
    # sat.to_continue_calculation(final_simulation_time)
    # sat.to_save('test3018', 'D:\\results')

    #file = open('D:\\Documents\\Конференция Леса Евразии\\Test3.txt', 'w')
    #for count_inx in range(0, len(sat.interpreted_segments_of_polygons)):
    #    for pol_inx in range(0, len(sat.interpreted_segments_of_polygons[count_inx])):
    #        for seg_inx in sat.interpreted_segments_of_polygons[count_inx][pol_inx]:
    #            file.write(str(sat.area_of_interest.splitted_polygons_list[pol_inx].segments_geo_coordinates_list[seg_inx][0]) + '\t' +
    #                       str(sat.area_of_interest.splitted_polygons_list[pol_inx].segments_geo_coordinates_list[seg_inx][1]) + '\n')
    #file.close()

    print(datetime.now() - start)
