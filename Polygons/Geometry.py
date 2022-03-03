import numpy as np


def polygon_contains_points(polygons_border_points, points):
    with np.errstate(all='ignore'):
        points_coordinates = np.swapaxes(points, 0, 1)
        x = points_coordinates[0]
        y = points_coordinates[1]
        borders_coordinates = np.swapaxes(polygons_border_points, 0, 1)
        polygons_x = borders_coordinates[0]
        polygons_y = borders_coordinates[1]

        x_table = np.swapaxes([x] * len(polygons_x), 0, 1)
        y_table = np.swapaxes([y] * len(polygons_y), 0, 1)
        polygons_x_table = [polygons_x] * len(x)
        polygons_y_table = [polygons_y] * len(y)
        shifted_polygons_x_table = np.roll(polygons_x_table, 1)
        shifted_polygons_y_table = np.roll(polygons_y_table, 1)

        logic_table = ((np.less_equal(polygons_y_table, y_table) * np.less(y_table, shifted_polygons_y_table)) +
                       (np.less_equal(shifted_polygons_y_table, y_table) * np.less(y_table, polygons_y_table))) * \
                      (np.greater(x_table, (shifted_polygons_x_table - polygons_x_table) *
                                  (y_table - polygons_y_table) / (shifted_polygons_y_table - polygons_y_table) +
                                  polygons_x_table))

        true_results_for_points = np.where(logic_table == True)[0]
        counts_of_true_results_for_points = np.array([len(np.where(true_results_for_points == i)[0])
                                                      for i in range(0, len(x))])
    return counts_of_true_results_for_points % 2 > 0
