import numpy as np
import pandas as pd
import statistics
import os
from copy import deepcopy
from datetime import timedelta
from calendar import isleap
from matplotlib import pyplot as plt
from sklearn.neighbors.kde import KernelDensity
from Time.TimeUnits import to_get_unit_in_seconds, seconds_to_unit, unit_in_symbol
from Time.DaysNumber import to_determine_days_number_in_not_leap_year, to_determine_date_by_days_number_in_not_leap_year


DAYS_IN_NOT_LEAP_YEAR = 365
DAYS_IN_LEAP_YEAR = 366


class SatellitesStatistic:
    def __init__(self, satellites_list):
        self.satellites_list = satellites_list
        self.time_of_capturing = []
        for satellite in self.satellites_list:
            self.time_of_capturing.append(satellite.time_of_interpreted_capturing)
        self.time_of_capturing = np.sort(np.unique(np.concatenate(np.array(self.time_of_capturing))))

        self.captured_segments_of_polygons = []
        for time_count in self.time_of_capturing:
            self.captured_segments_of_polygons.append([])
            for satellite in self.satellites_list:
                sat_time_inx = np.where(satellite.time_of_interpreted_capturing == time_count)[0]
                if len(sat_time_inx) > 0:
                    for capt_pol in satellite.interpreted_segments_of_polygons[sat_time_inx[0]]:
                        self.captured_segments_of_polygons[-1].append(capt_pol)
        for time_count_inx in range(0, len(self.captured_segments_of_polygons)):
            for capt_pol_inx in range(0, len(self.captured_segments_of_polygons[time_count_inx])):
                self.captured_segments_of_polygons[time_count_inx][capt_pol_inx] = \
                    list(np.unique(np.array(
                        self.captured_segments_of_polygons[time_count_inx][capt_pol_inx])))
        self.area_of_interest = self.satellites_list[0].area_of_interest
        self.step = np.timedelta64(self.satellites_list[0].step, 's')
        self.initial_simulation_time = self.satellites_list[0].initial_simulation_time
        self.final_simulation_time = self.satellites_list[0].final_simulation_time
        self.initial_annual_observation_period = self.satellites_list[0].initial_annual_observation_period
        self.final_annual_observation_period = self.satellites_list[0].final_annual_observation_period
        self.capturing_counts = None
        self.spans_time_counts = None
        self.solution_time_counts = None
        self.spans_periods_sec = None
        self.solution_periods_sec = None

    def to_calculate_counts(self, min_percent_solution):
        self.spans_time_counts = []
        self.solution_time_counts = []
        self.capturing_counts = []
        for pol_inx in range(0, len(self.area_of_interest.splitted_polygons_list)):
            self.capturing_counts.append([])
            for seg_inx in range(0, len(self.area_of_interest.splitted_polygons_list[pol_inx].
                                                segments_geo_coordinates_list)):
                self.capturing_counts[pol_inx].append(0)
        self.spans_time_counts.append(self.time_of_capturing[0])
        self.to_add_captured_segments(0)
        if self.to_calc_percentages_of_captured_segments_n_times(1) >= min_percent_solution:
            self.solution_time_counts.append(self.time_of_capturing[0])
        for time_count_inx in range(1, len(self.captured_segments_of_polygons[1:])):
            last_count_delta = np.timedelta64(self.time_of_capturing[time_count_inx] -
                                              self.time_of_capturing[time_count_inx - 1], 's')
            if last_count_delta > self.step:
                self.spans_time_counts.append(self.time_of_capturing[time_count_inx])
            self.to_add_captured_segments(time_count_inx)
            if self.to_calc_percentages_of_captured_segments_n_times(len(self.solution_time_counts) + 1) >= min_percent_solution:
                self.solution_time_counts.append(self.time_of_capturing[time_count_inx])

    def to_add_captured_segments(self, time_count_inx):
        for pol_inx in range(0, len(self.captured_segments_of_polygons[time_count_inx])):
            for seg in self.captured_segments_of_polygons[time_count_inx][pol_inx]:
                self.capturing_counts[pol_inx][seg] += 1

    def to_calc_percentages_of_captured_segments_n_times(self, n):
        """
        @Описание:
            Метод определяет ход выполнения тематической задачи – какая площадь иили сколько процентов площади тестовых
                полигонов попало в поле зрения ГСК, сколько раз
        :param result_in_percents: если result_in_percents=True, то результат возвращается в процентах, если
            result_in_percents=False, то результат возвращается в м^2
        :return: список, в каждом элементе которого содержится процент площади группы полигонов, захваченных в ПЗ ГСК n
            раз, где n равняется номеру элемента списка плюс один
        """
        captured_area = 0
        for pol_inx in range(0, len(self.capturing_counts)):
            for captured_seg_inx, times_of_captured in enumerate(self.capturing_counts[pol_inx]):
                if times_of_captured >= n:
                    captured_area += self.area_of_interest.splitted_polygons_list[pol_inx].segments_area_list[
                        captured_seg_inx]
        return captured_area / self.area_of_interest.full_area * 100

    def to_calculate_counts_data_obsolescence(self, min_percent_solution, days_to_data_aging):
        self.spans_time_counts = []
        self.solution_time_counts = []
        self.capturing_counts_sample = []
        self.capturing_counts_list = []
        self.time_of_capturing_list = []
        for pol_inx in range(0, len(self.area_of_interest.splitted_polygons_list)):
            self.capturing_counts_sample.append([])
            for seg_inx in range(0, len(self.area_of_interest.splitted_polygons_list[pol_inx].
                                                segments_geo_coordinates_list)):
                self.capturing_counts_sample[pol_inx].append(0)
        self.spans_time_counts.append(self.time_of_capturing[0])
        self.to_add_captured_segments_data_obsolescence(0)
        if self.to_calc_percentages_of_captured_segments_n_times_data_obsolescence(1) >= min_percent_solution:
            self.solution_time_counts.append(self.time_of_capturing[0])
        for time_count_inx in range(1, len(self.captured_segments_of_polygons[1:])):
            time_of_aging = self.time_of_capturing[time_count_inx] - np.timedelta64(days_to_data_aging, 'D')
            aging_counts = np.flip(np.where(np.array(self.time_of_capturing_list) < time_of_aging)[0])
            for i in aging_counts:
                del self.time_of_capturing_list[i]
                del self.capturing_counts_list[i]
            last_count_delta = np.timedelta64(self.time_of_capturing[time_count_inx] -
                                              self.time_of_capturing[time_count_inx - 1], 's')
            if last_count_delta > self.step:
                self.spans_time_counts.append(self.time_of_capturing[time_count_inx])
            self.to_add_captured_segments_data_obsolescence(time_count_inx)
            if self.to_calc_percentages_of_captured_segments_n_times_data_obsolescence(len(self.solution_time_counts) + 1) >= min_percent_solution:
                self.solution_time_counts.append(self.time_of_capturing[time_count_inx])

    def to_add_captured_segments_data_obsolescence(self, time_count_inx):
        for pol_inx in range(0, len(self.captured_segments_of_polygons[time_count_inx])):
            current_time = self.time_of_capturing[time_count_inx]
            self.capturing_counts_list.append(deepcopy(self.capturing_counts_sample))
            self.time_of_capturing_list.append(current_time)
            for seg in self.captured_segments_of_polygons[time_count_inx][pol_inx]:
                self.capturing_counts_list[-1][pol_inx][seg] += 1

    def to_calc_percentages_of_captured_segments_n_times_data_obsolescence(self, n, data_aging_years=None):
        """
        @Описание:
            Метод определяет ход выполнения тематической задачи – какая площадь иили сколько процентов площади тестовых
                полигонов попало в поле зрения ГСК, сколько раз
        :param result_in_percents: если result_in_percents=True, то результат возвращается в процентах, если
            result_in_percents=False, то результат возвращается в м^2
        :return: список, в каждом элементе которого содержится процент площади группы полигонов, захваченных в ПЗ ГСК n
            раз, где n равняется номеру элемента списка плюс один
        """
        capturing_counts = deepcopy(self.capturing_counts_sample)
        for count in range(0, len(self.capturing_counts_list)):
            for pol_inx in range(0, len(self.capturing_counts_list[count])):
                for seg_inx in range(0, len(self.capturing_counts_list[count][pol_inx])):
                    if self.capturing_counts_list[count][pol_inx][seg_inx] > 0:
                        capturing_counts[pol_inx][seg_inx] += 1
        captured_area = 0
        for pol_inx in range(0, len(self.capturing_counts_sample)):
            for captured_seg_inx, times_of_captured in enumerate(capturing_counts[pol_inx]):
                if times_of_captured >= n:
                    captured_area += self.area_of_interest.splitted_polygons_list[pol_inx].segments_area_list[
                        captured_seg_inx]
        return captured_area / self.area_of_interest.full_area * 100

    def to_determine_periods(self, to_skip_time_out_of_observation_period=False):
        self.spans_periods_sec = self.to_calculate_periods_sec(self.spans_time_counts,
                                                               self.initial_simulation_time,
                                                               to_skip_time_out_of_observation_period)
        self.solution_periods_sec = self.to_calculate_periods_sec(self.solution_time_counts,
                                                                  self.initial_simulation_time,
                                                                  to_skip_time_out_of_observation_period)
        # if len(self.solution_time_counts) > 0:
        #     self.solution_periods_sec = self.to_calculate_periods_sec(self.solution_time_counts,
        #                                                               self.initial_simulation_time,
        #                                                               to_skip_time_out_of_observation_period)
        # else:
        #     self.solution_periods_sec = self.to_calculate_periods_sec(self.solution_time_counts,
        #                                                               None,
        #                                                               to_skip_time_out_of_observation_period)

    def to_calculate_periods_sec(self, time_list, first_time=None, to_skip_time_out_of_observation_period=False):
        """
            @Описание:
                Метод определяет время от времени начала наблюдения до первого значения времени из заданного списка и между
                    соседними значениями времени из заданного списка в секундлах. Возможно вычисление периодов только по
                    допустимому годовому периоду наблюдения.
            :param time_list: список значений времени, по которому будут вычисляться (date_time).
            :param first_time: начальное время наблюдений, от которого будет отсчитываться первый период (datetime).
                Допустимо None. При этом сразу вычисляются периоды между значениями time_list. По умолчанию None.
            :param initial_annual_observation_period: номер первого дня в невисокосном году допустимого годового периода
                наблюдения (времени в году, когда допустима съемка) (int). Задается методом to_set_annual_observations_period.
                По умолчанию DEFAULT_INITIAL_ANNUAL_OBSERVATION_PERIOD.
            :param final_annual_observation_period: номер последнего дня в невисокосном году годового подустимого периода
                наблюдений (времени в году, когда допустима съемка) (int). Задается методом to_set_annual_observations_period.
                По умолчанию DEFAULT_FINAL_ANNUAL_OBSERVATION_PERIOD.
            :param to_skip_time_out_of_observation_period: пропускать промежутки времени, не входящие в годовые периоды
                   наблюдения (boolean). По умолчанию False.
            :return: список приодов в секундах (int, double)
            """
        # Если нету пары значений времени (начального времени first_time и одного значения из time_list или два значения
        #   из time_list), то возвращается [], а дальнейшие вычисления не производятся. Если есть, то задаётся первое
        #   значение времени, от которого будет вестись отсчет - previous_time.
        if first_time is not None and len(time_list) > 0 :
            time_data_frame = np.append(np.datetime64(first_time), time_list)
        else:
            time_data_frame = time_list
        if len(time_data_frame) > 1:
            previous_time = np.datetime64(time_data_frame[0])

        else:
            return []
        # В список periods записываются периоды
        periods = []
        time_array = np.array(time_data_frame)
        # Если не пропускается время между периодами наблюдения
        if not to_skip_time_out_of_observation_period:
            # Список periods заполняется значениями периодов
            for time in time_array[1:]:
                period = np.timedelta64(time - previous_time, 's').astype(int)
                periods.append(period)
                previous_time = time
        # Если пропускается время между периодами наблюдения
        else:
            # Определение индикатора observation_period_inside_one_year, который показывает, что
            #   final_annual_observation_period больше initial_annual_observation_period. Это означает, что что и начало и
            #   конец периода наблюдения находятся внутри одного года без перехода в следующий (если True, если False, то
            #   наоборот).
            observation_period_inside_one_year = self.final_annual_observation_period > \
                                                 self.initial_annual_observation_period
            # Вычисление даты начала годового периода наблюдения
            years_list = [int(str(t.astype('datetime64[Y]'))) for t in time_data_frame]
            current_year = years_list[0]
            start_of_current_observation_period = to_determine_date_by_days_number_in_not_leap_year(
                self.initial_annual_observation_period, current_year)
            end_of_current_observation_period = to_determine_date_by_days_number_in_not_leap_year(
                self.final_annual_observation_period, current_year)
            # Если previous_time оказывается не в периоде наблюдения, оно перемещается на начало следующего периода
            if observation_period_inside_one_year:
                if start_of_current_observation_period > previous_time > end_of_current_observation_period:
                    previous_time = start_of_current_observation_period
            else:
                if previous_time < start_of_current_observation_period:
                    previous_time = start_of_current_observation_period
                    if previous_time < end_of_current_observation_period:
                        if not isleap(current_year):
                            days_in_current_year = DAYS_IN_LEAP_YEAR
                        else:
                            days_in_current_year = DAYS_IN_NOT_LEAP_YEAR
                        previous_time += timedelta(days=days_in_current_year)
            # Если период наблюдений не пересекает границу соседних годов
            if observation_period_inside_one_year:
                for time_inx in range(1, len(time_array)):
                    current_year = years_list[time_inx]
                    year_of_previous_time = years_list[time_inx - 1]
                    skipped_time = 0
                    year_of_skipped_period = current_year
                    # Вычисление времени, не входящего в период наблюдения
                    while year_of_skipped_period > year_of_previous_time:
                        start_of_new_observation_period = to_determine_date_by_days_number_in_not_leap_year(
                            self.initial_annual_observation_period, year_of_skipped_period)
                        end_of_old_observation_period = to_determine_date_by_days_number_in_not_leap_year(
                            self.final_annual_observation_period, year_of_skipped_period - 1)
                        skipped_time += (
                                start_of_new_observation_period - end_of_old_observation_period).total_seconds()
                        year_of_skipped_period -= 1
                    # Вычисление периода без времени, не входящего в период наблюдения
                    period = np.timedelta64(time_array[time_inx] - previous_time, 's').astype(int) - skipped_time
                    periods.append(period)
                    previous_time = time_array[time_inx]
            # Если период наблюдений пересекает границу соседних годов
            else:
                for time_inx in range(1, len(time_array)):
                    day_of_time_in_year = to_determine_days_number_in_not_leap_year(previous_time)
                    if day_of_time_in_year >= self.initial_annual_observation_period:
                        year_of_not_observing_period = previous_time.year + 1
                    else:
                        year_of_not_observing_period = previous_time.year
                    end_of_current_observation_period = to_determine_date_by_days_number_in_not_leap_year(
                        self.final_annual_observation_period, year_of_not_observing_period)
                    skipped_time = 0
                    time_copy = time_array[time_inx]
                    # Вычисление времени, не входящего в период наблюдения
                    while time_copy >= end_of_current_observation_period:
                        day_of_time_copy_in_year = to_determine_days_number_in_not_leap_year(previous_time)
                        if day_of_time_copy_in_year >= self.initial_annual_observation_period:
                            year_of_skipped_period = previous_time.year
                        else:
                            year_of_skipped_period = previous_time.year - 1
                        start_of_new_observation_period = to_determine_date_by_days_number_in_not_leap_year(
                            self.initial_annual_observation_period, year_of_skipped_period)
                        end_of_old_observation_period = to_determine_date_by_days_number_in_not_leap_year(
                            self.final_annual_observation_period, year_of_skipped_period)
                        skipped_time += (start_of_new_observation_period - end_of_old_observation_period). \
                            total_seconds()
                        time_copy = to_determine_date_by_days_number_in_not_leap_year(
                            to_determine_days_number_in_not_leap_year(time_copy), time_copy.year - 1)
                    # Вычисление периода без времени, не входящего в период наблюдения
                    period = np.timedelta64(time_array[time_inx] - previous_time, 's').astype(int) - skipped_time
                    periods.append(period)
                    previous_time = time_array[time_inx]
        return periods

    def to_calculate_data(self, file_name, directory, data_time_unit='days', time_of_simulation_unit='years'):
        average_spans_list = []
        median_spans_list = []
        max_spans_list = []
        min_spans_list = []
        dispersion_of_spans_list = []
        standard_deviation_of_spans_list = []
        spans_times_list = []
        average_solutions_list = []
        median_solutions_list = []
        max_solutions_list = []
        min_solutions_list = []
        dispersion_of_solutions_list = []
        standard_deviation_of_solutions_list = []
        solutions_times_list = []
        for time_count_inx, time_count in enumerate(self.time_of_capturing):
            current_spans_periods = np.array(self.spans_periods_sec)[
                np.where(self.spans_time_counts <= time_count)]
            current_average_spans, \
            current_median_spans, \
            current_max_spans, \
            current_min_spans, \
            current_dispersion_of_spans, \
            current_standard_deviation_of_spans = \
                to_determine_median_average_max_min_dispersion_standard(current_spans_periods, data_time_unit)
            average_spans_list.append(current_average_spans)
            median_spans_list.append(current_median_spans)
            max_spans_list.append(current_max_spans)
            min_spans_list.append(current_min_spans)
            dispersion_of_spans_list.append(current_dispersion_of_spans)
            standard_deviation_of_spans_list.append(current_standard_deviation_of_spans)
            spans_times_list.append(len(current_spans_periods))
            current_solutions_periods = np.array(self.solution_periods_sec)[
                np.where(self.solution_time_counts <= time_count)]
            current_average_solutions, \
            current_median_solutions, \
            current_max_solutions, \
            current_min_solutions, \
            current_dispersion_of_solutions, \
            current_standard_deviation_of_solutions = \
                to_determine_median_average_max_min_dispersion_standard(current_solutions_periods, data_time_unit)
            average_solutions_list.append(current_average_solutions)
            median_solutions_list.append(current_median_solutions)
            max_solutions_list.append(current_max_solutions)
            min_solutions_list.append(current_min_solutions)
            dispersion_of_solutions_list.append(current_dispersion_of_solutions)
            standard_deviation_of_solutions_list.append(current_standard_deviation_of_solutions)
            solutions_times_list.append(len(current_solutions_periods))
        data_time_unit_symbol = unit_in_symbol(data_time_unit)
        time_of_simulation_unit_symbol = unit_in_symbol(time_of_simulation_unit)
        data_header = pd.MultiIndex.from_product(
            np.array([['Spans', 'Solutions'],
                      ['Average (' + data_time_unit_symbol  + ')',
                       'Median (' + data_time_unit_symbol  + ')',
                       'Max (' + data_time_unit_symbol  + ')',
                       'Min (' + data_time_unit_symbol  + ')',
                       'Dispersion (' + data_time_unit_symbol  + ')',
                       'Standard deviation (' + data_time_unit_symbol  + ')',
                       'Number']]))
        table = pd.DataFrame(np.array([np.array(average_spans_list),
                                       np.array(median_spans_list),
                                       np.array(max_spans_list),
                                       np.array(min_spans_list),
                                       np.array(dispersion_of_spans_list),
                                       np.array(standard_deviation_of_spans_list),
                                       np.array(spans_times_list),
                                       np.array(average_solutions_list),
                                       np.array(median_solutions_list),
                                       np.array(max_solutions_list),
                                       np.array(min_solutions_list),
                                       np.array(dispersion_of_solutions_list),
                                       np.array(standard_deviation_of_solutions_list),
                                       np.array(solutions_times_list)]).T, columns=data_header)
        time_header = pd.MultiIndex.from_product(np.array([['Time'], ['Date and time', 'Time of simulation (' +
                                                                      time_of_simulation_unit_symbol  + ')']]))
        time_of_simulation = (self.time_of_capturing -
                              np.datetime64(self.initial_simulation_time)).astype('timedelta64[s]').astype('int') / \
                             to_get_unit_in_seconds(time_of_simulation_unit)
        time_table = pd.DataFrame(np.array([self.time_of_capturing, time_of_simulation]).T, columns=time_header)
        table = pd.concat([time_table, table], axis=1)
        if not os.path.exists(directory):
            os.makedirs(directory)
        # сохранение таблицы
        table.to_excel(directory + '\\' + file_name + '.xlsx')

    def to_draw_histogram_of_spans(self, title='Histogram of spans periods', unit_of_time='days',
                                   x_size=10, y_size=6, dpi=100, title_font=16, ticks_label_font=12, axis_title_font=14,
                                   color_ind='blue', show=True, file_name=None, directory_address=None):
        to_draw_histogram(self.spans_periods_sec, title,
                          'Spans periods (' + unit_in_symbol(unit_of_time) + ')', 'Spans number', unit_of_time,
                          x_size, y_size, dpi, title_font, ticks_label_font, axis_title_font,
                          color_ind, show, file_name, directory_address)

    def to_draw_histogram_of_solutions(self, title='Histogram of solutions periods', unit_of_time='days',
                          x_size=10, y_size=6, dpi=100, title_font=16, ticks_label_font=12, axis_title_font=14,
                          color_ind='blue', show=True, file_name=None, directory_address=None):
        to_draw_histogram(self.solution_periods_sec, title,
                          'Spans periods (' + unit_in_symbol(unit_of_time) + ')', 'Solutions number', unit_of_time,
                          x_size, y_size, dpi, title_font, ticks_label_font, axis_title_font,
                          color_ind, show, file_name, directory_address)

    def to_draw_distribution_density_of_spans(self, title='Distribution density of spans period', unit_of_time='days',
                          x_size=10, y_size=6, dpi=100, title_font=16, ticks_label_font=12, axis_title_font=14,
                          color_ind='blue', show=True, file_name=None, directory_address=None):
        to_draw_distribution_density(self.spans_periods_sec, title,
                                     'Spans periods (' + unit_in_symbol(unit_of_time) + ')',
                                     'Probability density of period', unit_of_time,
                                     x_size, y_size, dpi, title_font, ticks_label_font, axis_title_font,
                                     color_ind, show, file_name, directory_address)

    def to_draw_distribution_density_of_spans_by_kde(self, freq, bandwidth, kernel='gaussian',
                                                     title='Distribution density of spans period', unit_of_time='days',
                                                     x_size=10, y_size=6, dpi=100, title_font=16, ticks_label_font=12,
                                                     axis_title_font=14, color_ind='blue', show=True, file_name=None,
                                                     directory_address=None):
        to_draw_distribution_density_by_kde(self.spans_periods_sec, freq, bandwidth, title,
                                            'Spans periods (' + unit_in_symbol(unit_of_time) + ')',
                                            'Probability density of period', kernel, unit_of_time,
                                            x_size, y_size, dpi, title_font, ticks_label_font, axis_title_font,
                                            color_ind, show, file_name, directory_address)

    def to_draw_distribution_density_of_solutions(self, title='Distribution density of solutions period', unit_of_time='days',
                          x_size=10, y_size=6, dpi=100, title_font=16, ticks_label_font=12, axis_title_font=14,
                          color_ind='blue', show=True, file_name=None, directory_address=None):
        to_draw_distribution_density(self.solution_periods_sec, title,
                                     'Spans periods (' + unit_in_symbol(unit_of_time) + ')',
                                     'Probability density of period', unit_of_time,
                                     x_size, y_size, dpi, title_font, ticks_label_font, axis_title_font,
                                     color_ind, show, file_name, directory_address)

    def to_draw_distribution_density_of_solutions_by_kde(self, freq, bandwidth, kernel='gaussian',
                                                         title='Distribution density of spans period', unit_of_time='days',
                                                         x_size=10, y_size=6, dpi=100, title_font=16, ticks_label_font=12,
                                                         axis_title_font=14, color_ind='blue', show=True, file_name=None,
                                                         directory_address=None):
        to_draw_distribution_density_by_kde(self.solution_periods_sec, freq, bandwidth, title,
                                            'Spans periods (' + unit_in_symbol(unit_of_time) + ')',
                                            'Probability density of period', kernel, unit_of_time,
                                            x_size, y_size, dpi, title_font, ticks_label_font, axis_title_font,
                                            color_ind, show, file_name, directory_address)

    def to_draw_probability_of_spans(self, title='Probability of spans period', unit_of_time='days',
                          x_size=10, y_size=6, dpi=100, title_font=16, ticks_label_font=12, axis_title_font=14,
                          color_ind='blue', show=True, file_name=None, directory_address=None):
        to_draw_probability(self.spans_periods_sec, title,
                            'Spans periods (' + unit_in_symbol(unit_of_time) + ')',
                            'Probability of period', unit_of_time,
                            x_size, y_size, dpi, title_font, ticks_label_font, axis_title_font,
                            color_ind, show, file_name, directory_address)

    def to_draw_probability_of_solutions(self, title='Probability of solutions period', unit_of_time='days',
                                         x_size=10, y_size=6, dpi=100, title_font=16, ticks_label_font=12,
                                         axis_title_font=14, color_ind='blue', show=True, file_name=None,
                                         directory_address=None):
        to_draw_probability(self.solution_periods_sec, title,
                            'Solutions periods (' + unit_in_symbol(unit_of_time) + ')',
                            'Probability of period', unit_of_time,
                            x_size, y_size, dpi, title_font, ticks_label_font, axis_title_font,
                            color_ind, show, file_name, directory_address)


def to_determine_median_average_max_min_dispersion_standard(numbers, time_unit='days'):
    """
    @Описание:
        По списку чисел вычисляет их среднее, медианное, максимальное, минимальное значение, их дисперсию и
            среднеквадратическое отклонение.
    :param numbers: список чисел, для которых проводится вычисления.
    :return:
        average_value: среднее значение по списку numbers (double)
        median_value: медианное значение по списку numbers (double)
        max_value: максимальное значение по списку numbers (double)
        min_value: минимальное значение по списку numbers (double)
        dispersion: дисперсия по списку numbers (double)
        standard_deviation: среднеквадратическое отклонение по списку numbers (double)
    """
    if len(numbers) > 0:
        unit_in_seconds = to_get_unit_in_seconds(time_unit)
        average_value = statistics.mean(numbers) / unit_in_seconds
        median_value = statistics.median(numbers) / unit_in_seconds
        max_value = max(numbers) / unit_in_seconds
        min_value = min(numbers) / unit_in_seconds
        if len(numbers) > 1:
            dispersion = statistics.variance(numbers) / unit_in_seconds
            standard_deviation = statistics.stdev(numbers) / unit_in_seconds
        else:
            dispersion = 0
            standard_deviation = 0
    else:
        average_value = 0
        median_value = 0
        max_value = 0
        min_value = 0
        dispersion = 0
        standard_deviation = 0
    return average_value, median_value, max_value, min_value, dispersion, standard_deviation


def to_draw_histogram(values, title, x_label, y_label, unit_of_time='days',
                      x_size=10, y_size=6, dpi=100, title_font=16, ticks_label_font=12, axis_title_font=14,
                      color_ind='blue', show=True, file_name=None, directory_address=None):
    histogram = to_make_histogram(values, unit_of_time)
    np.sum(histogram * np.arange(0, len(histogram)))
    to_draw_graph(histogram, title, x_label, y_label, x_size, y_size, dpi, title_font, ticks_label_font,
                  axis_title_font, color_ind, show, file_name, directory_address)


def to_draw_distribution_density(values, title, x_label, y_label, unit_of_time='days',
                                 x_size=10, y_size=6, dpi=100, title_font=16, ticks_label_font=12, axis_title_font=14,
                                 color_ind='blue', show=True, file_name=None, directory_address=None):
    histogram = to_make_histogram(values, unit_of_time)
    distribution_density = histogram / sum(histogram)
    to_draw_graph(distribution_density, title, x_label, y_label, x_size, y_size, dpi, title_font, ticks_label_font,
                  axis_title_font, color_ind, show, file_name, directory_address)


def to_draw_probability(values, title, x_label, y_label, unit_of_time='days',
                        x_size=10, y_size=6, dpi=100, title_font=16, ticks_label_font=12, axis_title_font=14,
                        color_ind='blue', show=True, file_name=None, directory_address=None):
    histogram = to_make_histogram(values, unit_of_time)
    distribution_density = histogram / sum(histogram)
    probability = np.zeros(len(distribution_density))
    probability[0] = distribution_density[0]
    for i in range(1, len(probability[1:])):
        probability[i] = distribution_density[i] + probability[i - 1]
    to_draw_graph(probability, title, x_label, y_label, x_size, y_size, dpi, title_font, ticks_label_font,
                  axis_title_font, color_ind, show, file_name, directory_address)


def to_draw_distribution_density_by_kde(values, freq, bandwidth,  title, x_label, y_label, kernel='gaussian',
                                        unit_of_time='days', x_size=10, y_size=6, dpi=100, title_font=16,
                                        ticks_label_font=12, axis_title_font=14, color_ind='blue', show=True,
                                        file_name=None, directory_address=None):
    step = to_get_unit_in_seconds(unit_of_time)
    # Определение максимума гистограммы по времени в секундах
    max_value = seconds_to_unit((max(values) // step) * step + step, unit_of_time)
    value_grid = np.linspace(0, max_value, freq * max_value)[:, np.newaxis]
    kde = KernelDensity(kernel=kernel, bandwidth=bandwidth).fit(seconds_to_unit(np.array(values), unit_of_time)
                                                                [:, np.newaxis])
    dens = np.exp(kde.score_samples(value_grid))
    to_draw_compressed_graph(dens, freq, title, x_label, y_label, x_size, y_size, dpi, title_font,
                             ticks_label_font, axis_title_font, color_ind, show, file_name, directory_address)

def to_make_histogram(values, unit_of_time='days'):
    # Перевод единиц измерения unit_of_time в секунды
    step = to_get_unit_in_seconds(unit_of_time)
    # Определение максимума гистограммы по времени в секундах
    max_value = (max(values) // step) * step + step
    # Перевод максимума гистограммы по времени из секунд в unit_of_time
    histogram_len = int(max_value // step)
    # Составление гистограммы
    histogram = np.zeros((histogram_len))
    for value in values:
        index = int(value // step)
        histogram[index] += 1
    return histogram


def to_draw_graph(values, title, x_label, y_label, x_size=10, y_size=6, dpi=100, title_font=16, ticks_label_font=12,
                  axis_title_font=14, color_ind='blue', show=True, file_name=None, directory_address=None):
    """
            @Описание:
                Метод принимает данные о гистограмме и строит по ним графическое представление с помощью пакета matplotlib.
                    Также для метода, задаются параметры графического представления построенной гистограммы: высота, ширина,
                    количество пикселей на дюйм, размеры шрифтов подписей обозначений на осях, заголовка и осей, а также
                    сами заголовок и названия осей. После постраения гистограмма сохраняется по выбранному адресу в png, pdf
                    или обоих форматах.
            :param hist_address: адрес по которому сохраняется гистограмма (с названием, но без формата) (String)
            :param hist_title: название гистограммы, которым он будет подписан в графическом представлении (String)
            :param hist_values: данные о гистограмме (значения каждого столбца в списке) (int)
            :param x_size: размер графического представления графика по оси x (в ширину) в дюймах (int, float)
            :param y_size: размер графического представления графика по оси y (в высоту) в дюймах (int, float)
            :param dpi: количество пикселей на дюм в графическом представлении графика (int)
            :param x_label: название оси x, которым будет подписываться эта ось на графическом представлении (String)
            :param y_label: название оси y, которым будет подписываться эта ось на графическом представлении (String)
            :param ticks_label_font: размер шрифта названия графика graph_title, которым он подписывается на графическом
                представлении (int)
            :param title_font: размер шрифта названия графика осей x_label, y_label на графическом представлении (int)
            :param axis_title_font: размер шрифта, обозначений засечек на осях x и y в графическом представлении (int)
            :param color_ind: текстовое обозначение цвета графика (String)
            :param output_in_png: если True, то выводить график в формате png (boolean). По умолчанию False
            :param output_in_pdf: если True, то выводить график в формате pdf (boolean). По умолчанию True
            :return: сохраняет по адресу hist_address графическое представление гистограммы в форматах png, pdf или обоих
                одновременно, построенный по зависимости hist_title, time_axis_original с размерами x_size, y_size,
                разрешениеь dpi на дюйм, с подписью графика graph_title, подписями осей x_label, y_label, с рамерами
                шрифтов названия, подписей осей и подписей осей ticks_label_font, title_font, axis_title_font,
                соответственно, цветом графика color.
            """
    graph = plt.figure(figsize=(x_size, y_size), dpi=dpi)
    ax = graph.add_subplot(111)

    plt.bar(range(0, len(values)), values, width=1, color=color_ind)
    # Изменение размера шрифта подписей оси x
    for label in ax.xaxis.get_ticklabels():
        label.set_fontsize(ticks_label_font)
    # Изменение размера шрифта подписей оси y
    for label in ax.yaxis.get_ticklabels():
        label.set_fontsize(ticks_label_font)
    plt.title(title, fontsize=title_font)
    plt.ylabel(y_label, fontsize=axis_title_font)
    plt.xlabel(x_label, fontsize=axis_title_font)
    plt.grid(True)
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['font.family'] = 'Calibri'
    # Позиции засечек на графике минимальных, чтобы весь график был виден для обеих осей
    y_ticks = ax.get_yticks()
    max_y_value = max(values)
    i = -1
    while ~(y_ticks[i] >= max_y_value > y_ticks[i - 1]):
        i = i - 1
    min_suitable_y_tick = y_ticks[i]
    x_ticks = ax.get_xticks()
    max_x_value = len(values)
    i = -1
    while ~(x_ticks[i] >= max_x_value > x_ticks[i - 1]):
        i = i - 1
    min_suitable_x_tick = x_ticks[i]
    plt.axis([0, min_suitable_x_tick, 0, min_suitable_y_tick])
    if show:
        plt.show()
    # Вывод зависимости в виде графика в файле png, если требуется
    if file_name is not None and directory_address is not None:
        graph.savefig("".join([directory_address, '\\', file_name, '.png']))


def to_draw_compressed_graph(values, comp_coef, title, x_label, y_label, x_size=10, y_size=6, dpi=100, title_font=16,
                             ticks_label_font=12, axis_title_font=14, color_ind='blue', show=True, file_name=None,
                             directory_address=None):
    """
            @Описание:
                Метод принимает данные о гистограмме и строит по ним графическое представление с помощью пакета matplotlib.
                    Также для метода, задаются параметры графического представления построенной гистограммы: высота, ширина,
                    количество пикселей на дюйм, размеры шрифтов подписей обозначений на осях, заголовка и осей, а также
                    сами заголовок и названия осей. После постраения гистограмма сохраняется по выбранному адресу в png, pdf
                    или обоих форматах.
            :param hist_address: адрес по которому сохраняется гистограмма (с названием, но без формата) (String)
            :param hist_title: название гистограммы, которым он будет подписан в графическом представлении (String)
            :param hist_values: данные о гистограмме (значения каждого столбца в списке) (int)
            :param x_size: размер графического представления графика по оси x (в ширину) в дюймах (int, float)
            :param y_size: размер графического представления графика по оси y (в высоту) в дюймах (int, float)
            :param dpi: количество пикселей на дюм в графическом представлении графика (int)
            :param x_label: название оси x, которым будет подписываться эта ось на графическом представлении (String)
            :param y_label: название оси y, которым будет подписываться эта ось на графическом представлении (String)
            :param ticks_label_font: размер шрифта названия графика graph_title, которым он подписывается на графическом
                представлении (int)
            :param title_font: размер шрифта названия графика осей x_label, y_label на графическом представлении (int)
            :param axis_title_font: размер шрифта, обозначений засечек на осях x и y в графическом представлении (int)
            :param color_ind: текстовое обозначение цвета графика (String)
            :param output_in_png: если True, то выводить график в формате png (boolean). По умолчанию False
            :param output_in_pdf: если True, то выводить график в формате pdf (boolean). По умолчанию True
            :return: сохраняет по адресу hist_address графическое представление гистограммы в форматах png, pdf или обоих
                одновременно, построенный по зависимости hist_title, time_axis_original с размерами x_size, y_size,
                разрешениеь dpi на дюйм, с подписью графика graph_title, подписями осей x_label, y_label, с рамерами
                шрифтов названия, подписей осей и подписей осей ticks_label_font, title_font, axis_title_font,
                соответственно, цветом графика color.
            """
    graph = plt.figure(figsize=(x_size, y_size), dpi=dpi)
    ax = graph.add_subplot(111)

    plt.bar(np.arange(0, len(values)) / comp_coef, values, width=1 / comp_coef, color=color_ind)
    # Изменение размера шрифта подписей оси x
    for label in ax.xaxis.get_ticklabels():
        label.set_fontsize(ticks_label_font)
    # Изменение размера шрифта подписей оси y
    for label in ax.yaxis.get_ticklabels():
        label.set_fontsize(ticks_label_font)
    plt.title(title, fontsize=title_font)
    plt.ylabel(y_label, fontsize=axis_title_font)
    plt.xlabel(x_label, fontsize=axis_title_font)
    plt.grid(True)
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['font.family'] = 'Calibri'
    # Позиции засечек на графике минимальных, чтобы весь график был виден для обеих осей
    y_ticks = ax.get_yticks()
    max_y_value = max(values)
    i = -1
    while ~(y_ticks[i] >= max_y_value > y_ticks[i - 1]):
        i = i - 1
    min_suitable_y_tick = y_ticks[i]
    x_ticks = ax.get_xticks()
    max_x_value = len(values) / comp_coef
    i = -1
    while ~(x_ticks[i] >= max_x_value > x_ticks[i - 1]):
        i = i - 1
    min_suitable_x_tick = x_ticks[i]
    plt.axis([0, min_suitable_x_tick, 0, min_suitable_y_tick])
    if show:
        plt.show()
    # Вывод зависимости в виде графика в файле png, если требуется
    if file_name is not None and directory_address is not None:
        graph.savefig("".join([directory_address, '\\', file_name, '.png']))