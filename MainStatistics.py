import numpy as np
import scipy
from Satellite import Satellite
from Statistics import Statistics
from Statistics.Statistics import SatellitesStatistic
from Weather import Cloudiness

def bootstrap(stat, percent, num):
    sel = np.array(stat.solution_time_counts)
    new_num = []
    new_av = []
    new_med = []
    new_max = []
    new_min = []
    new_std = []
    for i in range(num):
        new_sel = sel[np.where(np.random.random_sample(len(sel)) > percent / 100)]
        new_diff = []
        for j in range(1, len(new_sel)):
            d = float((new_sel[j] - new_sel[j - 1]) / np.timedelta64(1, 's'))
            while d / 24 / 60 / 60 >= 365 - 243 + 120:
                d -= (365 - 243 + 120) * 24 * 60 * 60
            new_diff.append(d)
        new_num.append(len(new_diff) + 1)
        new_av.append(np.mean(np.array(new_diff)))
        new_med.append(np.median(np.array(new_diff)))
        #new_max.append(np.max(np.array(new_diff)))
        #new_min.append(np.min(np.array(new_diff)))
        new_std.append(np.std(np.array(new_diff)))
    return [
            np.mean(np.array(new_num)),
            np.mean(np.array(new_av)) / 24 / 60 / 60,
            np.mean(np.array(new_med)) / 24 / 60 / 60,
            #np.mean(np.array(new_max)) / 24 / 60 / 60,
            #np.mean(np.array(new_min)) / 24 / 60 / 60,
            np.mean(np.array(new_std)) / 24 / 60 / 60]

if __name__ == '__main__':
    # sentinel_2a = Satellite.to_load_satellite('D:/Проекты/Оперативность/Бронницкое лесничество/Промежуточные данные/sentinel_2a_2031_s.file')
    # sentinel_2b = Satellite.to_load_satellite('D:/Проекты/Оперативность/Бронницкое лесничество/Промежуточные данные/sentinel_2b_2031_s.file')
    # sentinel_2a.to_consider_solar_angle(40)
    # sentinel_2b.to_consider_solar_angle(40)
    resurs_p3 = Satellite.to_load_satellite('D:/Проекты/Оперативность/Бронницкое лесничество/Промежуточные данные/resurs_p3_hs_2121_s.file')
    resurs_p3.to_consider_solar_angle(40)
    # cloud_distr = Cloudiness.to_load_cloud_distr_dataset('D:/Проекты/Оперативность/Бронницкое лесничество/Промежуточные данные/'
    #                                                      'Бронницкое_лесничество_cloud_cover.file')
    # sentinel_2a.to_consider_cloudiness(1, cloud_distr)
    # sentinel_2b.to_consider_cloudiness(1, cloud_distr)
    #resurs_p3.to_consider_cloudiness(1, cloud_distr)
    # sat_stat = SatellitesStatistic([sentinel_2a, sentinel_2b])
    sat_stat = SatellitesStatistic([resurs_p3])
    sat_stat.to_calculate_counts(99)
    #sat_stat.to_calculate_counts_data_obsolescence(100, 1)
    sat_stat.to_determine_periods(to_skip_time_out_of_observation_period=True)
    #sat_stat.to_draw_distribution_density_of_spans()
    #sat_stat.to_draw_distribution_density_of_spans_by_kde(100, 0.5)
    #sat_stat.to_draw_probability_of_spans()
    # sat_stat.to_draw_distribution_density_of_solutions()
    # sat_stat.to_draw_probability_of_solutions()
    print(bootstrap(sat_stat, 100 - 23.3, 10000))
    sat_stat.to_calculate_data('resurs_p3', 'D:/Проекты/Оперативность/Бронницкое лесничество/TLE')
    print()
