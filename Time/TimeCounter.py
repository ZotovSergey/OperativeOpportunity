import numpy as np

def to_get_time_counts_for_period(initial_time, final_time, step):
    """
    @Описание:
        Метод выводит все отсчеты времени между initial_time и final_time с интервалом step секунд, включая сами
            initial_time и final_time.
    :param initial_time: начальное время интервала (datetime)
    :param final_time: конечное время интервала (datetime)
    :param step: время между интервалами в секундах или долях секунд (float)
    :return: массив datetime64 значений времени между интервалами
    """
    return np.append(np.arange(initial_time, final_time, 1000 * step, dtype='datetime64[ms]'),
                     np.datetime64(final_time))