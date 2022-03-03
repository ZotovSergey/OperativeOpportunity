import numpy as np
from datetime import datetime, timedelta
from calendar import isleap

# Номер дня 29 февраля в високосном году
NUMBER_OF_29_FEB = 60
# Количество дней в невисокосном году
DAYS_IN_NOT_LEAP_YEAR = 365
# Количество дней в високосном году
DAYS_IN_LEAP_YEAR = 366

DATETIPE64_START = np.datetime64('1970-01-01T00:00:00')


def to_determine_days_number_in_year(date):
    """
    @Описание:
        Метод определяет номер дня в году для некоторой даты.
    :param date: дата, для которой определяяется номер дня в году (datetime).
    :return: номер дня в году для date.
    """
    # Вычисление происходит для того года, в котором и date
    return (date - datetime(date.year, 1, 1, 0, 0, 0)).days + 1


def to_determine_days_number_in_not_leap_year(date):
    """
    @Описание:
        Метод определяет номер дня в невисокосном году для некоторой даты.
    :param date: дата, для которой определяяется номер дня в году невисокосном (datetime).
    :return: номер дня в невисокосном году для date.
    """
    # Определяется номер в году date.year
    if type(date) == type(DATETIPE64_START):
        timestamp = ((date - DATETIPE64_START) / np.timedelta64(1, 's'))
        date = datetime.utcfromtimestamp(timestamp)

    days_number_in_year = to_determine_days_number_in_year(date)
    # Если date.year - високосный, то номер дня корректируется при необходимости
    if (isleap(date.year)) and (days_number_in_year >= NUMBER_OF_29_FEB):
        days_number_in_year -= 1
    return days_number_in_year

def days_in_year(year):
    """
    @Описание:
        Метод определяет количество дней в году.
    :param year: год от рождества Христова (int)
    :return: оличество дней в году year (int)
    """
    if not isleap(year):
        return DAYS_IN_LEAP_YEAR
    else:
        days_in_current_year = DAYS_IN_NOT_LEAP_YEAR

def to_determine_date_by_days_number_in_year(days_number_in_year, year):
    """
    @Описание:
        Метод определяет дату по году и номеру дня в этом году.
    :param days_number_in_year: номер дня в году year.
    :param year: год определяемой даты.
    :return: дата по номеру дня days_number_in_year и году year (datetime).
    """
    return datetime(year, 1, 1, 0, 0, 0) + timedelta(days=days_number_in_year - 1)


def to_determine_date_by_days_number_in_not_leap_year(days_number_in_not_leap_year, year):
    """
    @Описание:
        Метод определяет дату по году и номеру дня в невисокосном году.
    :param days_number_in_not_leap_year: номер дня в невиссокосном году.
    :param year: год определяемой даты.
    :return: дата по номеру дня days_number_in_year и году year (datetime).
    """
    # Количество дней корректируется, если год year - високосный и если необходимо
    if (isleap(year)) and (days_number_in_not_leap_year >= NUMBER_OF_29_FEB):
        days_number_in_any_year = days_number_in_not_leap_year + 1
    else:
        days_number_in_any_year = days_number_in_not_leap_year
    return to_determine_date_by_days_number_in_year(days_number_in_any_year, year)


def seconds_in_year(time):
    """
    @Описание:
        Вычисляет количество секунд от заданного времени до начала года для этого времени.
    :param time: время, от которого количество секунд до начала года для этого времени
    :return: количество секунд от заданного времени до начала года для этого времени.
    """
    return (time - datetime(time.year, 1, 1, 0, 0, 0)).total_seconds()
